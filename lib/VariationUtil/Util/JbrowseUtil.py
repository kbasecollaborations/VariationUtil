import gzip
import json
import logging
import os
import shutil
import subprocess
import uuid
from collections import Counter
import requests

from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class JbrowseUtil:
    def __init__(self, Config):
        callback_url = os.environ['SDK_CALLBACK_URL']
        ws_url = Config['ws_url']
        self.wsc = Workspace(ws_url)
        self.dfu = DataFileUtil(callback_url)
        self.gfu = GenomeFileUtil(callback_url)
        #service-wizard url
        self.sw_url = Config['sw_url']
        self.shock_url = Config['shock_url']
        scratch = Config['scratch']
        session = str(uuid.uuid4())
        self.session_dir = (os.path.join(scratch, session))
        os.mkdir(self.session_dir)
        pass

    def get_variation_service_url(self, sw_url):
        '''
        get the most recent VariationFileServ url from the service wizard.
        sw_url: service wizard url
        '''
        # TODO Fix the following dev thing to beta or release or future
        json_obj = {
            "method": "ServiceWizard.get_service_status",
            "id": "",
            "params": [{"module_name": "VariationFileServ", "version": "dev"}]
        }
        sw_resp = requests.post(url=sw_url, data=json.dumps(json_obj))
        vfs_resp = sw_resp.json()
        self.shock_url = self.shock_url.replace("https://", "")
        vfs_url = vfs_resp['result'][0]['url'] + "/jbrowse_query/" + self.shock_url + "/node"
        return vfs_url

    def _run_cmd(self, cmd):
        try:
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if stdout:
                logging.info("ret> ", process.returncode)
                logging.info("OK> output ", stdout)
            if stderr:
                logging.info("ret> ", process.returncode)
                logging.info("Error> error ", stderr.strip())

        except OSError as e:
            logging.info("OSError > ", e.errno)
            logging.info("OSError > ", e.strerror)
            logging.info("OSError > ", e.filename)

    def create_refseqs_data_from_assembly(self, assembly_ref):
        '''

        :param assembly_json:
        :return:
        '''
        refseqs_data = []
        # 1) Download assembly contig info and parse contig length information
        data = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': assembly_ref
        }])[0]['data']
        for key in data['contigs']:
            refseqs_data.append(
                {"end": data['contigs'][key]["length"],
                 "length": data['contigs'][key]["length"],
                 "name": data['contigs'][key]["contig_id"],
                 "seqChunkSize": 20000,
                 "start": 0
                 }
            )
        return refseqs_data


    def prepare_genome_features_track(self, genome_ref, vfs_url):
        """
        Builds track for genome features

        :param genome_ref:
        :return:
        """
        shock_handles = list()
        gff_track = ""

        # 1) Download gff using genomefileutil
        gff_file_info = self.gfu.genome_to_gff({'genome_ref': genome_ref})
        gff_file = gff_file_info["file_path"]

        # 2) sort gff
        outfile = gff_file + "_sorted"
        sorted_gff_cmd = " ".join(["sort -k1,1 -k4,4n",
                                  gff_file, ">", outfile])
        self._run_cmd(sorted_gff_cmd)

        # 3) compress gff
        zip_cmd = "bgzip " + outfile
        self._run_cmd(zip_cmd)

        # 4) index gff
        index_gff_cmd = "tabix -p gff " + gff_file + "_sorted.gz"
        self._run_cmd(index_gff_cmd)

        gff_gz_file_path = gff_file + "_sorted.gz"
        gff_index_file_path = gff_file + "_sorted.gz.tbi"

        # 5) Upload gff and gff index to shock
        if os.path.exists(gff_gz_file_path):
            gff_shock_ref = self.dfu.file_to_shock(
                {'file_path': gff_gz_file_path, 'make_handle': 1}
            )
        if os.path.exists(gff_index_file_path):
            gff_index_shock_ref = self.dfu.file_to_shock(
                {'file_path': gff_index_file_path, 'make_handle': 1}
            )

        # 6 Create gff track text that will be used for genome features track
        gff_track = '''
        {
            "label": "Genome Features",
            "key": "GenomeFeatures",
            "storeClass": "JBrowse/Store/SeqFeature/GFF3Tabix",
            "urlTemplate":"<vfs_url>/<gff_shock_ref>",
            "tbiUrlTemplate": "<vfs_url>/<gff_index_shock_ref>",
            "type": "JBrowse/View/Track/CanvasFeatures"
        }
        '''
        gff_track = gff_track.replace("<gff_shock_ref>",
                                      gff_shock_ref['handle']['id'])
        gff_track = gff_track.replace("<gff_index_shock_ref>",
                                      gff_index_shock_ref['handle']['id'])
        gff_track = gff_track.replace("<vfs_url>", vfs_url)
        gff_track_dict = json.loads(gff_track)

        # 7) Capture shock handles
        shock_handles.append(gff_shock_ref['handle'])
        shock_handles.append(gff_index_shock_ref['handle'])

        # 8) return shock handles and gff track info
        return {"shock_handle_list": shock_handles, "track_item": gff_track_dict}



    def prepare_snp_frequency_track(self, vcf_filepath, assembly_ref, binsize, vfs_url):
        """

        :param vcf_filepath:
        :param assembly_ref:
        :param binsize:
        :return:
        """
        BEDGRAPHTOBIGWIG="/kb/deployment/bin/bedGraphToBigWig"
        shock_handles = list()

        chr_length_dict = {}
        chr_length_data = ""
        chr_length_path = None
        counts = Counter()

        # 1) Download assembly contig info and parse contig length information
        data = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': assembly_ref
        }])[0]['data']

        contigs = data["contigs"]
        for contig in contigs:
            contig_data = data["contigs"][contig]
            chr_length_data += str(contig_data['contig_id']) + '\t' + str(contig_data['length']) + '\n'
            c_id = str(contig_data['contig_id'])
            c_length = str(contig_data['length'])
            chr_length_dict[c_id] = c_length

        # 2) Write contig lengths to a file (needed later)
        if chr_length_data is not None:
            chr_length_path = os.path.join(self.session_dir,
                                            "chr_length.txt")
            with open(chr_length_path, "w") as f:
                f.write(chr_length_data)

        # 3) Read and parse vcf file (must be bgzip compressed)
        #    Caclculate number of SNPs in each bin and write in bedgraph format
        reader = gzip.open(vcf_filepath, "rt")
        logging.info("Generating bedgraph file\n")
        for record in reader:
            if record[0] == "#":
                continue
            rs = record.split("\t")
            CHR, POS = rs[0], rs[1]
            bin_pos = int(int(POS) / binsize)
            bin_id = str(CHR) + "\t" + str(bin_pos)
            counts[bin_id] += 1
        bedgraph_file = os.path.join(self.session_dir, "vcf_bedgraph.txt")
        try:
            with open(bedgraph_file, "w") as fout:
                for j, k in counts.items():
                    chromosome, bin_num = j.split("\t")
                    bin_start = int(bin_num) * binsize
                    bin_end = bin_start + binsize
                    chr_length = chr_length_dict[chromosome]
                    if bin_end <= int(chr_length):
                        fout.write(chromosome + "\t" + str(bin_start) + "\t" + str(bin_end) + "\t" + str(k) + "\n")
                    else:
                        fout.write(chromosome + "\t" + str(bin_start) + "\t" + str(chr_length) + "\t" + str(k) + "\n")
        except IOError:
            logging.info("Unable to write " + bedgraph_file, + " file on disk.")

        # 4) Sort bedgraph file by chromosome id and co-ordinates
        sorted_bedgraph_file = bedgraph_file + "_sorted"
        sort_cmd = "sort -k1,1 -k2,2n " + bedgraph_file + "> " + sorted_bedgraph_file
        self._run_cmd(sort_cmd)

        # 5) Convert sorted bedgraph to bigwig format using utility bedgraphTOBigWig tool
        output_bigwig_file = bedgraph_file + "_bigwig.bw"
        cmd = BEDGRAPHTOBIGWIG + " " + sorted_bedgraph_file + " " + chr_length_path + " " + output_bigwig_file
        logging.info("Generating bigwig ..\n" + cmd + "\n")
        self._run_cmd(cmd)

        # 6) upload bigwig file to shock
        logging.info("Uploading Bigwig file to shock")
        if os.path.exists(output_bigwig_file):
            bigwig_shock_ref = self.dfu.file_to_shock(
                {'file_path': output_bigwig_file, 'make_handle': 1}
            )
        # 7) Append shock handle to genomic_indexes
        shock_handles.append(bigwig_shock_ref['handle'])

        # 8) Build snp frequency track
        output_bigwig_shock = bigwig_shock_ref['handle']['id']
        snp_frequency_track = '''
        {
            "label": "Variation Densityy", 
            "key": "Variation_density", 
            "storeClass": "JBrowse/Store/SeqFeature/BigWig", 
            "urlTemplate": "<vfs_url>/<bigwig_shock_id>", 
            "type": "JBrowse/View/Track/Wiggle/XYPlot"
        } 
        '''
        snp_frequency_track = snp_frequency_track.replace("<bigwig_shock_id>", output_bigwig_shock)
        snp_frequency_track = snp_frequency_track.replace("<vfs_url>", vfs_url)
        snp_frequency_track_dict = json.loads(snp_frequency_track)
        # 9) Return shock handles and track info
        return {"shock_handle_list": shock_handles, "track_item": snp_frequency_track_dict}


    def prepare_snp_track(self, vcf_shock_id, vcf_index_shock_id, vfs_url):
        """

        :param vcf_shock_id:
        :param vcf_index_shock_id:
        :return:
        """
        shock_handles = list()

        snp_track ='''
            {
                "label": "Variation", 
                "key": "Variation", 
                "storeClass": "JBrowse/Store/SeqFeature/VCFTabix", 
                "urlTemplate": "<vfs_url>/<vcf_shock_id>", 
                "tbiUrlTemplate": "<vfs_url>/<vcf_index_shock_id>", 
                "type": "JBrowse/View/Track/HTMLVariants"
            }
        '''
        snp_track = snp_track.replace("<vcf_shock_id>", vcf_shock_id)
        snp_track = snp_track.replace("<vcf_index_shock_id>", vcf_index_shock_id)
        snp_track = snp_track.replace("<vfs_url>", vfs_url)
        snp_track_dict = json.loads(snp_track)
        # shock handles should be empty list in return when built from shock ids
        return {"shock_handle_list": shock_handles, "track_item": snp_track_dict}

    def build_jbrowse_data_folder(self, jbrowse_path):
        shock_handles = list()
        data_folder_shock_ref = self.dfu.file_to_shock({'file_path': jbrowse_path,
                                            'pack': 'zip', 'make_handle': 1})
        shock_handles.append(data_folder_shock_ref['handle'])
        return {"shock_handle_list": shock_handles}

    def build_jbrowse(self, jbrowse_src, jbrowse_path, refseqs_data, genomic_indexes, tracklist_items):
        """

        :param jbrowse_src:
        :param jbrowse_path:
        :param genomic_indexes:
        :param tracklist_items:
        :return:
        """
        jbrowse_report = {}

        # 1) Copy the jbrowse source code to build report
        destination = shutil.copytree(jbrowse_src, jbrowse_path)

        # 2) Put tracklist.json in jbrowse data path
        tracklist_path = os.path.join(jbrowse_path, "data", "trackList.json")
        trackdata = {
            'formatVersion': 1,
            'tracks': tracklist_items
        }
        with open(tracklist_path, "w") as f:
            f.write(json.dumps(trackdata))

        # 3) Put refseq.json in jbrowse seq path
        refseqs_json_path = os.path.join(jbrowse_path, "data", "seq", "refSeqs.json")
        with open(refseqs_json_path, "w") as f:
            f.write(json.dumps(refseqs_data))

        #Build jbrowse data folder to support jbrowse widget in narrative
        res = self.build_jbrowse_data_folder(jbrowse_path)
        data_folder_index = res['shock_handle_list']
        genomic_indexes = genomic_indexes + data_folder_index

        # Build jbrowse report dict
        jbrowse_report["jbrowse_data_path"] = jbrowse_path
        jbrowse_report["genomic_indexes"] = genomic_indexes

        return jbrowse_report

    def prepare_jbrowse_report(self, jbrowse_params):
        """
        Build genomic indexes, prepare jbrowse report
        :param input_params:
        :return:
        """
        # Service wizard
        sw_url = self.sw_url
        # Variation file service url for serving jbrowse track files
        vfs_url = self.get_variation_service_url(sw_url)

        print(vfs_url)

        genomic_indexes = list()
        tracklist_items = list()
        refseqs_data = None

        # 1) Build refseqs_data
        #    This is used to build refseqs.json file for jbrowse
        #    Jbrowse report can not be built if assembly ref doesn't exist
        if 'assembly_ref' in jbrowse_params:
            assembly_ref = jbrowse_params['assembly_ref']
            refseqs_data = self.create_refseqs_data_from_assembly(assembly_ref)
        else:
            raise ValueError ("assembly ref not found")
            return

        # 2) Build genome features track
        if 'genome_ref' in jbrowse_params:
            genome_ref = jbrowse_params['genome_ref']
            output = self.prepare_genome_features_track(genome_ref, vfs_url)
            shock_handles, track_item = output["shock_handle_list"], output["track_item"]
            genomic_indexes = genomic_indexes + shock_handles
            tracklist_items.append(track_item)
        else:
            print ("Skipping genome features track")


        # 3) Build SNP frequency track
        cond1 = 'vcf_path' in jbrowse_params
        cond2 = 'assembly_ref' in jbrowse_params
        cond3 = 'binsize' in jbrowse_params
        if cond1 and cond2 and cond3:
            vcf_path = jbrowse_params['vcf_path']
            assembly_ref = jbrowse_params['assembly_ref']
            binsize = jbrowse_params["binsize"]
            output = self.prepare_snp_frequency_track(vcf_path, assembly_ref, binsize, vfs_url)
            shock_handles, track_item = output["shock_handle_list"], output["track_item"]
            if shock_handles:
                genomic_indexes = genomic_indexes + shock_handles
            tracklist_items.append(track_item)
        else:
            print ("Skipping SNP frequency track")

        # 4) Build SNP track
        cond1 = 'vcf_shock_id' in jbrowse_params
        cond2 = 'vcf_index_shock_id' in jbrowse_params
        if cond1 and cond2:
            vcf_shock_id = jbrowse_params['vcf_shock_id']
            vcf_index_shock_id = jbrowse_params['vcf_index_shock_id']
            output = self.prepare_snp_track(vcf_shock_id, vcf_index_shock_id, vfs_url)
            shock_handles, track_item = output["shock_handle_list"], output["track_item"]
            genomic_indexes = genomic_indexes + shock_handles
            tracklist_items.append(track_item)
        else:
            print ("Skipping SNP track")
        # 5) Build jbrowse directory with index.html
        # jbrowse directory later on gets uploaded as html report
        jbrowse_src = "/kb/module/deps/jbrowse"
        jbrowse_path = os.path.join(self.session_dir, "jbrowse")

        if tracklist_items:
             jbrowse_report = self.build_jbrowse(jbrowse_src,
                                                 jbrowse_path,
                                                 refseqs_data,
                                                 genomic_indexes,
                                                 tracklist_items)
        else:
            raise ValueError ("No tracks found")
            return

        return jbrowse_report

if __name__ == "__main__":
    Config = {
        "ws_url": "https://appdev.kbase.us/services/ws",
        "scratch": "/kb/module/work/tmp",
    }
    jbrowse_params = {
        "genome_ref": "41807/77/1",
        "assembly_ref": "41807/76/1",
        "binsize": 10000,
        "vcf_shock_id": "429099d1-30a4-45a1-a59e-bd98edd25a93",
        "vcf_index_shock_id": "aa475ffe-a3b6-4813-9532-d8c8de27282f"
    }
    jb = JbrowseUtil(Config)
    jbrows_report = jb.prepare_gbrowse_report (jbrowse_params)
    print (jbrows_report)
