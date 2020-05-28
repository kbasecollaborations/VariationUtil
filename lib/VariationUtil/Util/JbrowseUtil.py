import logging
from installed_clients.WorkspaceClient import Workspace
import json
import gzip
import uuid
import os
from collections import Counter
import subprocess
import shutil
import logging
from installed_clients.GenomeFileUtilClient import GenomeFileUtil

class JbrowseUtil:
    def __init__(self):
        pass

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

    def get_gff_track(self, gff_shock_ref, gff_index_shock_ref):

        gff_track = '''
        {
            "label": "Genome Features",
            "key": "GenomeFeatures",
            "storeClass": "JBrowse/Store/SeqFeature/GFF3Tabix",
            "urlTemplate":"https://appdev.kbase.us/dynserv/682063b283a644bbcb27ca7a49919b8093608d05.VariationFileServ/shock/<gff_shock_ref>",
            "tbiUrlTemplate": "https://appdev.kbase.us/dynserv/682063b283a644bbcb27ca7a49919b8093608d05.VariationFileServ/shock/<gff_index_shock_ref>",
            "type": "JBrowse/View/Track/CanvasFeatures"
        }
        '''
        gff_track = gff_track.replace("<gff_shock_ref>", gff_shock_ref)
        gff_track = gff_track.replace ("<gff_index_shock_ref>", gff_index_shock_ref )

        return gff_track


    def prepare_gff(self, genome_ref):
        gfu = GenomeFileUtil(self.callback_url)
        gff_file_info = gfu.genome_to_gff({'genome_ref': genome_ref})
        gff_file =  gff_file_info["file_path"]
        #gff_file = "/kb/module/work/Populus_trichocarpa.gff"

        sorted_gff_cmd = "sort -k1,1 -k4,4n " + gff_file + " > " +  gff_file + "_sorted"
        self._run_cmd(sorted_gff_cmd)

        zip_cmd = "bgzip "  + gff_file + "_sorted"
        self._run_cmd(zip_cmd)

        index_gff_cmd = "tabix -p gff "   + gff_file + "_sorted.gz"
        self._run_cmd(index_gff_cmd)

        gff_gz_file_path = gff_file + "_sorted.gz"
        gff_index_file_path = gff_file + "_sorted.gz.tbi"

        if os.path.exists(gff_gz_file_path):
            gff_shock_ref = self.dfu.file_to_shock(
                {'file_path': gff_gz_file_path, 'make_handle': 1}
        )
        if os.path.exists(gff_index_file_path):
            gff_index_shock_ref = self.dfu.file_to_shock(
                {'file_path': gff_index_file_path, 'make_handle': 1}
            )

        return {"gff_shock_ref":   gff_shock_ref , "gff_index_shock_ref": gff_index_shock_ref}

    def create_refseqs_json_from_assembly(self ):
        '''

        :param assembly_json:
        :return:
        '''
        refseqs_json_data = []
        with open (self.assembly_json_file) as json_data:
            data = json.load(json_data)
        for key in data['contigs']:
            refseqs_json_data.append(
                {"end": data['contigs'][key]["length"],
                 "length": data['contigs'][key]["length"],
                 "name": data['contigs'][key]["contig_id"],
                 "seqChunkSize": 20000,
                 "start": 0
                 }
            )
        self.refseqs_json_path = os.path.join(self.session_dir, "refSeqs.json")
        with open (self.refseqs_json_path, "w") as f:
            f.write(json.dumps(refseqs_json_data))

        if os.path.exists(self.refseqs_json_path):
            return self.refseqs_json_path
        else:
            raise ValueError("File not found: " + self.refseqs_json_path)

    def create_chr_length_file(self):
        self.chr_length_dict = {}
        chr_length_data = ''
        with open(self.assembly_json_file) as json_data:
            data = json.load(json_data)

        contigs = data["contigs"]
        for contig in contigs:
            contig_data = data["contigs"][contig]
            chr_length_data += str(contig_data['contig_id']) + '\t' + str(contig_data['length']) + '\n'
            c_id = str(contig_data['contig_id'])
            c_length = str(contig_data['length'])
            self.chr_length_dict[c_id] = c_length
        self.chr_length_path = os.path.join(self.session_dir,
                                            "chr_length.txt")
        with open(self.chr_length_path, "w") as f:
            f.write(chr_length_data)

        if os.path.exists(self.chr_length_path):
            return self.chr_length_path
        else:
            raise ValueError("File not found: " + self.chr_length_path)


    def create_bedgraph_from_vcf(self):
        vcf_filepath = self.vcf_filepath
        reader = gzip.open(vcf_filepath, "rt")
        counts = Counter()
        logging.info("Generating bedgraph file\n")
        for record in reader:
            if record[0]=="#":
                 continue
            rs = record.split ("\t")
            CHR,POS= rs[0], rs[1]
            bin_pos = int(int(POS) / self.binsize)
            bin_id = str(CHR) + "\t" + str(bin_pos)
            counts[bin_id] += 1
        bedgraph_file = os.path.join(self.session_dir , "vcf_bedgraph.txt")
        try:
                with open(bedgraph_file, "w") as fout:
                    for j, k in counts.items():
                        #logging.info (str(j) + "\t" + str(k))
                        chromosome, bin_num = j.split("\t")
                        bin_start = int(bin_num) * self.binsize
                        bin_end = bin_start + self.binsize


                        chr_length = self.chr_length_dict[chromosome]
                        if bin_end <= int(chr_length):
                            fout.write(chromosome + "\t" + str(bin_start) + "\t" + str(bin_end) + "\t" + str(k) + "\n")
                        else:
                            fout.write(chromosome + "\t" + str(bin_start) + "\t" + str(chr_length) + "\t" + str(k) + "\n")

        except IOError:
                logging.info("Unable to write " + bedgraph_file, + " file on disk.")

        self.output_bigwig_file = bedgraph_file + "_bigwig.bw"
        sorted_bedgraph_file = bedgraph_file + "_sorted"
        sort_cmd = "sort -k1,1 -k2,2n " + bedgraph_file + "> " + sorted_bedgraph_file
        self._run_cmd(sort_cmd)
        cmd = self.bedGraphToBigWig  + " " + sorted_bedgraph_file + " " + self.chr_length_path + " " + self.output_bigwig_file
        logging.info("Generating bigwig ..\n" + cmd + "\n")
        self._run_cmd(cmd)


        if os.path.exists(self.output_bigwig_file):
            return self.output_bigwig_file
        else:
            logging.info ("error in generating: " + output_bigwig_file)





    def prepare_jbrowse_report(self, input_params):

        genomic_indexes = list()
        # input_params={'ws_url': 'https://appdev.kbase.us/services/ws', 'assembly_ref': '1745/511/24', 'scratch': '/kb/module/work/tmp', 'vcf_filepath': '/kb/module/work/tmp/791434bb-7ed5-435e-8a75-ad757dc1b30a/variation.vcf.gz', 'binsize': 10000, 'bedGraphToBigWig': '/kb/deployment/bin/bedGraphToBigWig', 'vcf_shock_id': 'fb361afb-c2dd-4ae3-9929-031344287270', 'vcf_index_shock_id': 'caad1a17-aeb4-4050-bcae-2d9eaa7d5cc1'}

       # self.assembly_json_file = input_params['assembly_json_file']
        self.assembly_ref = input_params['assembly_ref']
        ws_url = input_params['ws_url']
        self.callback_url = input_params['callback_url']
        self.wsc = Workspace(ws_url)
        self.dfu = input_params['dfu']
        if 'genome_ref' in input_params:
            genome_ref = input_params['genome_ref']
            gff_info= self.prepare_gff(genome_ref)
            gff_shock_ref_handle = gff_info['gff_shock_ref']['handle']
            gff_index_shock_ref_handle = gff_info['gff_index_shock_ref']['handle']
            genomic_indexes.append(gff_shock_ref_handle)
            genomic_indexes.append(gff_index_shock_ref_handle)
            self.gff_shock = gff_shock_ref_handle['id']
            self.gff_index_shock = gff_index_shock_ref_handle['id']



        data = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': self.assembly_ref
        }])[0]['data']
        self.scratch = input_params['scratch']
        session = str(uuid.uuid4())
        self.session_dir = (os.path.join(self.scratch, session))
        os.mkdir (self.session_dir)
        self.assembly_json_file = self.session_dir + "/assembly_json_file"
        with open (self.assembly_json_file, "w") as f:
            f.write(json.dumps(data))

        self.vcf_filepath = input_params['vcf_filepath']
        self.binsize = input_params['binsize']
        self.bedGraphToBigWig = input_params['bedGraphToBigWig']
        self.vcf_shock_id = input_params['vcf_shock_id']
        self.vcf_index_shock_id = input_params['vcf_index_shock_id']

        output = {}
        output['refseqs_json_path'] = self.create_refseqs_json_from_assembly()
        self.create_chr_length_file()
        self.create_bedgraph_from_vcf()
        output['output_bigwig_file'] = self.output_bigwig_file

        logging.info("Uploading Bigwig file to shock")
        if os.path.exists(self.output_bigwig_file):
            bigwig_shock_ref = self.dfu.file_to_shock(
                {'file_path': self.output_bigwig_file, 'make_handle': 1}
            )
        self.output_bigwig_shock = bigwig_shock_ref['handle']['id']
        logging.info (self.output_bigwig_shock)
        genomic_indexes.append(bigwig_shock_ref['handle'])

        jbrowse_src = "/kb/module/deps/jbrowse"
        jbrowse_dest = os.path.join(self.session_dir, "jbrowse")
        destination = shutil.copytree(jbrowse_src, jbrowse_dest)
        logging.info("After copying file:")
        logging.info(os.listdir(destination))

        jbrowse_path = os.path.join(self.session_dir, "jbrowse")
        jbrowse_seq_path = os.path.join(self.session_dir, "jbrowse", "data", "seq")
        jbrowse_data_path = os.path.join(self.session_dir, "jbrowse", "data")
        dest = shutil.copy(self.refseqs_json_path, jbrowse_seq_path  )
        logging.info ("dest is " + dest)
        logging.info("After copying refseqs json seq path:")
        logging.info(os.listdir(jbrowse_seq_path))

        logging.info("After copying refseqs json data path:")
        logging.info(os.listdir(jbrowse_data_path))
        logging.info("Jbrowse seq path")
        logging.info(jbrowse_seq_path)
        with open (self.refseqs_json_path) as f:
            data = f.read()
        logging.info ("Refsesa data")
        logging.info (data)



        #dest = shutil.move(self.output_bigwig_file, jbrowse_data_path + "/vcf.bw")
        #logging.info (dest)
        tracklist_path = os.path.join(jbrowse_data_path, "trackList.json")
        with open (tracklist_path, "r") as f:
            data = f.read()

        data=data.replace("<vcf_shock_id>", self.vcf_shock_id)
        data=data.replace("<vcf_index_shock_id>", self.vcf_index_shock_id)
        data = data.replace("<output_bigwig_shock>", self.output_bigwig_shock)
        data_j = json.loads(data)
        tracks = data_j['tracks']
 
        if "genome_ref" in input_params:
            gff_track = self.get_gff_track(self.gff_shock, self.gff_index_shock)
            gff_track_obj = json.loads(gff_track)
            tracks.append(gff_track_obj)

        trackdata = {
        'formatVersion':1,
        'tracks':tracks
        }
        with open (tracklist_path, "w") as f:
            f.write(json.dumps(trackdata))

        logging.info (jbrowse_dest)
        jbrowse_report = {}
        jbrowse_report["jbrowse_data_path"] = jbrowse_dest
        jbrowse_report["genomic_indexes"] = genomic_indexes
        return jbrowse_report


#        chr_length_data += str(data['contigs'][key]["contig_id"]) + "\t" + str(data['contigs'][key]["length"]) + "\n"




if __name__ == "__main__":
    #---Assembly file as json.

    #create refseqs.json from assembly file
    input_params = {
         'ws_url': 'https://appdev.kbase.us/services/ws', 'assembly_ref': '1745/511/24',
         'scratch': '/kb/module/work/tmp',
         'vcf_filepath': '/kb/module/work/tmp/387ea94b-7789-41c4-b932-8a76d8c9d782/variation.vcf.gz', 'binsize': 10000,
         'bedGraphToBigWig': '/kb/deployment/bin/bedGraphToBigWig',
         'vcf_shock_id': 'fb361afb-c2dd-4ae3-9929-031344287270',
         'vcf_index_shock_id': 'caad1a17-aeb4-4050-bcae-2d9eaa7d5cc1'}

    vu = JbrowseUtil()
    refseqs_json_path  =  vu.prepare_jbrowse_report(input_params)
    logging.info (refseqs_json_path)
    #create chr length from assembly file
    #---Path to vcf file
    #create bedgraph from vcf file
    #create bigwig from bedgraph and length file
    #get_tracklist_variation from shock urls
    #get_tracklist_density from bigwig file path
    #Returns a directory of file path outside of jbrowse
