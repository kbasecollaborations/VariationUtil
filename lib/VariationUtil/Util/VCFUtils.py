
import gzip
import binascii
from itertools import islice
import subprocess
import uuid
import os
import re
import logging
from collections import Counter

#logging.basicConfig(format='%(created)s %(levelname)s: %(message)s')
logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                    level=logging.INFO)


#logging.basicConfig(
#    format='%(asctime)s %(levelname)-8s %(message)s',
#    level=logging.INFO,
#    datefmt='%Y-%m-%d %H:%M:%S')




class VCFUtils:
    def __init__(self, params, config):
        self.params = params
        self.scratch = config['scratch']

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            raise RuntimeError(f"Specified path is not valid:{path}")
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def guess_effective_staging_path (self, filepath):
        # Guess effective staging file path
        staging_file_path = self.params["vcf_staging_file_path"]
        if (staging_file_path.startswith("/kb/module")):
            effective_staging_filepath = staging_file_path
        else:
            effective_staging_filepath = os.path.join("/staging", staging_file_path)

        if (os.path.exists(effective_staging_filepath)):
            return (effective_staging_filepath)
        else:
            raise RuntimeError ("File not found at " + effective_staging_filepath)

    def is_gz_file(self, filepath):
        with open(filepath, 'rb') as test_f:
            return binascii.hexlify(test_f.read(2)) == b'1f8b'

    def is_ascii_file(self, filepath):
        #TODO: Put proper test for ascii file
        return True

    def looks_like_vcf_file (self, filepath):
        """
        Read up to first 10,000 lines of a uncompressed VCF file
        to see if the file looks like a vcf file

        :param filepath:
        :return:
        """
        filetype = None
        N=100000
        if self.is_gz_file(filepath):
            with gzip.open(filepath, 'rt') as infile:
                lines_gen = islice(infile, N)
                for line in lines_gen:
                    if line.startswith("#CHROM"):
                        filetype = "gzip"

        else:
            with open (filepath, "r") as infile:
                lines_gen = islice(infile, N)
                try:
                    for line in lines_gen:
                        if line.startswith("#CHROM"):
                            filetype = "text"
                except:
                    filetype = "unsupported"

        if filetype is None:
            raise RuntimeError ("Could not find valid VCF header line in input file " + filepath)
        if filetype is "unsupported":
            raise RuntimeError("Unsupported file format in " + filepath)
        return filetype


    def run_command(self, command):
        logging.info ("Running command " + command)
        cmdProcess = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
        for line in cmdProcess.stdout:
            logging.info(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            logging.info('return code: ' + str(cmdProcess.returncode))
            if cmdProcess.returncode != 0:
                raise ValueError('Error in running command with return code: ' + command +
                                 str(cmdProcess.returncode) + '\n')

        logging.info("command " + command + " ran successfully")
        return "success"

    def stage_vcf_file (self, filepath, filetype, session_directory, filename):
        """

        :param filepath:
        :return:
        """
        destination_path = os.path.join (session_directory, filename)
        if filetype == "gzip":
            command = "zcat " +  filepath  + " |bgzip -c > " +  destination_path
        if filetype == "text":
            command = "cat " +  filepath  + "|bgzip -c > " +  destination_path

        job_status = self.run_command(command)
        if job_status == "success":
            return destination_path
        else:
            raise RuntimeError ("error in creating bgzipped staging file " + filepath)

    def index_vcf_file (self, filepath):
        default_index_suffix = ".tbi"
        command = "tabix -p vcf " + filepath
        job_status = self.run_command(command)
        if (job_status == "success"):
            destination_path = filepath + default_index_suffix
            if (os.path.exists ( destination_path )):
                return destination_path
            else:
                raise RuntimeError ("Index created but can not file index file " + destination_path)
        else:
            raise RuntimeError ("Problem creating index file for " + filepath)


    def check_and_fix_headers (self, filepath):
        changed_headers = 0
        dest_path = filepath +  "_header.txt"
        command = "tabix -H " + filepath + " > " + dest_path
        job_status = self.run_command(command)
        if job_status == "success":
            new_header_data = ""
            with open(dest_path, "r") as f:
                for line in f:
                    if (line.startswith("#CHROM")):
                        header_line = line
                        header_line_elements = header_line.split("\t")
                        new_header = list()
                        # Take care of samples with "/"
                        pattern1 = re.compile(r".*/")
                        for sample in header_line_elements:
                            # can use multiple_patterns here one in each line
                            # sample_new = pattern2.sub("", sample_new)
                            sample_new = pattern1.sub("", sample)
                            if sample_new != sample:
                                changed_headers += 1
                                logging.info("Replacing old " + sample + " to " + sample_new)
                                new_header.append(sample_new)
                            else:
                                new_header.append(sample)
                        new_header_data += "\t".join(new_header)
                    else:
                        new_header_data += line
        if (changed_headers > 0):
            new_header_path = filepath + "_new_header.txt"
            with open(new_header_path, "w") as wf:
                wf.write(new_header_data)

            reheader_vcf_path = filepath.replace(".vcf.gz", "") + "_reheader.vcf.gz"
            command = "tabix -r " + new_header_path + " " + filepath +  " > " +  reheader_vcf_path
            job_status = self.run_command(command)
            if job_status == "success":
                if (os.path.exists(reheader_vcf_path)):
                    reheader_index_file_path = self.index_vcf_file (reheader_vcf_path)
                    return (reheader_vcf_path, reheader_index_file_path)
        else:
            return

    def _parse_header(self, record, category):
        '''
        parses vcf header which looks like the following

        ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
        ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
        ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
        ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
        ##FILTER=<ID=q10,Description="Quality below 10">
        ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
        '''
        #TODO: Improve code
        returninfo = {"Category": category}
        info = re.sub(".*<", "", record)
        info = info.replace(">", "")
        infolist = info.replace('"', '').rstrip().split(",")
        for fields in infolist:
            data = fields.split("=")
            key = data.pop(0)
            val = "=".join (data)
            val = val.replace("\"", "")
            returninfo[key] = val
        return returninfo

    def parse_vcf_data(self, vcf_filepath):
        # TODO: Take care of vcf local file path
        # TODO: Attach chromosome length in main code
        #
        #vcf_filepath = params['vcf_local_file_path']
        reader = gzip.open(vcf_filepath, "rt")
        # version = get_version(vcf_filepath)
        # genotypes = get_genotypes(vcf_filepath)
        version = ""
        genotypes = ""
        #
        #TODO: Take care of length
        # 'length': int(record.affected_end - record.affected_start),
        #  Remove passvariants
        counter = 0
        chromosomes = []
        contigs = {}
        header = list()
        totalvars = 0

        for record in reader:
            if record[0] == "#":
                if record.startswith("##fileformat"):
                    version = record.replace("##fileformat=", "").rstrip()
                if record.startswith("##INFO=<"):
                    info = self._parse_header(record, "INFO")
                    header.append(info)
                if record.startswith("##FORMAT=<"):
                    info = self._parse_header(record, "FORMAT")
                    header.append(info)
                if record.startswith("##FILTER=<"):
                    info = self._parse_header(record, "FILTER")
                    header.append(info)

                if (record.startswith("#CHROM")):
                    # This is the chrome line
                    record = record.rstrip()
                    values = record.split("\t")
                    genotypes = values[9:]
                continue
            counter = counter + 1
            values = record.split("\t")
            totalvars += 1
            CHROM, POS = values[0], values[1]
            if CHROM not in chromosomes:
                chromosomes.append(CHROM)
                contigs[CHROM] = {
                    'contig_id': CHROM,
                    'totalvariants': 1,
                    'passvariants': 0,
                }
            else:
                contigs[CHROM]['totalvariants'] += 1
        vcf_info = {
            'version': version,
            'contigs': contigs,
            'total_variants': totalvars,
            'genotype_ids': genotypes,
            'chromosome_ids': chromosomes,
            'file_ref': vcf_filepath,
            'header': header
        }

        return vcf_info

    def sanitize_vcf(self):

        #Generate a session id to keep all files in one common folder
        session_id = str(uuid.uuid4())
        session_directory = os.path.join (self.scratch, session_id)
        self._mkdir_p(session_directory)

        #Guess staging file path. Takes care of the path issue in staging
        effective_staging_path = self.guess_effective_staging_path (self.params["vcf_staging_file_path"])

        # Quicly guess filetype (Just read first 100,000 lines) and fail if
        # file is not valid vcf. This is needed for really large files with millions
        # of rows.
        logging.info ("Guessing VCF file type and compression")
        filetype = self.looks_like_vcf_file(effective_staging_path)
        logging.info ("Input file type is " + filetype)

        # Create bgzip compressed VCF variation and index
        # bgzip compressed variation can be used with large number of
        # tools that work with vcf files. Default index is .tbi
        logging.info ("Compressing VCF file using bgzip")
        bgzip_filename = "variation.vcf.gz"
        bgzip_filepath = self.stage_vcf_file (effective_staging_path, filetype,
                                              session_directory, bgzip_filename  )
        logging.info ("Created compressed file" + bgzip_filepath)

        logging.info ("Creating index for compressed vcf")
        bgzip_index_filepath = self.index_vcf_file (bgzip_filepath)
        logging.info ("Created index " + bgzip_index_filepath)

        # Check if the column header in the VCF file has characters like / fix them
        logging.info ("Checking if headers in VCF need to be fixed")

        reheader_result = self.check_and_fix_headers(bgzip_filepath)
        if not reheader_result:
            final_vcf, final_index = bgzip_filepath, bgzip_index_filepath
        else:
            final_vcf, final_index = result
        return(final_vcf, final_index)


if __name__ == '__main__':
    session_id = uuid.uuid4()
    params={
        "vcf_staging_file_path" : '/kb/module/work/xu.vcf.gz',
    }
    config ={
        "scratch":"/kb/module/work"
    }
    vv = VCFUtils (params, config)
    final_vcf, final_index = vv.sanitize_vcf()
    vcf_info = vv.parse_vcf_data(final_vcf)
    filename = "/kb/module/work/vcf_info.json"
    with open(filename, "w") as f:
        f.write(vcf_info)
    print (final_vcf)
    print (final_index)

