import binascii
import gzip
import logging
import os
import re
import subprocess
import uuid
from itertools import islice


class VCFUtils:

    def __init__(self, Config):
        self.scratch = Config['scratch']

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

    def guess_effective_staging_path(self, staging_path):
        """
        Handles a bug in narrative staging file path.
        by prepending "/staging"

        :param filepath:
        :return: effective_staging_path
        """
        # Guess effective staging file path
        if (staging_path.startswith("/kb/module")):
            effective_staging_path = staging_path
        else:
            effective_staging_path = os.path.join("/staging",
                                                      staging_path)

        if (os.path.exists(effective_staging_path)):
            return (effective_staging_path)
        else:
            raise RuntimeError("File not found at "
                               + effective_staging_path)

    def is_gz_file(self, filepath):
        with open(filepath, 'rb') as test_f:
            return binascii.hexlify(test_f.read(2)) == b'1f8b'

    def looks_like_vcf_file(self, filepath):
        """
        Read up to first 100,000 lines of a VCF file
        to see if the file looks like a vcf file.
        This is a simple test to fail early if
        a line that starts with "#CHROM" is not present
        in the file

        :param filepath: input .gz or ascii format vcf file
        :return: "gzip", "text" or "unsupported"
        """
        filetype = None
        N = 100000
        if self.is_gz_file(filepath):
            with gzip.open(filepath, 'rt') as infile:
                lines_gen = islice(infile, N)
                for line in lines_gen:
                    if line.startswith("#CHROM"):
                        filetype = "gzip"

        else:
            with open(filepath, "r") as infile:
                lines_gen = islice(infile, N)
                try:
                    for line in lines_gen:
                        if line.startswith("#CHROM"):
                            filetype = "text"
                except:
                    filetype = "unsupported"

        if filetype is None:
            raise RuntimeError("Could not find valid VCF header line in "
                               + filepath)
        if filetype is "unsupported":
            raise RuntimeError("Unsupported file format in "
                               + filepath)
        return filetype

    def run_command(self, command):
        """
        This function runs a third party command line tool
        eg. bgzip etc.
        :param command: command to be run
        :return: success
        """
        logging.info("Running command " + command)
        cmdProcess = subprocess.Popen(command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT,
                                      shell=True)
        for line in cmdProcess.stdout:
            logging.info(line.decode("utf-8").rstrip())
            cmdProcess.wait()
            logging.info('return code: ' + str(cmdProcess.returncode))
            if cmdProcess.returncode != 0:
                raise ValueError('Error in running command with return code: '
                                 + command
                                 + str(cmdProcess.returncode) + '\n')
        logging.info("command " + command + " ran successfully")
        return "success"

    def bgzip_vcf_file(self, filepath, filetype, session_directory, filename):
        """
        Reads and compresses input vcf file (filepath) to a user specified
        destination (session_directory + filename)
        if filetype is gzip, it is decompressed and then recompressed using bgzip
        algorithm
        This also acts as a validation of vcf file because tabix indexing of vcf file
        will only happen for a valid vcf file
        :param filepath: user input vcf file (ascii or .gzip format)
        :return: path of compressed bgzip file
        """
        destination_path = os.path.join(session_directory,
                                        filename)
        if filetype == "gzip":
            command = " ".join(["zcat",
                                filepath,
                                "|bgzip -c >",
                                destination_path])
        if filetype == "text":
            command = " ".join(["cat",
                                filepath,
                                "|bgzip -c >",
                                destination_path])
        if filetype == "unsupported":
            raise RuntimeError("Unsupported format in vcf file")

        job_status = self.run_command(command)
        if job_status == "success":
            return destination_path
        else:
            raise RuntimeError("error in creating bgzipped staging file "
                               + filepath)

    def index_vcf_file(self, filepath):
        """
        Runs tabix command line tool to
        create .tbi index
        :param filepath: path of bgzipped vcf file
        :return: destination_path: path of index file
        """
        default_index_suffix = ".tbi"
        command = "tabix -p vcf " + filepath
        job_status = self.run_command(command)
        if (job_status == "success"):
            destination_path = filepath + default_index_suffix
            if (os.path.exists(destination_path)):
                return destination_path
            else:
                raise RuntimeError("Index created but can not find index file "
                                   + destination_path)
        else:
            raise RuntimeError("Problem creating index file for "
                               + filepath)

    def check_and_fix_headers(self, filepath):
        """
        Some software put the entire file path instead of sample
        name in the VCF column header. This function checks the
        column headers and removes everything before / if it looks
        like a path.
        tabix -H only dumps the header from the vcf file
        tabix -r optimizes the reheading process
        :param filepath:
        :return:
        """
        changed_headers = 0
        dest_path = filepath + "_header.txt"
        command = " ".join(["tabix -H",
                            filepath,
                            " > ",
                            dest_path])

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
            command = " ".join(["tabix -r",
                                new_header_path,
                                filepath,
                                ">",
                                reheader_vcf_path])

            job_status = self.run_command(command)
            if job_status == "success":
                if (os.path.exists(reheader_vcf_path)):
                    reheader_index_file_path = self.index_vcf_file(reheader_vcf_path)
                    return (reheader_vcf_path, reheader_index_file_path)
        else:
            return



    def validate_compress_and_index_vcf(self, params):
        """
        Parses VCF file, validates, compresses, indexes
        and returns vcf file and index file path
        :return:
        """

        staging_path = params['vcf_staging_file_path']

        # 1) Generate a session id to keep all intermediate files
        # in one common folder
        session_id = str(uuid.uuid4())
        session_directory = os.path.join(self.scratch, session_id)
        self._mkdir_p(session_directory)

        # 2) Guess staging file path. Takes care of the path issue in staging
        effective_staging_path = self.guess_effective_staging_path(staging_path)

        # 3) Quickly guess filetype (Just read first 100,000 lines) and fail if
        # file is not valid vcf. This is needed for really large files with millions
        # of rows.
        logging.info("Guessing VCF file type and compression")
        filetype = self.looks_like_vcf_file(effective_staging_path)
        logging.info("Input file type is " + filetype)

        # 4) Create bgzip compressed VCF variation
        # bgzip compressed variation can be used with large number of
        # tools that work with vcf files. Default index is .tbi
        logging.info("Compressing VCF file using bgzip")
        bgzip_filename = "variation.vcf.gz"
        bgzip_filepath = self.bgzip_vcf_file(effective_staging_path,
                                             filetype,
                                             session_directory,
                                             bgzip_filename)
        logging.info("compressed vcf is in " + bgzip_filepath)

        # 5) Create vcf index using tabix
        logging.info("Creating index")
        bgzip_index_filepath = self.index_vcf_file(bgzip_filepath)
        logging.info("Created index " + bgzip_index_filepath)

        # 6) Check if the column header in the VCF file has characters like "/" fix them
        logging.info("Checking if headers in VCF need to be fixed")
        reheader_result = self.check_and_fix_headers(bgzip_filepath)
        if not reheader_result:
            final_vcf, final_index = bgzip_filepath, bgzip_index_filepath
        else:
            final_vcf, final_index = reheader_result
        return (final_vcf, final_index)


if __name__ == '__main__':
    session_id = uuid.uuid4()
    params = {
        "vcf_staging_file_path": '/kb/module/test/sample_data/small_poplar/small_poplar.vcf.gz',
    }
    config = {
        "scratch": "/kb/module/work"
    }
    vv = VCFUtils(params, config)
    final_vcf, final_index = vv.sanitize_vcf()
    vcf_info = vv.parse_vcf_data(final_vcf)
    filename = "/kb/module/work/vcf_info.json"
    with open(filename, "w") as f:
        f.write(vcf_info)
    print(final_vcf)
    print(final_index)
