import uuid
import shutil
import os
import subprocess
import logging
import time

from installed_clients.DataFileUtilClient import DataFileUtil

logging.basicConfig(format='%(created)s %(levelname)s: %(message)s')
def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class VCFToVariation:
    def __init__(self, callback_url, scratch):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)

    def _parse_vcf_data(self, vcf_filepath):
        # TODO: move to pyvcf for information extraction
        contigs = []
        genotypes = []
        version = ''
        with(gzip.open if vcf_filepath.endswith('.gz') else open)(vcf_filepath, 'rt') as vcf:
            line = vcf.readline()
            tokens = line.split('=')

            if not(tokens[0].startswith('##fileformat')):
                log("Invalid VCF.  ##fileformat line in meta is improperly formatted.")
                raise ValueError("Invalid VCF.  ##fileformat line in meta is improperly formatted.")
            version = float(tokens[1][-4:].rstrip())
            log("VCF version: {}".format(version))
            for line in vcf:
                if line.startswith("#CHROM"):
                    log("#CHROM encountered, exiting loop.")
                    genotypes = line.split()[9:]
                    log("Number Genotypes in vcf: {}".format(len(genotypes)))
                    break
                tokens = line.split("=")

                if tokens[0].startswith('##contig'):
                    contigs.append(tokens[2][:-2])
        return version, contigs, genotypes

    def _validate_vcf(self, ctx, params):
        if 'genome_ref' not in params:
            raise ValueError('Genome reference not in input parameters: \n\n'+params)
        if 'vcf_staging_file_path' not in params:
            raise ValueError('VCF staging file path not in input parameters: \n\n' + params)
        if 'variation_object_name' not in params or params['variation_object_name'] is None or params['variation_object_name'] == '':
            var_obj_name = 'Variation_'+str(uuid.uuid4())

        # TODO: shutil.move from staging location to scratch
        vcf_filepath = params['vcf_staging_file_path']

        validation_output_dir = os.path.join(self.scratch, 'validation_' + str(uuid.uuid4()))
        os.mkdir(validation_output_dir)

        vcf_version, contigs, genotypes = self._parse_vcf_data(vcf_filepath)

        if vcf_version >= 4.1:
            print("Using vcf_validator_linux...")
            validator_cmd = ["vcf_validator_linux"]
            validator_cmd.append("-i")
            validator_cmd.append(vcf_filepath)
            validator_cmd.append("-o")
            validator_cmd.append(validation_output_dir)
        else:
            print("Using vcftools to validate...")
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_filepath)
            print("VCF version below 4.1.  No validation logging.")

        print("Validator command: {}".format(validator_cmd))
        p = subprocess.Popen(validator_cmd,
                             cwd=self.scratch,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)

        validator_output = []
        while True:
            line = p.stdout.readline()
            if not line:
                break
            validator_output.append(line)

        p.wait()

        validation_output_filename = [f for f in os.listdir(validation_output_dir) if f.endswith('.txt')][0]
        validation_output_filepath = os.path.join(validation_output_dir, validation_output_filename)

        if not validation_output_filename:
            print('Validator did not generate log file!')
            raise Exception("Validator did not generate a log file.")

        log("Validator output filepath: {}".format(validation_output_filepath))

        log("Return code from validator {}".format(p.returncode))

        return validation_output_filepath, p.returncode

    def import_vcf(self, ctx, params):
        self._validate_vcf(ctx, params)

        return 'cool'