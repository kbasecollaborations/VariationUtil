import uuid
import shutil
import os
import subprocess
import logging
import time

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

logging.basicConfig(format='%(created)s %(levelname)s: %(message)s')
def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class VCFToVariation:
    def __init__(self, callback_url, scratch):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)
        self.wsc = Workspace(callback_url)
        self.wsc_appdev = Workspace("https://appdev.kbase.us/services/ws")

    def _parse_vcf_data(self, vcf_filepath):
        # TODO: move to pyvcf for information extraction
        contigs = []
        genotypes = []
        chromosomes = []
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
                if line.startswith('#'):
                    if line.startswith("#CHROM"):
                        # log("#CHROM encountered, exiting loop.")
                        genotypes = line.split()[9:]
                        log("Number Genotypes in vcf: {}".format(len(genotypes)))
                        # break

                    tokens = line.split("=")
                    if tokens[0].startswith('##contig'):
                        contigs.append(tokens[2][:-2])
                else:
                    tokens = line.split('\t')
                    if tokens[0] not in chromosomes:
                        chromosomes.append(tokens[0])

        return version, contigs, genotypes, chromosomes

    def _validate_vcf_to_sample(self, vcf_genotypes, sample_ids):
        check = True

        for geno in vcf_genotypes:
            if geno not in sample_ids:
                check = False

        return check

    def _chk_if_vcf_ids_in_assembly(self, vcf_chromosomes, assembly_chromosomes):
        check = True

        for chromo in vcf_chromosomes:
            if chromo not in assembly_chromosomes:
                check = False

        return check

    def validate_vcf(self, ctx, params):
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

        vcf_version, contigs, vcf_genotypes, vcf_chromosomes = self._parse_vcf_data(vcf_filepath)

        # vcftools (vcf-validator) supports VCF v4.0-4.2
        # https://github.com/vcftools/vcftools

        # EBIvariation/vcf-validator (vcf_validator_linux) supports VCF v4.14.3
        # https://github.com/EBIvariation/vcf-validator

        # vcftools is only to validate VCF v4.0

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

        validation_output_filename = os.path.join(validation_output_dir, 'vcf_validation.txt')

        try:
            if validator_output[0][:6] == '[info]':
                # validation by vcf_validator_linux
                vo = validator_output[1].split(' ')
                if os.path.exists(vo[6]):
                    shutil.move(vo[6], validation_output_filename)
                else:
                    raise Error('No output from vcf validator, check installation')
            else:
                if validator_output:
                    with open(validation_output_filename, 'w') as f:
                        for line in validator_output:
                            f.write(str(line))
                        f.close()
                else:
                    raise Error('No output from vcftools, check installation')
        except IndexError:
            raise IndexError('Validation output failed to provide any details! Please check VCF validator:' + str(validator_cmd[0]))

        if not os.path.exists(validation_output_filename):
            print('Validator did not generate log file!')
            raise Exception("Validator did not generate a log file.")

        log("Validator output filepath: {}".format(validation_output_filename))

        log("Return code from validator {}".format(p.returncode))

        return validation_output_filename, vcf_chromosomes, vcf_genotypes

    def _validate_assembly_ids(self, ctx, params, vcf_chromosomes):
        # All chromosome ids from the vcf should be in assembly
        # but not all assembly chromosome ids should be in vcf

        # TODO: validate really with DFU vs ref assembly

        subset = self.wsc_appdev.get_object_subset([{
            'included': ['/assembly_ref'],
            'ref': params['genome_ref']
        }])

        assembly_ref = subset[0]['data']['assembly_ref']

        assembly_chromosome_ids_call = self.wsc_appdev.get_object_subset([{
            'included': ['/contigs'],
            'ref': assembly_ref
        }])

        assembly_chromosomes = assembly_chromosome_ids_call[0]['data']['contigs'].keys()

        if not self._chk_if_vcf_ids_in_assembly(vcf_chromosomes, assembly_chromosomes):
            raise Error('VCF chromosome ids do not correspond to chromosome IDs from the assembly master list')

        return assembly_chromosomes

    def _validate_sample_ids(self, ctx, params, vcf_genotypes):
        # All samples within the VCF file need to be in sample attribute list
        # TODO: validate against real concurrently uploaded sample ids

        sample_ids_subset = self.wsc_appdev.get_object_subset([{
            'included': ['/instances'],
            'ref': params['sample_attribute_ref']
        }])

        sample_ids = sample_ids_subset[0]['data']['instances'].keys()

        if not self._validate_vcf_to_sample(vcf_genotypes, sample_ids):
            raise Error('VCF file genotype ids do not correspond to Sample IDs in the Sample Attribute master list')

        return sample_ids

    def import_vcf(self, ctx, params):
        file_valid_result, vcf_chr_ids, vcf_genos = self.validate_vcf(ctx, params)

        assembly_chr_ids = self._validate_assembly_ids(ctx, params, vcf_chr_ids)

        sample_ids = self._validate_sample_ids(ctx, params, vcf_genos)

        return 'cool'