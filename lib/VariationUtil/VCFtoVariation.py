import uuid
import shutil
import os
import subprocess
import logging
import time

os.system("find / -name 'AT-var-import.vcf'")
exit()

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

        # TODO: .move from staging location to scratch

        try:
            if ctx['test_env']:
                vcf_filepath = params['vcf_staging_file_path']
            else:
                vcf_filepath = self._save_staging_to_local(ctx, params)
        except KeyError:
            vcf_filepath = self._save_staging_to_local(ctx, params)

        validation_output_dir = os.path.join(self.scratch, 'validation_' + str(uuid.uuid4()))
        os.mkdir(validation_output_dir)

        vcf_version, vcf_contigs, vcf_genotypes, vcf_chromosomes = self._parse_vcf_data(vcf_filepath)

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
            validator_cmd.append("-l")
            validator_cmd.append('error')
            print("VCF version "+str(vcf_version)+".")
        else:
            print("Using vcftools to validate...")
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_filepath)
            print("VCF version 4.0.")

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
            validator_output.append(line.decode("utf-8"))

        p.wait()

        validation_output_filename = os.path.join(validation_output_dir, 'vcf_validation.txt')
        file_output_chk = []

        try:
            if validator_output[0][:6] == '[info]':
                # validation by vcf_validator_linux
                validation_output_filename = validator_output[1].split(' ')[6].strip('\n')
                vo = validator_output[2].split(' ')
                file_output_chk = ''.join(vo[9:]).strip('\n')

                if not os.path.exists(validation_output_filename):
                    raise ValueError(validation_output_filename+' does not exist!')

                if not file_output_chk == 'isvalid':
                    print('\n'.join(validator_output))
                    raise ValueError('\n'.join(validator_output))

                #TODO: more detailed validation parsing for vcf_validator_linux
            else:
                if validator_output:
                    with open(validation_output_filename, 'w') as f:
                        for line in validator_output:
                            f.write(str(line))
                        f.close()
                    print('\n'.join(validator_output))
                    raise ValueError('\n'.join(validator_output))
                else:
                    with open(validation_output_filename, 'w') as f:
                        f.write("vcftools used to validate vcf file:\n"+vcf_filepath+"\n\File is validate as of vcf spec v4.0")
                        f.close()

                # TODO: more detailed validation parsing for vcftools
        except IndexError:
            # if vcf file < v4.1, and valid it will produce index error on line 132
            if validator_output:
                with open(validation_output_filename, 'w') as f:
                    for line in validator_output:
                        f.write(str(line))
                    f.close()
                print('\n'.join(validator_output))
                raise ValueError('\n'.join(validator_output))
            else:
                with open(validation_output_filename, 'w') as f:
                    f.write("vcftools used to validate vcf file:\n" + vcf_filepath + "\n\File is validate as of vcf spec v4.0")
                    f.close()

        if not os.path.exists(validation_output_filename):
            print('Validator did not generate log file!')
            raise SystemError("Validator did not generate a log file.")

        log("Validator output filepath: {}".format(validation_output_filename))

        log("Return code from validator {}".format(p.returncode))

        vcf_info = {
            'version': vcf_version,
            'contigs': vcf_contigs,
            'genotype_ids': vcf_genotypes,
            'chromosome_ids': vcf_chromosomes
        }

        return validation_output_filename, vcf_info

    def _save_staging_to_local(self, ctx, params):
        dl_vcf = self.dfu.download_staging_file({'staging_file_subdir_path' : params['vcf_staging_file_path']})
        vcf_local_file_path = dl_vcf.get('copy_file_path')

        return vcf_local_file_path

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
            raise ValueError('VCF chromosome ids do not correspond to chromosome IDs from the assembly master list')

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
        file_valid_result, vcf_info = self.validate_vcf(ctx, params)

        assembly_chr_ids = self._validate_assembly_ids(ctx, params, vcf_info['chromosome_ids'])

        sample_ids = self._validate_sample_ids(ctx, params, vcf_info['genotype_ids'])

        return 'cool'