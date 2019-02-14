import uuid
import shutil
import os
import sys
import subprocess
import logging
import time
import binascii
import vcf
import gzip
import json
from pprint import pprint as pp

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

logging.basicConfig(format='%(created)s %(levelname)s: %(message)s')
def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

class VCFToVariation:
    def __init__(self, callback_url, scratch):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)
        #self.wsc = Workspace(callback_url)
        self.wsc = Workspace("https://appdev.kbase.us/services/ws")

    def _parse_vcf_data(self, ctx, params):
        try:
            if ctx['test_env']:
                # in the testing environment, vcf_staging_file_path set to local
                # test vcf file
                vcf_filepath = params['vcf_staging_file_path']
            else:
                # otherwise extract the file location from the input ui parameters
                vcf_filepath = self._stage_input(ctx, params)
        except KeyError:
            vcf_filepath = self._stage_input(ctx, params)

        # file is validated by this point, can assume vcf_filepath is valid
        reader = vcf.Reader(open(vcf_filepath, 'r'))

        version = float(reader.metadata['fileformat'][4:6])
        genotypes = reader.samples
        chromosomes = []
        contigs = {}
        totalvars = 0

        for record in reader:
            totalvars += 1
            if record.CHROM not in chromosomes:
                chromosomes.append(record.CHROM)

            if record.CHROM not in contigs.keys():
                passvar = 1 if not record.FILTER else 0

                contigs[record.CHROM] = {
                    'contig_id': record.CHROM,
                    'totalvariants': 1,
                    'passvariants': passvar,
                    'length': str(record.affected_end-record.affected_start),
                }
            else:
                contigs[record.CHROM]['totalvariants'] += 1
                if not record.FILTER:
                    contigs[record.CHROM]['passvariants'] += 1

        vcf_info = {
            'version': version,
            'contigs': contigs,
            'total_variants': totalvars,
            'genotype_ids': genotypes,
            'chromosome_ids': chromosomes,
            'file_ref': vcf_filepath
        }

        return vcf_info

    def _validate_vcf_to_sample(self, vcf_genotypes, sample_ids):
        check = True
        gid = None

        vgenotypes = [x.upper().strip() for x in vcf_genotypes]
        sids = [x.upper().strip() for x in sample_ids]

        for geno in vgenotypes:
            if geno not in sids:
                check = False
                gid = geno
                break

        return check, gid

    def _chk_if_vcf_ids_in_assembly(self, vcf_chromosomes, assembly_chromosomes):
        check = True

        for chromo in vcf_chromosomes:
            if chromo not in assembly_chromosomes:
                check = False
                break

        return check

    def _get_vcf_version(self, vcf_filepath):
        with(gzip.open if is_gz_file(vcf_filepath) else open)(vcf_filepath, 'rt') as vcf:
            line = vcf.readline()
            tokens = line.split('=')

            if not (tokens[0].startswith('##fileformat')):
                log("Invalid VCF.  ##fileformat line in meta is improperly formatted.")
                raise ValueError("Invalid VCF.  ##fileformat line in meta is improperly formatted. "
                                 "Check VCF file specifications: https://samtools.github.io/hts-specs/")

            vcf_version = float(tokens[1][-4:].rstrip())

            return vcf_version

    def validate_vcf(self, ctx, params):
        if 'genome_ref' not in params:
            raise ValueError('Genome reference not in input parameters: \n\n'+params)
        if 'vcf_staging_file_path' not in params:
            raise ValueError('VCF staging file path not in input parameters: \n\n' + params)

        try:
            if ctx['test_env']:
                # in the testing environment, vcf_staging_file_path set to local
                # test vcf file
                vcf_filepath = params['vcf_staging_file_path']
            else:
                # otherwise extract the file location from the input ui parameters
                vcf_filepath = self._stage_input(ctx, params)
        except KeyError:
            vcf_filepath = self._stage_input(ctx, params)

        vcf_version = self._get_vcf_version(vcf_filepath)

        # setup directorys for validation output
        validation_output_dir = os.path.join(self.scratch, 'validation_' + str(uuid.uuid4()))
        os.mkdir(validation_output_dir)

        # vcftools (vcf-validator) supports VCF v4.0-4.2
        # https://github.com/vcftools/vcftools

        # EBIvariation/vcf-validator (vcf_validator_linux) supports VCF v4.1-4.3
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
        elif vcf_version >= 4.0:
            print("Using vcftools to validate...")
            validator_cmd = ["vcf-validator"]
            validator_cmd.append(vcf_filepath)
            print("VCF version 4.0.")
        else:
            raise ValueError('VCF Version not in file, or fileformat line malformatted, or not version >=4.0. file format line must be the '
                             'first line of vcf file and in appropriate syntax. Check VCF file specifications: '
                             'https://samtools.github.io/hts-specs/')

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
            if line.decode("utf-8").strip().startswith('[info]'):
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

        return validation_output_filename

    def _stage_input(self, ctx, params):
        # extract file location from input ui parameters
        staging_dir = '/staging'
        vcf_local_file_path = os.path.join(staging_dir, params['vcf_staging_file_path'])

        if not os.path.exists(vcf_local_file_path):
            raise OSError('Staging input file does not exists, or is not readable')
        
        # TODO: use data file utils here, upload vcf to shock, use dfu.
        if is_gz_file(vcf_local_file_path):
            # /staging is read only, therefore have to copy before uncompressing
            copy = shutil.copy(vcf_local_file_path, os.path.join(self.scratch,params['vcf_staging_file_path']))
            # TODO: alter DFU to accept option output dir?
            unpack = self.dfu.unpack_file({'file_path': copy})
            return unpack['file_path']
        else:
            return vcf_local_file_path

    def _validate_assembly_ids(self, ctx, params, vcf_info):
        # All chromosome ids from the vcf should be in assembly
        # but not all assembly chromosome ids should be in vcf


        subset = self.wsc.get_object_subset([{
            'included': ['/assembly_ref'],
            'ref': params['genome_ref']
        }])

        vcf_info['assembly_ref'] = subset[0]['data']['assembly_ref']

        assembly_chromosome_ids_call = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': vcf_info['assembly_ref']
        }])

        assembly_chromosomes = assembly_chromosome_ids_call[0]['data']['contigs'].keys()
        vcf_chromosomes = vcf_info['chromosome_ids']

        if not self._chk_if_vcf_ids_in_assembly(vcf_chromosomes, assembly_chromosomes):
            raise ValueError('VCF chromosome ids do not correspond to chromosome IDs from the assembly master list')

        return assembly_chromosomes

    def _validate_sample_ids(self, ctx, params, vcf_info):
        # All samples within the VCF file need to be in sample attribute list
        vcf_genotypes = vcf_info['genotype_ids']

        sample_ids_subset = self.wsc.get_object_subset([{
            'included': ['/instances'],
            'ref': params['sample_attribute_ref']
        }])

        sample_ids = sample_ids_subset[0]['data']['instances'].keys()

        chk, chkerror = self._validate_vcf_to_sample(vcf_genotypes, sample_ids)

        if not chk:
            raise ValueError('VCF file genotype id: '+str(chkerror)+' is not listed in the sample meta data')

        return sample_ids

    def _construct_contig_info(self, ctx, params, vcf_info):
        """
            KBaseGwasData.Variations type spec

            /*
               Contig variation data
                 contig_id - contig identifier
                 totalvariants - total number of variants in each contig
                 passvariants - total number of variants that pass quality variation filter in contig
                 length - length of contig from assembly data
             */

             typdef structure {
               string contig_id;
               int totalvariants;
               int passvariants;
               int length; // from assembly
             } contig_info;
        """
        contigs = []

        contig_infos = vcf_info['contigs']

        for variant in contig_infos:
            contigs.append(contig_infos[variant])

        return contigs

    def _construct_variation(self, ctx, params, contigs_info, vcf_info):
        """
            KBaseGwasData.Variations type spec
             /*
               Variation object data structure
                 num_genotypes - number of total genotypes within variant file
                 num_variants - number of total variants within variant file
                 contigs - list of contig ids and variant information
                 attribute_ref - KBase reference to attribute mapping workspace object
                 genome_ref - KBase reference to genome workspace object
                 assembly_ref - KBase reference to assemebly workspace object
                 vcf_handle_ref - VCF handle reference to VCF file

                 @optional genome_ref
             */
             typedef structure {
               int numgenotypes;
               int numvariants;
               list<contig_info> contigs;
               attribute_ref population; // KBaseExperiments.AttributeMapping
               genome_ref genome_ref; // KBaseGenomes.Genome
               assembly_ref assemby_ref; // KBaseGenomeAnnotations.Assembly
               vcf_handle_ref vcf_handle_ref;
             } Variations;

            :param ctx: KBase context reference
            :param params: KBase ui input parameters
            :param population: previoiusly constructed sample population data
            :return: constructed variation object (dictionary)
        """

        if vcf_info['file_ref'].startswith(self.scratch):
            vcf_shock_file_ref = self.dfu.file_to_shock({'file_path': vcf_info['file_ref'], 'make_handle': 1})
        else:
            new_vcf_file = os.path.join(self.scratch, os.path.basename(vcf_info['file_ref']))
            new_vcf = shutil.copy(vcf_info['file_ref'], new_vcf_file)

            vcf_shock_file_ref = self.dfu.file_to_shock({'file_path': new_vcf, 'make_handle': 1})

        if not vcf_shock_file_ref['shock_id']:
            raise ValueError('Unable to upload VCF to Shock!')

        variation_obj = {
            'numgenotypes': int(len(vcf_info['genotype_ids'])),
            'numvariants': int(vcf_info['total_variants']),
            'contigs': contigs_info,
            'population': params['sample_attribute_ref'],
            'genome_ref': params['genome_ref'],
            'assembly_ref': vcf_info['assembly_ref'],
            'vcf_handle_ref': vcf_shock_file_ref['handle']
        }

        return variation_obj

    def _save_var_obj(self, ctx, params, var):
        print('Saving Variation to workspace...\n')
        sys.stdout.flush()

        if 'variation_object_name' not in params or params['variation_object_name'] is None or params['variation_object_name'] == '':
            var_obj_name = 'Variation_'+str(uuid.uuid4())

        if var:
            if not 'variation_object_name' in params:
                var_obj_name = 'variation_'+str(uuid.uuid4())
            else:
                var_obj_name = params['variation_object_name']

            """
            var_obj_info = self.dfu.save_objects({
                'id': self.dfu.ws_name_to_id(params['workspace_name']),
                'objects': [{
                    'type': 'KBaseGwasData.Variations',
                    'data': var,
                    'name': var_obj_name
                }]
            })[0]
            """
            # varjson = json.loads(var)

            with open(os.path.join(self.scratch, 'var-data.json'), 'w') as f:
                json.dump(var, f,indent=4)

            pp(var)

            return '1/2/3'
        else:
            raise ValueError('Variation object blank, cannot not save to workspace!')

    def import_vcf(self, ctx, params):
        # VCF validation
        # VCF file validation
        file_valid_result = self.validate_vcf(ctx, params)
        # VCF file parsing
        vcf_info = self._parse_vcf_data(ctx, params)
        # Validate vcf chromosome ids against assembly chromosome ids
        self._validate_assembly_ids(ctx, params, vcf_info)
        # Validate vcf genotypes against sample meta data ids
        self._validate_sample_ids(ctx, params, vcf_info)

        # Variation object construction
        # construct contigs_info
        contigs_info = self._construct_contig_info(ctx, params, vcf_info)
        # construct variation
        var = self._construct_variation(ctx, params, contigs_info, vcf_info)

        # Save variation object to workspace
        var_wksp_obj = self._save_var_obj(ctx, params, var)

        return var_wksp_obj
