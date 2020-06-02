import hashlib
import logging
import os
import gzip
import re
import time

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


def md5_sum_local_file(fname):
    md5hash = hashlib.md5()
    with open(fname, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5hash.update(chunk)
    return md5hash.hexdigest()


def compare_md5_local_with_shock(fname, shock_file_ref):
    local_md5 = md5_sum_local_file(fname)

    shock_md5 = shock_file_ref['handle']['remote_md5']

    if local_md5 != shock_md5:
        raise ValueError(f'Local md5 {local_md5} does not match shock md5 {shock_md5}')

    if not shock_file_ref['shock_id']:
        raise ValueError('Unable to upload assembly index to Shock!')


class VCFToVariation:
    def __init__(self, Config):
        self.scratch = Config['scratch']
        ws_url = Config['ws_url']
        callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(callback_url)
        self.wsc = Workspace(ws_url)
        self.au = AssemblyUtil(callback_url)
        self.vcf_info = dict()

    def _parse_header(self, record, category):
        """
        parses vcf header which looks like the following
        and get details for the IDs like DP, q10
        This information is useful in filtering
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
        ##FILTER=<ID=q10,Description="Quality below 10">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        [
            {
            },
            {
            },
            {
            }
        ]
        """
        returninfo = {"Category": category}
        info = re.sub(".*<", "", record)
        info = info.replace(">", "")
        infolist = info.replace('"', '').rstrip().split(",")
        for fields in infolist:
            data = fields.split("=")
            key = data.pop(0)
            val = "=".join(data)
            val = val.replace("\"", "")
            returninfo[key] = val
        return returninfo

    def parse_vcf_data(self, vcf_filepath):
        """
        parses vcf file including headers and prepares
        information that will be uploaded to KBase workspace
        :param vcf_filepath:
        :return:
        """
        reader = gzip.open(vcf_filepath, "rt")
        version = ""
        genotypes = ""
        #
        counter = 0
        chromosomes = []
        contigs = {}
        header = list()
        totalvars = 0

        for record in reader:

            # Handle header lines and parse information
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

            # Handle the actual VCF content and parese information
            counter = counter + 1
            CHROM, POS, *r = record.split("\t")
            totalvars += 1
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

    def _validate_vcf_to_sample(self, vcf_genotypes, sample_ids):
        genos_not_found = []
        vgenotypes = [x.upper().strip() for x in vcf_genotypes]
        sids = [x.upper().strip() for x in sample_ids]

        for geno in vgenotypes:
            if geno not in sids:
                genos_not_found.append(geno)

        if not genos_not_found:
            return True
        else:
            return genos_not_found

    def _chk_if_vcf_ids_in_assembly(self, vcf_chromosomes, assembly_chromosomes):
        """
        Check if all chromosome ids in vcf are also present in assembly
        :param vcf_chromosomes:
        :param assembly_chromosomes:
        :return: returns list of chromosome ids,
                present in vcf, but absent from assembly
        """
        chromos_not_in_assembly = []

        for chromo in vcf_chromosomes:
            if chromo not in assembly_chromosomes:
                chromos_not_in_assembly.append(chromo)

        if not chromos_not_in_assembly:
            return True
        else:
            return chromos_not_in_assembly

    def _validate_assembly_ids(self, vcf_info):
        """
        All chromosome ids from the vcf should be in assembly
        but not all assembly chromosome ids need to be in vcf
        :param params:
        :return: list of all assembly chromosome ids
        """
        assembly_chromosome_ids_call = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': vcf_info['assembly_ref']
        }])
        assembly_chromosomes = assembly_chromosome_ids_call[0]['data']['contigs'].keys()
        vcf_chromosomes = vcf_info['chromosome_ids']
        chk_assembly_ids = self._chk_if_vcf_ids_in_assembly(vcf_chromosomes, assembly_chromosomes)

        if isinstance(chk_assembly_ids, list):
            failed_ids = ' '.join(chk_assembly_ids)
            print(f'VCF contig ids: {failed_ids} are not present in assembly.')
            raise ValueError(f'VCF contig ids: {failed_ids} are not present in assembly.')

        return assembly_chromosomes

    def _validate_sample_ids(self, vcf_info, sample_attribute_ref):
        # All samples within the VCF file need to be in sample attribute list
        # Sample attribute ref is not mandatory anymore

        vcf_genotypes = vcf_info['genotype_ids']
        sample_ids_subset = self.wsc.get_object_subset([{
            'included': ['/instances'],
            'ref': sample_attribute_ref
        }])
        sample_ids = sample_ids_subset[0]['data']['instances'].keys()
        validate_genotypes = self._validate_vcf_to_sample(vcf_genotypes, sample_ids)
        if isinstance(validate_genotypes, list):
            failed_genos = ' '.join(validate_genotypes)
            print(f'VCF genotypes: {failed_genos} are not present in sample attribute mapping.')
            raise ValueError(f'VCF genotypes: {failed_genos} are not present in sample attribute mapping.')
        else:
            return sample_ids

    def _construct_contig_info(self, vcf_info):
        """
           From KBaseGwasData.Variations type spec
            /*
               Contig variation data
                 contig_id - contig identifier
                 totalvariants - total number of variants in each contig
                 length - length of contig from assembly data
             */
             typdef structure {
               string contig_id;
               int totalvariants;
               int length; // from assembly
             } contig_info;

        """
        assembly_chromosome_dict = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': vcf_info['assembly_ref']
        }])[0]['data']['contigs']

        contigs = []

        contig_infos = vcf_info['contigs']

        for contig_id in contig_infos:
            length_contig = assembly_chromosome_dict[contig_id].get("length")
            contig_infos[contig_id]["length"] = length_contig
            contigs.append(contig_infos[contig_id])

        return contigs

    def _construct_variation_object_json(self, vcf_info):

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
                 samples

                 @optional genome_ref
             */
             typedef structure {
               int numgenotypes;
               int numvariants;
               list<contig_info> contigs;
               attribute_ref samples; // KBaseExperiments.AttributeMapping
               genome_ref genome_ref; // KBaseGenomes.Genome
               assembly_ref assembly_ref; // KBaseGenomeAnnotations.Assembly
               vcf_handle_ref vcf_handle_ref;
             } Variations;

            :param params: KBase ui input parameters
            :return: constructed variation object (dictionary)
        """
        logging.info("Uploading VCF file to shock")
        vcf_shock_file_ref = None
        vcf_index_shock_file_ref = None
        if os.path.exists(vcf_info['vcf_compressed']):
            vcf_shock_file_ref = self.dfu.file_to_shock({
                'file_path': vcf_info['vcf_compressed'],
                'make_handle': 1
            })
        # compare_md5_local_with_shock(bgzip_file_path, vcf_shock_file_ref)

        logging.info("Uploading VCF index file to shock")
        if os.path.exists(vcf_info['vcf_index']):
            vcf_index_shock_file_ref = self.dfu.file_to_shock({
                'file_path': vcf_info['vcf_index'],
                'make_handle': 1
            })
        # compare_md5_local_with_shock(index_file_path, vcf_index_shock_file_ref)

        variation_obj_data = {
            'numgenotypes': int(len(vcf_info['genotype_ids'])),
            'numvariants': int(vcf_info['total_variants']),
            'contigs': vcf_info['contigs_info'],
            'samples': vcf_info['genotype_ids'],
            "header": vcf_info['header'],
            'assembly_ref': vcf_info['assembly_ref'],
            'vcf_handle_ref': vcf_shock_file_ref['handle']['hid'],
            'vcf_handle': vcf_shock_file_ref['handle'],
            'vcf_index_handle_ref': vcf_index_shock_file_ref['handle']['hid'],
            'vcf_index_handle': vcf_index_shock_file_ref['handle'],
        }
        if 'genome_ref' in vcf_info:
            variation_obj_data['genome_ref'] = vcf_info['genome_ref']
        if 'sample_attribute_ref' in vcf_info:
            variation_obj_data['sample_attribute_ref'] = vcf_info['sample_attribute_ref']
        return variation_obj_data

    def generate_variation_object_data (self, params):
        # VCF validation
        # VCF file parsing

        # Copy vcf_compressed, vcf_index,

        vcf_info = self.parse_vcf_data(params['vcf_compressed'])
        vcf_info['vcf_compressed'] = params['vcf_compressed']
        vcf_info['vcf_index'] = params['vcf_index']

        assembly_ref = params['assembly_ref']
        if 'genome_ref' in params:
            genome_ref = params['genome_ref']

        logging.info("Parsing vcf started")
        vcf_info['assembly_ref'] = assembly_ref
        if 'genome_ref' in params:
            vcf_info['genome_ref'] = genome_ref

        logging.info("Comparing assembly ids")
        # Validate vcf chromosome ids against assembly chromosome ids
        result = self._validate_assembly_ids(vcf_info)
        # Variation object construction
        # construct contigs_info
        if result:
            logging.info("Creating contig info")
            vcf_info['contigs_info'] = self._construct_contig_info(vcf_info)

        logging.info("Validating sample ids")
        # Validate vcf genotypes against sample meta data ids
        # provided in sample_attribute_ref
        if 'sample_attribute_ref' in params:
            result = self._validate_sample_ids(vcf_info, params['sample_attribute_ref'])
            if result:
                vcf_info['sample_attribute_ref'] = params['sample_attribute_ref']

        # construct variation
        variation_object_json = self._construct_variation_object_json(vcf_info)

        return variation_object_json

