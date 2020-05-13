import uuid
import shutil
import os
import subprocess
import logging
import time
import binascii
import vcf
import gzip
import hashlib
from pprint import pprint as pp

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.GenericsAPIClient import GenericsAPI


#logging.basicConfig(format='%(created)s %(levelname)s: %(message)s')
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


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
    def __init__(self, config, scratch, callback_url ):
        self.scratch = config['scratch']
        self.ws_url = config['workspace-url']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(self.callback_url)
        self.wsc = Workspace(self.ws_url)
        self.scratch = scratch
        self.callback_url = callback_url
        self.au = AssemblyUtil(self.callback_url)
        self.gapi = GenericsAPI(self.callback_url)


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
        chromos_not_in_assembly = []

        for chromo in vcf_chromosomes:
            if chromo not in assembly_chromosomes:
                chromos_not_in_assembly.append(chromo)

        if not chromos_not_in_assembly:
            return True
        else:
            return chromos_not_in_assembly

    def _create_sample_attribute_file(self, vcf_file, sample_attribute_mapping_file):
        """
        function for creating sample attribute mapping file.
        """
        try:
            with open (vcf_file, 'r') as vcf_handle:
                Lines = vcf_handle.readlines()

                for line in Lines:
                    if(line.startswith("#CHROM")):
                       header = line.lstrip().split("\t")

                       try:
                          with open (sample_attribute_mapping_file, 'w') as attribute_mapping_handle:
                              attribute_mapping_handle.write("Attribute\tAttribute ontology ID\tUnit\tUnit ontology ID")

                              for i in range(9,len(header)):
                                  attribute_mapping_handle.write("\t"+header[i])
                              #attribute_mapping_handle.write("\n")


                              attribute_mapping_handle.write("label\t\t\t")
                              for j in range(9,len(header)):
                                  attribute_mapping_handle.write("\t"+header[j])
                              #attribute_mapping_handle.write("\n")
                       except IOError:
                           print("Could not write to file:", sample_attribute_mapping_file)

        except IOError:
               print("Could not read file:", vcf_file)

    def _validate_assembly_ids(self, params):
        # All chromosome ids from the vcf should be in assembly
        # but not all assembly chromosome ids should be in vcf


        if ('genome_ref' in params):
            subset = self.wsc.get_object_subset([{
                'included': ['/assembly_ref'],
                'ref': params['genome_or_assembly_ref']
            }])

            self.vcf_info['assembly_ref'] = subset[0]['data']['assembly_ref']

        if ('assembly_ref' in params):
            self.vcf_info['assembly_ref'] = params['assembly_ref']

        assembly_chromosome_ids_call = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': self.vcf_info['assembly_ref']
        }])

        assembly_chromosomes = assembly_chromosome_ids_call[0]['data']['contigs'].keys()
        vcf_chromosomes = self.vcf_info['chromosome_ids']

        chk_assembly_ids =  self._chk_if_vcf_ids_in_assembly(vcf_chromosomes, assembly_chromosomes)

        if isinstance(chk_assembly_ids, list):
            failed_ids = ' '.join(chk_assembly_ids)
            print(f'VCF contig ids: {failed_ids} are not present in assembly.')
            raise ValueError(f'VCF contig ids: {failed_ids} are not present in assembly.')


        return assembly_chromosomes

    def _validate_sample_ids(self, params):
        # All samples within the VCF file need to be in sample attribute list


        vcf_genotypes = self.vcf_info['genotype_ids']

        sample_ids_subset = self.wsc.get_object_subset([{
            'included': ['/instances'],
            'ref': params['sample_attribute_ref']
        }])

        sample_ids = sample_ids_subset[0]['data']['instances'].keys()

        validate_genotypes = self._validate_vcf_to_sample(vcf_genotypes, sample_ids)

        if isinstance(validate_genotypes, list):
            failed_genos = ' '.join(validate_genotypes)
            print(f'VCF genotypes: {failed_genos} are not present in sample attribute mapping.')
            raise ValueError(f'VCF genotypes: {failed_genos} are not present in sample attribute mapping.')

        return sample_ids

    def _construct_contig_info(self, params):
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

        assembly_chromosome_dict = self.wsc.get_object_subset([{
            'included': ['/contigs'],
            'ref': self.vcf_info['assembly_ref']
        }])[0]['data']['contigs']


        contigs = []

        contig_infos = self.vcf_info['contigs']


        for contig_id in contig_infos:
            length_contig = assembly_chromosome_dict[contig_id].get("length")
            contig_infos[contig_id]["length"] = length_contig
            contigs.append(contig_infos[contig_id])

        return contigs
   

    def _index_assembly(self, assembly_file):
        if not os.path.exists(assembly_file):
           logging.info (assembly_file + " does not exist")

        logging.info("indexing assembly file")

        assembly_index_cmd = ["samtools", "faidx", assembly_file]
        print(assembly_index_cmd)
        p = subprocess.Popen(assembly_index_cmd,
                             cwd=self.scratch,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=False)

        out, err = p.communicate()

        logging.info("indexing of assembly file done!")

        return assembly_file + ".fai"

    def _download_assembly(self, assembly_ref):
        file = self.au.get_assembly_as_fasta({
          'ref': assembly_ref
        })
       #file = "/kb/module/work/tmp/Athaliana_TAIR10.assembly.fa"
        return file
 
    def _construct_variation(self, params, contigs_info):
        
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
               attribute_ref population; // KBaseExperiments.AttributeMapping
               genome_ref genome_ref; // KBaseGenomes.Genome
               assembly_ref assemby_ref; // KBaseGenomeAnnotations.Assembly
               vcf_handle_ref vcf_handle_ref;
             } Variations;

            :param params: KBase ui input parameters
            :param population: previoiusly constructed sample population data
            :return: constructed variation object (dictionary)
        """

        logging.info("Uploading VCF file to shock")
        bgzip_file_path = params['vcf_local_file_path']
        vcf_shock_file_ref = self.dfu.file_to_shock(
            {'file_path': bgzip_file_path, 'make_handle': 1}
        )
        #compare_md5_local_with_shock(bgzip_file_path, vcf_shock_file_ref)

        logging.info("Uploading VCF index file to shock")
        index_file_path = params['vcf_index_file_path']
        vcf_index_shock_file_ref = self.dfu.file_to_shock(
            {'file_path': index_file_path, 'make_handle': 1}
        )
        #compare_md5_local_with_shock(index_file_path, vcf_index_shock_file_ref)

      #  assembly_file_path = self._download_assembly(self.vcf_info['assembly_ref'])['path']
      #  assembly_index_file_path = self._index_assembly(assembly_file_path)
      #  assembly_index_shock_file_ref = self.dfu.file_to_shock(
      #      {'file_path': assembly_index_file_path, 'make_handle': 1}
      #  )
        #compare_md5_local_with_shock(assembly_index_file_path, assembly_index_shock_file_ref)
        variation_obj = {
            'numgenotypes': int(len(self.vcf_info['genotype_ids'])),
            'numvariants': int(self.vcf_info['total_variants']),
            'contigs': contigs_info,
            'population': params['sample_attribute_ref'],
            'samples': self.vcf_info['genotype_ids'],
            "header": self.vcf_info['header'],

            # TODO: TYPE SPEC CHANGE: need to change type spec to assembly_ref instead of assemby_ref
            'assemby_ref': self.vcf_info['assembly_ref'],
            'vcf_handle_ref': vcf_shock_file_ref['handle']['hid'],
            'vcf_handle' : vcf_shock_file_ref['handle'],
            'vcf_index_handle_ref': vcf_index_shock_file_ref['handle']['hid'],
            'vcf_index_handle': vcf_index_shock_file_ref['handle'],
       #     'assembly_index_handle_ref': assembly_index_shock_file_ref['handle']['hid'],
       #     'assembly_index_handle': assembly_index_shock_file_ref['handle']
        }
        if 'genome_ref' in params:
            variation_obj['genome_ref'] =  params['genome_ref']

        return variation_obj

    def _save_var_obj(self, params, var):
        """
        :param params:
        :param var:
        :return:
            DataFileUtils object_info:
                objid - the numerical id of the object.
                name - the name of the object.
                type - the type of the object.
                save_date - the save date of the object.
                ver - the version of the object.
                saved_by - the user that saved or copied the object.
                wsid - the id of the workspace containing the object.
                workspace - the name of the workspace containing the object.
                chsum - the md5 checksum of the object.
                size - the size of the object in bytes.
                meta - arbitrary user-supplied metadata about the object.
        """

        logging.info('Saving Variation to workspace...\n')

        if var:
            if not 'variation_object_name' in params:
                var_obj_name = 'variation_'+str(uuid.uuid4())
            else:
                var_obj_name = params['variation_object_name']

            var_obj_info = self.dfu.save_objects({
                'id': self.dfu.ws_name_to_id(params['workspace_name']),
                'objects': [{
                    'type': 'KBaseGwasData.Variations',
                    'data': var,
                    'name': var_obj_name
                }]
            })[0]

            return var_obj_info
        else:
            raise ValueError('Variation object blank, cannot not save to workspace!')

    def _validate_sample_attribute_ref(self, params):

        #params["sample_attribute_ref"] = ''  #just for testing
        if not params['sample_attribute_ref']:
           sample_attribute_mapping_file = os.path.join(self.scratch ,"sample_attribute.tsv")   #hardcoded for testing
           self._create_sample_attribute_file(params['vcf_local_file_path'], sample_attribute_mapping_file)
          
           logging.info("Uploading sample attribute file to shock")
           vcf_sample_attribute_shock_file_ref = self.dfu.file_to_shock(
               {'file_path': sample_attribute_mapping_file, 'make_handle': 1}
           )
           shock_id = vcf_sample_attribute_shock_file_ref['shock_id']
           ws_id = self.dfu.ws_name_to_id(params['workspace_name'])
           import_params = {
                  'input_shock_id' : shock_id,
                  'output_ws_id': ws_id,
                  'output_obj_name': 'Sample_attribute'}

           ret = self.gapi.file_to_attribute_mapping(import_params)
           params['sample_attribute_ref'] = ret['attribute_mapping_ref']

    def import_vcf(self, params, vcf_info):
        # VCF validation
        # VCF file validation
        #file_valid_result = self.validate_vcf(params)
        logging.info("Validating sample attributes")
        self._validate_sample_attribute_ref(params)
        # VCF file parsing
        logging.info("Parsing vcf started")
        self.vcf_info = vcf_info
        logging.info("Comparing assembly ids")
        # Validate vcf chromosome ids against assembly chromosome ids
        self._validate_assembly_ids(params)

        logging.info("Validating sample ids")
        # Validate vcf genotypes against sample meta data ids
        self._validate_sample_ids(params)

        # Variation object construction
        # construct contigs_infoa
        logging.info("Creating contig info")
        contigs_info = self._construct_contig_info(params)

        logging.info("Adding chromoosome length")
        # construct variation
        var = self._construct_variation(params, contigs_info)

        logging.info("Saving variation object to workspace")
        # Save variation object to workspace
        var_wksp_obj = self._save_var_obj(params, var)

        return [var_wksp_obj, var]
