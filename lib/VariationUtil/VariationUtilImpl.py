# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import json
from VariationUtil.Util.VCFToVariation import VCFToVariation
from VariationUtil.Util.VCFUtils import VCFUtils
from VariationUtil.Util.VariationToVCF import VariationToVCF
from VariationUtil.Util.htmlreportutils import htmlreportutils
from VariationUtil.Util.StrainInfo import StrainInfo
from VariationUtil.Util.JbrowseUtil import JbrowseUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil

#END_HEADER


class VariationUtil:
    '''
    Module Name:
    VariationUtil

    Module Description:
    A KBase module: VariationUtil
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbasecollaborations/VariationUtil"
    GIT_COMMIT_HASH = "5c21f7b209448d534b4f4c1477d027046eb0247b"

    #BEGIN_CLASS_HEADER

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        # TODO: Make sure we need to define config just once
        # TODO: Change the code tp match this style
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.scratch = config['scratch']
        self.config['ws_url'] = config['workspace-url']

        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shared_folder = config['scratch']
        self.hr = htmlreportutils()
        self.ws_url = config['workspace-url']
        self.wsc = Workspace(self.ws_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.shock_url = config['shock-url']
        self.sw_url = config['srv-wiz-url']
        pass
        #END_CONSTRUCTOR
    def save_variation_from_vcf(self, ctx, params):
        """
        Save a variation (and trait?) object to Kbase given a reference genome, object output name,
        Variant Call Format (VCF) file, and sample attribute file.
        :param params: instance of type "save_variation_input" (## funcdef
           save_variation_from_vcf ## required input params:
           genome_or_assembly_ref: KBaseGenomes.Genome or
           KBaseGenomeAnnotations.Assembly object reference *** variation
           input data *** vcf_staging_file_path: path to location data
           associated with samples variation_object_name: output name for
           KBase variation object *** sample input data ***
           sample_attribute_ref: x/y/z reference to kbase sample attribute
           optional params: NA output report: report_name report_ref HTML
           visualization: Manhattan plot *** Visualization *** plot_maf:
           generate histogram of minor allele frequencies plot_hwe: generate
           histogram of Hardy-Weinberg Equilibrium p-values) -> structure:
           parameter "workspace_name" of String, parameter
           "genome_or_assembly_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "vcf_staging_file_path" of type "filepath"
           (KBase file path to staging files), parameter
           "variation_object_name" of String, parameter
           "sample_attribute_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "save_variation_output" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: report
        #BEGIN save_variation_from_vcf

        # Get workspace id
        ws_id = self.dfu.ws_name_to_id(params['workspace_name'])

        genome_ref = None
        assembly_ref = None

        # 1) Find whether the input is a genome or assembly
        #    and get genome_ref and assembly_ref

        genome_or_assembly_ref = params['genome_or_assembly_ref']
        obj_type = self.wsc.get_object_info3({
            'objects':[{
                'ref': genome_or_assembly_ref
                      }]})['infos'][0][2]
        if ('KBaseGenomes.Genome' in obj_type):
            genome_ref = genome_or_assembly_ref
            subset = self.wsc.get_object_subset([{
                    'included': ['/assembly_ref'],
                    'ref': genome_ref
                }])
            assembly_ref = subset[0]['data']['assembly_ref']
        elif ('KBaseGenomeAnnotations.Assembly' in obj_type):
            assembly_ref = genome_or_assembly_ref
        else:
            raise ValueError(obj_type + ' is not the right input for this method. '
                                      + 'Valid input include KBaseGenomes.Genome or '
                                      + 'KBaseGenomeAnnotations.Assembly ')


        # 2)  Validate VCF, compress, and build VCF index
        logging.info("Validating VCF, Compressing VCF and Indexing VCF")
        VCFUtilsConfig = {
            "scratch": self.scratch
        }
        VCFUtilsParams = {
            'vcf_staging_file_path': params['vcf_staging_file_path']
        }
        VCU = VCFUtils(VCFUtilsConfig)
        vcf_compressed, vcf_index, vcf_strain_ids = VCU.validate_compress_and_index_vcf(VCFUtilsParams)

        if vcf_index is not None:
            logging.info("vcf compressed :" + str(vcf_compressed))
            logging.info("vcf index :" + str(vcf_index))
            logging.info("vcf strain ids :" + str(vcf_strain_ids))
        else:
            raise ValueError("No result obtained after compression and indexing step")

        # Get strain info
        # TODO: Remove hard coded stuff
        StrainInfoConfig = self.config
        StrainInfoParams = {
            "ws_id": ws_id,
            "vcf_strain_ids": vcf_strain_ids,
            "sample_set_ref": params["sample_set_ref"],
            "sample_attribute_name": params["sample_attribute_name"]
        }
        si = StrainInfo(StrainInfoConfig)
        sample_attribute_ref, strains = si.sample_strain_info(StrainInfoParams)
        print (sample_attribute_ref)
        print (strains)


        # 3) Create json for variation object. In a following step genomic_indexes will be
        # added to this json before it is saved as Variation object

        VCFToVariationConfig = {
            "ws_url": self.ws_url,
            "scratch": self.scratch
        }
        VCFToVariationParams = {
            "vcf_compressed": vcf_compressed,
            "vcf_index": vcf_index,
            "assembly_ref": assembly_ref
        }
        if genome_ref is not None:
            VCFToVariationParams['genome_ref'] = genome_ref

        vtv = VCFToVariation(VCFToVariationConfig)
        variation_object_data = vtv.generate_variation_object_data(VCFToVariationParams)
        # Append sample information
        if sample_attribute_ref:
            variation_object_data['sample_attribute_ref'] = sample_attribute_ref
        else:
            raise ValueError(f'sample attribute ref not found')
        if strains:
            variation_object_data['strains'] = strains
        else:
            raise ValueError(f'strains not found')
        if 'sample_set_ref' in params:
            variation_object_data['sample_set_ref'] = params['sample_set_ref']
        else:
            raise ValueError(f'sample_set_ref not found in params')


        # 4)
        JbrowseConfig = {
            "ws_url": self.ws_url,
            "scratch": self.scratch,
            "sw_url": self.sw_url,
            "shock_url":self.shock_url
        }
        JbrowseParams = {
            "vcf_path": vcf_compressed,
            "assembly_ref": assembly_ref,
            "binsize": 10000,
            "vcf_shock_id": variation_object_data['vcf_handle']['id'],
            "vcf_index_shock_id":variation_object_data['vcf_index_handle']['id']
        }
        if genome_ref is not None:
            JbrowseParams["genome_ref"] = genome_ref

        jb = JbrowseUtil(JbrowseConfig)
        jbrowse_report = jb.prepare_jbrowse_report(JbrowseParams)


        # 5) Now we have the genomic indices and we have all the information needed to save
        # the variation object
        # TODO: Take out the genomic_indexes field from the object spec
        #  TODO: Take out the vcf_handle stuff not needed

        variation_object_data['genomic_indexes'] = jbrowse_report['genomic_indexes']

        # 6) We need to create a list of all handles needed and build the handles
        #    part of variation object
        #handles = list ()
        #handles.append(variation_object_data['vcf_handle'])
        #handles.append(variation_object_data['vcf_index_handle'])

        #for g in jbrowse_report['genomic_indexes']:
        #    handles.append(g)

        #variation_object_data['handles'] = handles
        #variation_object_data['handle'] = jbrowse_report['genomic_indexes'][0]


        #print (json.dumps(variation_object_data))

        var_obj = self.dfu.save_objects({
            'id': self.dfu.ws_name_to_id(params['workspace_name']),
            'objects': [{
                'type': 'KBaseGwasData.Variations',
                'data': variation_object_data,
                'name': params['variation_object_name']
            }]
        })[0]

        var_obj_ref = str(var_obj[6]) + "/" + str(var_obj[0]) + "/" + str(var_obj[4])
        print (var_obj_ref)


        # 5) Build jbrowse html report
        workspace = params['workspace_name']
        created_objects = []
        created_objects.append({
            "ref": var_obj_ref,
            "description": "Variation Object"
            })
        report = self.hr.create_html_report(jbrowse_report['jbrowse_data_path'],
                                            workspace,
                                            created_objects)
        report['variation_ref'] = var_obj_ref
        print(report)
        #END save_variation_from_vcf

        # At some point might do deeper type checking...
        if not isinstance(report, dict):
            raise ValueError('Method save_variation_from_vcf return value ' +
                             'report is not type dict as required.')
        # return the results
        return [report]

    def export_variation_as_vcf(self, ctx, params):
        """
        Export KBase variation object as Variant Call Format (VCF) file
        :param params: instance of type "export_variation_input" (## funcdef
           export_variation_as_vcf ## required input params: Variation object
           reference optional params: NA output report: Shock id pointing to
           exported vcf file) -> structure: parameter "input_var_ref" of type
           "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "export_variation_output" -> structure:
           parameter "shock_id" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN export_variation_as_vcf

        vtv = VariationToVCF(self.callback_url, self.shared_folder)
        output = vtv.export_as_vcf(params)

        #END export_variation_as_vcf

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method export_variation_as_vcf return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def get_variation_as_vcf(self, ctx, params):
        """
        Given a reference to a variation object, and output name: return a Variant Call Format (VCF)
        file path and name.
        :param params: instance of type "get_variation_input" (## funcdef
           get_variation_as_vcf ## required input params: Variation object
           reference output file name optional params: NA output report: path
           to returned vcf name of variation object) -> structure: parameter
           "variation_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "filename" of String
        :returns: instance of type "get_variation_output" -> structure:
           parameter "path" of type "filepath" (KBase file path to staging
           files), parameter "variation_name" of String
        """
        # ctx is the context object
        # return variables are: file
        #BEGIN get_variation_as_vcf
        vtv = VariationToVCF(self.callback_url, self.shared_folder)
        file = vtv.variation_to_vcf(params)

        #END get_variation_as_vcf

        # At some point might do deeper type checking...
        if not isinstance(file, dict):
            raise ValueError('Method get_variation_as_vcf return value ' +
                             'file is not type dict as required.')
        # return the results
        return [file]

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
