# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def import_variation(self, ctx, import_variation_params):
        """
        :param import_variation_params: instance of type
           "import_variation_params" (TODO: take assembly ref as input take
           care of metagenomic use case required params: genome_ref:
           KBaseGenomes.Genome object reference *** variation input data ***
           vcf_staging_file_path: path to location data associated with
           samples variation_object_name: output name for KBase variation
           object *** sample input data *** sample_attribute_ref: x/y/z
           reference to kbase sample attribute optional params: ***
           Visualization *** plot_maf: generate histogram of minor allele
           frequencies plot_hwe: generate histogram of Hardy-Weinberg
           Equilibrium p-values) -> structure: parameter "workspace_name" of
           String, parameter "genome_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "vcf_staging_file_path" of type "filepath"
           (KBase file path to staging files), parameter
           "variation_object_name" of String, parameter
           "sample_attribute_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "import_variation_results" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: returnVal
        #BEGIN import_variation
        print (import_variation_params);

        returnVal = {
			'report_ref': "report_ref",
			'report_name': 'report_name'
		}

        #END import_variation

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method import_variation return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
