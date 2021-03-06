# -*- coding: utf-8 -*-
import os
import time
import unittest
from configparser import ConfigParser

from VariationUtil.VariationUtilImpl import VariationUtil
from VariationUtil.Util.VariationToVCF import VariationToVCF
from VariationUtil.Util.VCFToVariation import VCFToVariation
from VariationUtil.VariationUtilServer import MethodContext
from VariationUtil.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace


class VariationUtilTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('VariationUtil'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'VariationUtil',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1,
                        'test_env': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = VariationUtil(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
#        cls.VCFtoVar = VCFToVariation(cls.cfg, "/kb/module/work/tmp", cls.callback_url)
        cls.vcf_test_dir = '/kb/module/test/sample_data/vcf'

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    #@unittest.skip('x')
    def test_save_variation_genome_ref(self):
        ret = self.serviceImpl.save_variation_from_vcf(self.ctx, {
            'workspace_name': 'pranjan77:narrative_1594999233120',
            'genome_or_assembly_ref': '52584/2/1',
            'sample_set_ref':'52586/7/1',
            'vcf_staging_file_path': '/kb/module/test/sample_data/small_poplar/test.vcf',
            'sample_attribute_ref': None,
            'variation_object_name': 'poplar_g2'
        })

    #'sample_attribute_ref': '39465/3/1',

    @unittest.skip('x')
    def test_save_variation_assembly_ref(self):
        ret = self.serviceImpl.save_variation_from_vcf(self.ctx, {
            'workspace_name': 'pranjan77:narrative_1594999233120',
            'genome_or_assembly_ref': '52584/2/1',
            'vcf_staging_file_path': '/kb/module/test/sample_data/small_poplar/small_poplar.vcf.gz',
            'variation_object_name': 'poplar_a'
        })

    @unittest.skip('x')
    def test_get_variation_as_vcf_1(self):
        ret = self.serviceImpl.get_variation_as_vcf(self.ctx, {
            'variation_ref': '42434/2/7',
            'filename':"small_poplar.vcf"
        })
        print (ret)


    """
    def test_vcf_validator_linux_pass(self):
        file_validation = self.VCFtoVar.validate_vcf(self.ctx, {'workspace_name': 'pranjan77:narrative_1549050842078',
                                                             'genome_ref': '24237/5/8',
                                                             'vcf_staging_file_path' : '/kb/module/test/sample_data/vcf/v4.3/pass/complexfile_passed_000.vcf',
                                                             'sample_attribute_ref' : '24237/17/1',
                                                             'variation_object_name' : 'arabidopsis_variation'})

    def test_vcf_validator_linux_fail(self):
        with self.assertRaises(ValueError):
            self.VCFtoVar.validate_vcf(self.ctx, {'workspace_name': 'pranjan77:narrative_1549050842078',
                                                     'genome_ref': '24237/5/8',
                                                     'vcf_staging_file_path' : '/kb/module/test/sample_data/vcf/v4.3/fail/failed_body_alt_000.vcf',
                                                     'sample_attribute_ref' : '24237/17/1',
                                                     'variation_object_name' : 'arabidopsis_variation'
            })

    def test_vcftools_pass(self):
        file_validation = self.VCFtoVar.validate_vcf(self.ctx, {'workspace_name': 'pranjan77:narrative_1549050842078',
                                                             'genome_ref': '24237/5/8',
                                                             'vcf_staging_file_path' : '/kb/module/test/sample_data/vcf/v4.0/pass/valid-4.0.vcf',
                                                             'sample_attribute_ref' : '24237/17/1',
                                                             'variation_object_name' : 'arabidopsis_variation'})

    
    def test_vcftools_fail(self):
        with self.assertRaises(ValueError):
            self.VCFtoVar.validate_vcf(self.ctx, {'workspace_name': 'pranjan77:narrative_1549050842078',
                                                     'genome_ref': '24237/5/8',
                                                     'vcf_staging_file_path' : '/kb/module/test/sample_data/vcf/v4.0/fail/invalid-4.0.vcf',
                                                     'sample_attribute_ref' : '24237/17/1',
                                                     'variation_object_name' : 'arabidopsis_variation'
            })
"""
