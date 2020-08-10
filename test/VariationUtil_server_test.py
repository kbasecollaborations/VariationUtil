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

    def test_save_variation_assembly_ref(self):
        ret = self.serviceImpl.save_variation_from_vcf(self.ctx, {
            'workspace_name': 'man4ish_gupta:narrative_1597069509290',
            'genome_or_assembly_ref': '52931/7/1',
            'sample_set_ref':'52931/5/1',
            'sample_attribute_name':'sample_attr_new',
            'vcf_staging_file_path': '/kb/module/test/sample_data/small_poplar/test3.vcf',
            'variation_object_name': 'poplar_g4'
        })



