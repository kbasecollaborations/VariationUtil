import os
import json
import unittest
import requests
import hashlib
from VariationUtil.Util.VCFUtils import VCFUtils

#from VariationUtil.Util.VCFToVariation import VCFToVariation



small_poplar_vcf_gz = "/kb/module/test/sample_data/small_poplar/small_poplar.vcf.gz"

params = {
    "vcf_staging_file_path": small_poplar_vcf_gz,
    "vcf_local_file_path": small_poplar_vcf_gz
}
config = {
    "scratch": "/kb/module/work"
}
def md5_sum_local_file(fname):
    md5hash = hashlib.md5()
    with open(fname, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5hash.update(chunk)
    return md5hash.hexdigest()



class TestApi(unittest.TestCase):


    @unittest.skip('x')
    def test_VCFUtils_stage_vcf(self):
        print("running test: test stage and sanitize vcf")
        vv = VCFUtils(params, config)
        result = vv.sanitize_vcf()
        print (result[0])
        print (result[1])
        print (md5_sum_local_file(result[0]))
        print (md5_sum_local_file(result[1]))
        self.assertTrue (md5_sum_local_file(result[0])=='88ba8571cd42ef2b3936d89109386ba2')
        self.assertTrue (md5_sum_local_file(result[1])=='4a1e2d1c284d59d5ed2ece54c66e5993')

    @unittest.skip('x')
    def test_VCFUtils_parse_vcf(self):
        print ("running test: test parse vcf ")
        vv = VCFUtils (params, config)
        vcf_info = vv.parse_vcf_data(params['vcf_local_file_path'])
        self.assertTrue (vcf_info['version']=="VCFv4.1")
        self.assertTrue (len(vcf_info['contigs'])==1)
        self.assertTrue (len(vcf_info['genotype_ids'])==917)
        self.assertTrue (len(vcf_info['chromosome_ids'])==1)
        self.assertTrue (len(vcf_info['header'])==2)


