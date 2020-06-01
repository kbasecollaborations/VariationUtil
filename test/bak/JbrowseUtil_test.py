import os
import json
import unittest
import requests
import hashlib

from VariationUtil.Util.JbrowseUtil import JbrowseUtil



Config = {
    "ws_url": "https://appdev.kbase.us/services/ws",
    "scratch": "/kb/module/work/tmp",
}

jbrowse_params = {
    "genome_ref": "41807/77/1",
    "assembly_ref": "41807/76/1",
    "binsize": 10000,
    "vcf_shock_id": "429099d1-30a4-45a1-a59e-bd98edd25a93",
    "vcf_index_shock_id": "aa475ffe-a3b6-4813-9532-d8c8de27282f"
}

jb = JbrowseUtil(Config)
jbrows_report = jb.prepare_jbrowse_report (jbrowse_params)
print (jbrows_report)

