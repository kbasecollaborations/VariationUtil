import os
from installed_clients.DataFileUtilClient import DataFileUtil
from VariationUtil.Util.SampleServiceUtil import SampleServiceUtil
import logging

class StrainInfo:
    def __init__(self, config):
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.sampleservice_util = SampleServiceUtil(config)

    def _sampleset_to_strain_info(self, sample_set_ref, vcf_strain_ids):
        '''
        :param sample_set_ref:
        :param vcf_strain_ids:
        :return: StrainInfo
        order of StrainInfo should be same as order of vcf_strain_ids
        '''
        sample_set = self.dfu.get_objects({"object_refs": [sample_set_ref]})['data'][0]['data']
        samples = sample_set['samples']
        sample_dict = {}
        for sample in samples:
            name = sample['name']
            sample_dict[name] = {
                "name":name,
                "sample_id": sample['id'],
                "version": sample['version']
            }
        StrainInfo = []
        missing_strains = []
        duplicated_strains = []
        seen_strain = {}
        for strain in vcf_strain_ids:
            if strain in seen_strain:
                duplicated_strains.append(strain)
            else:
                seen_strain[strain]=1

            if strain not in sample_dict:
                missing_strains.append(strain)
            else:
                StrainInfo.append(sample_dict[strain])

        dup_strains = ", ".join (duplicated_strains)
        if duplicated_strains:
            raise ValueError(f'duplicated strain ids need to be fixed in vcf file - {dup_strains}')
        if missing_strains:
            strains_not_found = ", ". join (missing_strains)
            raise ValueError (f'Missing strains from sample set {strains_not_found}')

        return (StrainInfo)

    def _sample_set_to_attribute_mapping(self, axis_ids, sample_set_ref, obj_name, ws_id):
        am_data = self.sampleservice_util.sample_set_to_attribute_mapping(sample_set_ref)
        unmatched_ids = set(axis_ids) - set(am_data['instances'].keys())
        if unmatched_ids:
            name = "Column"
            raise ValueError(f"The following {name} IDs from the uploaded matrix do not match "
                             f"the supplied {name} attribute mapping: {', '.join(unmatched_ids)}"
                             f"\nPlease verify the input data or upload an excel file with a"
                             f"{name} mapping tab.")

        logging.info('start saving AttributeMapping object: {}'.format(obj_name))
        info = self.dfu.save_objects({
            "id": ws_id,
            "objects": [{
                "type": "KBaseExperiments.AttributeMapping",
                "data": am_data,
                "name": obj_name
            }]
        })[0]
        sample_attribute_ref = str(info[6]) + "/" + str(info[0]) + "/" + str(info[4])
        return (sample_attribute_ref)

    def sample_strain_info(self, params):
        vcf_strain_ids = params["vcf_strain_ids"]
        sample_set_ref = params["sample_set_ref"]
        ws_id = params["ws_id"]
        obj_name = params["sample_attribute_name"]

        sample_attribute_ref = self._sample_set_to_attribute_mapping(vcf_strain_ids, sample_set_ref, obj_name, ws_id)
        strains = self._sampleset_to_strain_info (sample_set_ref, vcf_strain_ids)
        return (sample_attribute_ref, strains)


