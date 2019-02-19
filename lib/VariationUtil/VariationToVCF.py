import os
import shutil

from installed_clients.DataFileUtilClient import DataFileUtil

class VariationToVCF:
    def __init__(self, callback_url, scratch):
        self.scratch = scratch
        self.dfu = DataFileUtil(callback_url)

    def export_as_vcf(self, ctx, params):
        if 'input_var_ref' not in params:
            raise ValueError('Cannot export Variation- no input_var_ref field defined.')

        file = self.variation_to_vcf(ctx, {'variation_ref': params['input_var_ref']})

        export_dir = os.path.join(self.scratch, file['variation_name'])
        os.makedirs(export_dir)

        try:
            shutil.move(file['path'], os.path.join(export_dir, os.path.basename(file['path'])))
        except shutil.Error as e:
            exit(e)

        dfupkg = self.dfu.package_for_download({
             'file_path': export_dir,
             'ws_refs': [params['input_var_ref']]
        })

        return {'shock_id': dfupkg['shock_id']}

    def variation_to_vcf(self, ctx, params):
        self.validate_params(params)

        print('downloading ws object data: '+params["variation_ref"])

        variation_obj = self.dfu.get_objects({'object_refs': [params['variation_ref']]})['data'][0]
        ws_type = variation_obj['info'][2]
        obj_name = variation_obj['info'][1]

        if 'filename' in params:
            # TODO: check for file extentsion? add if not there (but maybe not because this call is internal
            output_filename = params['filename']
        else:
            output_filename = obj_name + '.vcf'

        # TODO: validate newly created vcf with vcf-validator

        output_vcf_file_path = os.path.join(self.scratch, output_filename)

        if 'KBaseGwasData.Variations' in ws_type:
            self.process_vcf(output_vcf_file_path, variation_obj['data'])
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseGwasData.Variations')

        return {'path': output_vcf_file_path, 'variation_name': obj_name}

    def process_vcf(self, output_vcf_file_path, data):
        self.dfu.shock_to_file({
            # TODO: does shock work with vcf_handle_ref? how does handle service work with vcf
            'handle_id': data['vcf_handle_ref'],
            'file_path': output_vcf_file_path,
        })

    def validate_params(self, params):
        for key in ['variation_ref']:
            if key not in params:
                raise ValueError('required "' + key + '" field was not defined')
