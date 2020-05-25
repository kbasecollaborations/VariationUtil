from VariationUtil.Util.JbrowseUtil import JbrowseUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace



class report:
    def __init__(self,  params):
        self.ws_url = params['ws_url']
        self.scratch = params['scratch']
        self.callback_url = params['callback_url']
        self.vcf_local_file_path = params['vcf_local_file_path']
        self.dfu = DataFileUtil(self.callback_url)
        self.wsc = Workspace(self.ws_url)



        pass

    def prepare_gff(self, gff_file):
        sorted_gff_cmd = "sort -k1,1 -k4,4n " + gff_file + " > " + "sorted_" + gff_file
        self.run_cmd(sorted_gff_cmd)

        zip_cmd = "bgzip " + "sorted_" + gff_file
        self.run_cmd(zip_cmd)

        index_gff_cmd = "tabix -p gff " + "sorted_" + gff_file + ".gz"
        self.run_cmd(index_gff_cmd)

        return {"gff_file_path": "sorted_" + gff_file + ".gz", "index_file_path": "sorted_" + gff_file + ".gz.tbi"}


    def prepare_report(self, var_data):
        #var_data = var_obj[1]
        jbrowse_input_params = {}
        jbrowse_input_params["ws_url"] = self.ws_url
        jbrowse_input_params['assembly_ref'] = var_data['assembly_ref']
        jbrowse_input_params['scratch'] = self.scratch
        jbrowse_input_params['vcf_filepath'] = self.vcf_local_file_path
        jbrowse_input_params['binsize'] = 10000
        jbrowse_input_params['bedGraphToBigWig'] = "/kb/deployment/bin/bedGraphToBigWig"
        jbrowse_input_params["vcf_shock_id"] = var_data['vcf_handle']['id']
        jbrowse_input_params["vcf_index_shock_id"] = var_data["vcf_index_handle"]["id"]
        jbrowse_input_params['callback_url'] = self.callback_url
        jbrowse_input_params['dfu'] = self.dfu
        jbrowse_input_params['wsc'] = self.wsc

        if 'genome_ref' in var_data:
            print ("There is genome_data")
            jbrowse_input_params['genome_ref'] = var_data['genome_ref']


        jb = JbrowseUtil()
        jbrowse_report = jb.prepare_jbrowse_report(jbrowse_input_params)
        return jbrowse_report







