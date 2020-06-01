import uuid
import os
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
import shutil 


class htmlreportutils:

    def __init__(self):
        callback_url = os.environ['SDK_CALLBACK_URL']
        self.dfu = DataFileUtil(callback_url)
        self.report = KBaseReport(callback_url)
        pass

    def create_html_report(self, output_dir, workspace_name, objects_created):
        '''
         function for creating html report
        '''

        report_name = 'VariationReport' + str(uuid.uuid4())

        report_shock_id = self.dfu.file_to_shock({'file_path': output_dir,
                                            'pack': 'zip'})['shock_id']

        html_file = {
            'shock_id': report_shock_id,
            'name': 'index.html',
            'label': 'index.html',
            'description': 'Variation HTML report'
            }
        
        report_info = self.report.create_extended_report({
                        'objects_created': objects_created,
                        'direct_html_link_index': 0,
                        'html_links': [html_file],
                        'report_object_name': report_name,
                        'workspace_name': workspace_name
                    })
        return {
            'report_name': report_info['name'],
            'report_ref': report_info['ref']
        }
