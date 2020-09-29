import uuid
import os
from installed_clients.WorkspaceClient import Workspace
import shutil 


class VariationReport:

    def __init__(self, config):
        self.ws_url = config['ws_url']
        self.scratch = config['scratch']
        pass

    def create_variation_report (self, params):
        '''
        Create a table report with
        contig_id, length, number_variation, density/mb
        :param variation_ref:
        '''
        ws = Workspace(self.ws_url)

        subset = ws.get_object_subset([{
            'included': ['/numgenotypes', 'numvariants'],
            'ref': params['variation_ref']
        }])

        numgenotypes = subset[0]['data']['numgenotypes']
        numvariants = subset[0]['data']['numvariants']

        variation_table = """
        <table>
           <thead>
               <tr>
                   <td>Number of strains/genotypes</td>
                   <td> ##numgenotypes##</td>
               </tr>
            </thead>
            <tbody>
                <tr>
                    <td>Number of variants</td>
                    <td>##numvariants##</td>
                </tr>
            </tbody>
        </table>
        """
        variation_table = variation_table.replace("##numgenotypes##", str(numgenotypes))
        variation_table = variation_table.replace("##numvariants##", str(numvariants))

        session = str(uuid.uuid4())
        htmlreport_dir = (os.path.join(self.scratch, session))
        os.mkdir(htmlreport_dir)
        index_html_path = os.path.join(htmlreport_dir, "index.html")
        with open(index_html_path, "w") as f:
            f.write(variation_table)
        return (htmlreport_dir)