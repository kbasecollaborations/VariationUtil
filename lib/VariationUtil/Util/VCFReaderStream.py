import gzip
import json
import os
import re


class VCFReaderStream (list):
    def __init__(self, vcf_filepath):
        self.vcf_filepath = vcf_filepath
        # Limits in byte / Limits in number of samples
        # TODO: Move this filesize limit to a different location
        self.vcf_filesize_limit = 1000000
        self.chr = dict()
        
    def __iter__(self):
        return self.createGenerator()

    def is_file_ok_for_populating_genos(self):
        '''
        The workspace object size limit creates a problem for large vcf
        with too many samples. This function checks whether file size limit is ok.
        Single sample vcf may be ok. So that is always true
        TODO: We may have to change the function depending on the performance and errors.

        '''

        reader = gzip.open(self.vcf_filepath, "rt")
        for record in reader:
            if record.startswith('#CHROM'):
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, *SAMPLES = record.split("\t")
                if len(SAMPLES)==1:
                    return True
        if os.stat(self.vcf_filepath).st_size > self.vcf_filesize_limit:
            return False
        else:
            return True






    def parse_annotation(self, info_string):
        """
        parses annotation string into a structure
        [gene_id, transcript_id,
        :param ann_string:
        :return:
        """
        #
        # Fields are delimited by ;
        # annotation field starts with ANN=
        # Each allele-effect in annotation field is separated by ,
        annotation_info = list()
        ANNOT = {
            "synonymous_variant":1,
            "missense_variant":1,
            "frameshift_variant":1,
            "stop_gained":1,
            "stop_lost":1
        }

        info = info_string.split(";")

        ann_string = None
        for j in info:
            if j.startswith("ANN="):
                ann_string = j

        if ann_string is None:
            return None
       # TODO: Add more variant effects that
       # TODO: may not be affecting the protein coding region

        allele_effect_list = ann_string.split(",")
        for a in allele_effect_list:
            eff = a.split("|")
            allele = eff[0].replace("ANN=","")
            annot = eff[1]
            if annot not in ANNOT:
                continue
            gene_id = eff[3]
            transcript_id = eff[6]
            base = eff[9]
            prot = eff[10]
            annotation_info.append([allele, annot, gene_id, transcript_id, base, prot])

        if annotation_info:
            return annotation_info
        else:
            return None


    def createGenerator(self):

      populate_genos = self.is_file_ok_for_populating_genos()
      reader = gzip.open(self.vcf_filepath, "rt")
      for record in reader:
         if record[0]=='#':
             continue

         CHROM, POS, ID, REF, ALT, QUAL , FILTER, INFO, FORMAT , *GENOS = record.rstrip().split("\t")
         alleles = ALT.split(",")
         annotation = self.parse_annotation(INFO)
         v = {"var":[CHROM,POS,REF,alleles]}

         if annotation is not None:
             v['annot'] = annotation

         if populate_genos:
             v['geno'] = [re.sub(r':.*', '', i) for i in GENOS]

         yield v




