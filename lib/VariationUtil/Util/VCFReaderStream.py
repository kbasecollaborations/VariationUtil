import gzip
import json

class VCFReaderStream (list):
    def __init__(self, vcf_filepath):
        self.vcf_filepath = vcf_filepath 
        self.chr = dict()
        
    def __iter__(self):
        return self.createGenerator()

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
        #  may not be affecting the protein coding region

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
      reader = gzip.open(self.vcf_filepath, "rt")
      for record in reader:
         if record[0]=='#':
             continue
         CHROM, POS, ID, REF, ALT, QUAL , FILTER, INFO, FORMAT , *_ = record.split("\t")
         alleles = ALT.split(",")
         annotation = self.parse_annotation(INFO)
         if annotation is not None:
             i = {"var":[CHROM,POS,REF,alleles], "annot": annotation}
         else:
             i = {"var":[CHROM,POS,REF,alleles]}
         yield i




