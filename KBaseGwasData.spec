module KBaseGwasData {


/* 
    KBase style object reference X/Y/Z to a KBaseExperiments.AttributeMapping reference
    @id ws KBaseExperiments.AttributeMapping
*/
typedef string attribute_ref;

/* 
    KBase style object reference X/Y/Z to a KBaseGenomeAnnotations.Assembly reference
    @id ws KBaseGenomeAnnotations.Assembly
*/
typedef string assembly_ref;

/* 
    KBase style object reference X/Y/Z to a KBaseGenomes.Genome reference
    @id ws KBaseGenomes.Genome
*/
typedef string genome_ref;

/* 
    handle reference for vcf file storage
    @id handle
*/
typedef string vcf_handle_ref;

/*
    handle reference for vcf index file storage
    @id handle
*/
typedef string vcf_index_handle_ref;

/*
    handle reference
    @id handle
*/

typedef string handle_ref;


/*
    Structure to show linked shock id
*/

typedef structure {
    string file_name;
    handle_ref hid;
    string id;
    string remote_md5;
    string type;
    string url;
} LinkedFile;



/* 
    Contig variation data
    contig_id - contig identifier
    totalvariants - total number of variants in each contig
    length - length of contig from assembly used in the alignment
*/


typedef structure {
    string contig_id;
    int totalvariants;
    int length;
} contig_info;

/*
    A structure to store the header information from VCF eg.
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">

    VCF format has optional fields and this structure will mirror
    the definitions in the header field of VCF
    This can help in downstream filtering apps

    ID - ID of the field eg. GQ
    Category - One of FORMAT, INFO or FILTER
    Description - Basic description like Genotype quality for GQ
    Number - Is this field a number 1 or 0
    Type - Float, int, string
    @optional Category Number Type
*/
typedef structure {
    string ID;
    string Category;
    string Description;
    string Number;
    string Type;
} headerinfo;



/*
    Variation object data structure
    num_genotypes - number of total genotypes within variant file
    num_variants - number of total variants within variant file
    contigs - list of contig ids and variant information
    samples - list of ids in VCF file for which variation information is stored
    sample_attribute_ref - Metadata for the samples is defined in this attribute mapping
    genome_ref - KBase reference to genome workspace object
    assembly_ref - KBase reference to assemebly workspace object
    vcf_handle_ref - Handle reference to bgzipped VCF file
    vcf_index_handle_ref - Handle reference to tabix indexed VCF index file
    genomic_indexes - List to store Linked index files for gff and fasta
    header - header from VCF
    @optional genome_ref sample_attribute_ref genomic_indexes header
*/
typedef structure {  
    int numgenotypes;
    int numvariants;
    list<contig_info> contigs;
    list <string> samples;
    attribute_ref sample_attribute_ref;
    genome_ref genome_ref;
    assembly_ref assembly_ref;
    vcf_handle_ref vcf_handle_ref;
    vcf_index_handle_ref vcf_index_handle_ref;
    LinkedFile vcf_handle;
    LinkedFile vcf_index_handle;
    list <LinkedFile> genomic_indexes;
    list <headerinfo> header;
} Variations;


/* 
  KBase style object reference X/Y/Z to a KBaseGwasData.Variations reference
    @id ws KBaseGwasData.Variations
*/
typedef string variation_ref;

/* 
  KBase style object reference X/Y/Z to a KBaseMatrices.TraitMatrix reference
    @id ws KBaseMatrices.TraitMatrix
*/
typedef string trait_ref;

/* 
  KBase style object reference X/Y/Z to a KBaseGenomes.Genome reference
    @id ws KBaseGenomes.Genome
*/
typedef string genome_ref;

/*
A boolean. 0 = false, other = true.
*/
typedef int bool;

/*
  SNP to Gene relationship
     gene_id - gene close to the snp
     distance_from_snp - minimum distance of gene from the snp
     snp_relative_location - inside_gene, 5', 3'
*/

typedef structure {
  string gene_id;
  int position;
  list<string> ontology_id;
} gene_info;

typedef structure {
  list<gene_info> genes;
} gene_list;

/* 
  A tuple containing revelant analytical results from association algorithms
*/

typedef list<tuple<string contig_id, string variant_id, int position, float pvalue, float pve>> association_results;

/*
  GWAS Association details
    traits - trait matricies for association study
    association_results - results from various association studies performed with sibling traits and parent variations
*/

typedef structure {
  string traits;
  association_results association_results;
} association_detail;

/*
  GWAS Associations
    description - text describing the group of association studies
    trait_ref - reference to trait_ref
    variation_id - reference to VCF file data
*/

typedef structure {
  string description;
  trait_ref trait_ref;
  variation_ref variation_id;
  list<association_detail> association_details;
} Associations;

};



