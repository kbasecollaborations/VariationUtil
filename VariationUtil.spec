/*
A KBase module: VariationUtil
*/

module VariationUtil {
    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int boolean;

    /* 
	An X/Y/Z style reference
    */
    typedef string obj_ref;
    /*            
        KBase file path to staging files
    */
    typedef string filepath;

    /*
		TODO:
			take assembly ref as input
			take care of metagenomic use case

        required params:
            genome_ref: KBaseGenomes.Genome object reference
	    
	    	*** variation input data ***
			vcf_staging_file_path: path to location data associated with samples
        	variation_object_name: output name for KBase variation object 
			
			*** sample input data ***
			sample_attribute_ref: x/y/z reference to kbase sample attribute
			
        optional params:
            *** Visualization ***
            plot_maf: generate histogram of minor allele frequencies
            plot_hwe: generate histogram of Hardy-Weinberg Equilibrium p-values
    */
    
    typedef structure {
        string workspace_name;
        obj_ref genome_ref;
        filepath vcf_staging_file_path;
		string variation_object_name;
		obj_ref sample_attribute_ref;
    } import_variation_params;

    typedef structure {
        string report_name;
        string report_ref;
    } import_variation_results;


    funcdef import_variation(import_variation_params) 
        returns (import_variation_results) authentication required;
};
