/*
A KBase module: VariationUtil
*/

module VariationUtil {
    /*
		TODO:
			take assembly ref as input
			take care of metagenomic use case
    */

    /* typedefs for input/output values */

    /* A boolean - 0 for false, 1 for true. @range (0, 1) */
    typedef int boolean;

    /* An X/Y/Z style reference*/
    typedef string obj_ref;

    /* KBase file path to staging files */
    typedef string filepath;

    /*
        ## funcdef save_variation_from_vcf ##

        required input params:
            genome_ref: KBaseGenomes.Genome object reference
	    
	    	*** variation input data ***
			vcf_staging_file_path: path to location data associated with samples
        	variation_object_name: output name for KBase variation object 
			
			*** sample input data ***
			sample_attribute_ref: x/y/z reference to kbase sample attribute
			
        optional params:
            NA

        output report:
            report_name
            report_ref

            HTML visualization: Manhattan plot

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
    } save_variation_input;

    typedef structure {
        string report_name;
        string report_ref;
    } save_variation_output;

    /*
        Save a variation (and trait?) object to Kbase given a reference genome, object output name,
        Variant Call Format (VCF) file, and sample attribute file.

        TODO: Viewer for Variation/Trait object?
    */

    funcdef save_variation_from_vcf(save_variation_input params)
		returns (save_variation_output report) authentication required;

	/*
        ## funcdef export_variation_as_vcf ##

        required input params:
            Variation object reference

        optional params:
            NA

        output report:
            Shock id pointing to exported vcf file
    */

	typedef structure {
	    obj_ref input_var_ref;
	} export_variation_input;

	typedef structure {
	    string shock_id;
	} export_variation_output;

	/*
	    Export KBase variation object as Variant Call Format (VCF) file
	*/

	funcdef export_variation_as_vcf(export_variation_input params)
	    returns (export_variation_output output) authentication required;

	/*
        ## funcdef get_variation_as_vcf ##

        required input params:
            Variation object reference
            output file name

        optional params:
            NA

        output report:
            path to returned vcf
            name of variation object
    */

	typedef structure {
	    obj_ref variation_ref;
	    string filename;
	} get_variation_input;

	typedef structure {
	    filepath path;
	    string variation_name;
	} get_variation_output;

	/*
	    Given a reference to a variation object, and output name: return a Variant Call Format (VCF)
	    file path and name.
    */

    funcdef get_variation_as_vcf(get_variation_input params)
	    returns (get_variation_output file) authentication required;
};
