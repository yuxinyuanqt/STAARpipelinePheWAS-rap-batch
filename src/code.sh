#!/bin/bash
# staarpipeline_v0.9.7
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.

main() {

    echo "Value of outfile: '$outfile'"
    echo "Value of test_type: '$test_type'"
    echo "Value of pheno_file: '$pheno_file'"
    echo "Value of grm_file: '$grm_file'"
    echo "Value of nullobj_file: '$nullobj_file'"
    echo "Value of agds_file: '$agds_file'"
    echo "Value of annotation_name_catalog_file: '$annotation_name_catalog_file'"
    echo "Value of phenotype: '$phenotype'"
    echo "Value of pheno_id: '$pheno_id'"
    echo "Value of covariates: '$covariates'"
    echo "Value of het_vars: '$het_vars'"
    echo "Value of random_time_slope: '$random_time_slope'"
    echo "Value of user_cores: '$user_cores'"
    echo "Value of arrayid: '$arrayid'"
    echo "Value of min_mac: '$min_mac'"
    echo "Value of max_maf: '$max_maf'"
    echo "Value of min_rv_num: '$min_rv_num'"
    echo "Value of max_rv_num: '$max_rv_num'"
    echo "Value of max_rv_num_prefilter: '$max_rv_num_prefilter'"
    echo "Value of sliding_window_length: '$sliding_window_length'"
    echo "Value of qc_label_dir: '$qc_label_dir'"
    echo "Value of variant_type: '$variant_type'"
    echo "Value of geno_missing_imputation: '$geno_missing_imputation'"
    echo "Value of annotation_dir: '$annotation_dir'"
    echo "Value of use_annotation_weights: '$use_annotation_weights'"
    echo "Value of annotation_name: '$annotation_name'"
    echo "Value of p_filter_cutoff: '$p_filter_cutoff'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    if [ -n "$pheno_file" ]
    then
        dx download "$pheno_file" -o pheno_file &
    pheno_file2="pheno_file"
    else
    pheno_file2="NO_PHENO_FILE"
    fi
    if [ -n "$grm_file" ]
    then
        dx download "$grm_file" -o grm_file &
    grm_file2="grm_file"
    else
    grm_file2="NO_GRM_FILE"
    fi
    if [ -n "$nullobj_file" ]
    then
        dx download "$nullobj_file" -o nullobj_file &
    nullobj_file2="nullobj_file"
    else
    nullobj_file2="NO_NULL_OBJ"
    fi
    if [ -n "$agds_file" ]
    then
        dx download "$agds_file" -o agds_file.gds &
    agds_file2="agds_file.gds"
    else
    agds_file2="NO_AGDS_FILE"
    fi
    if [ -n "$annotation_name_catalog_file" ]
    then
        dx download "$annotation_name_catalog_file" -o annotation_name_catalog_file &
    annotation_name_catalog_file2="annotation_name_catalog_file"
    else
    annotation_name_catalog_file2="NO_ANNOTATION_NAME_CATALOG_FILE"
    fi
    if [ -z "$phenotype" ]
    then
    phenotype="NO_PHENO"
    fi
    if [ -z "$covariates" ]
    then
    covariates="NO_COVAR"
    fi
    if [ -z "$het_vars" ]
    then
    het_vars="NO_HET_VARS"
    fi
    if [ -z "$random_time_slope" ]
    then
    random_time_slope="NO_RANDOM_TIME_SLOPE"
    fi

    echo "Rscript --vanilla staarpipelinePheWAS.R $outfile $test_type $pheno_file2 $grm_file2 $nullobj_file2 $agds_file2 $annotation_name_catalog_file2 $phenotype $pheno_id $covariates $het_vars $random_time_slope $user_cores $arrayid $min_mac $max_maf $min_rv_num $max_rv_num $max_rv_num_prefilter $sliding_window_length $qc_label_dir $variant_type $geno_missing_imputation $annotation_dir $use_annotation_weights $annotation_name $p_filter_cutoff"
    dx-docker run -v /home/dnanexus/:/home/dnanexus/ -w /home/dnanexus/ yuxinyuan/staarpipeline:v0.9.9 Rscript --vanilla staarpipelinePheWAS.R $outfile $test_type $pheno_file2 $grm_file2 $nullobj_file2 $agds_file2 $annotation_name_catalog_file2 $phenotype $pheno_id $covariates $het_vars $random_time_slope $user_cores $arrayid $min_mac $max_maf $min_rv_num $max_rv_num $max_rv_num_prefilter $sliding_window_length $qc_label_dir $variant_type $geno_missing_imputation $annotation_dir $use_annotation_weights $annotation_name $p_filter_cutoff
    mkdir -p out/results
    mv ${outfile}.Rdata out/results
    dx-upload-all-outputs
}

