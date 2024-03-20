process clair3_snv_calling {
    label 'clair3'
    label 'time_high' // not used
    label 'cpu_low'

    // same as in samtools, store the output to ./data/Example_out_test_run/test_run/
    publishDir path: "${params.outdir}/${params.sampleid}/${task.process.replaceAll(':', '_')}/", mode: 'copy'
    
    input: // input include four things: bam file, bam index and genome reference and its index
    path sorted_bam
    path bam_index
    path genomeref
    path genomeref_index

    // similar to deepvariant, we output the indel_snv_vcf file with its index and the phased_indel_snv_vcf with its index and the 
    // - phased_bam file with its index
    output:
    path "${params.sampleid}_clair3_merge_output.vcf.gz", emit: indel_snv_vcf
    path "${params.sampleid}_clair3_merge_output.vcf.gz.tbi", emit: indel_snv_vcf_index
    path "${params.sampleid}_clair3_phased_merge_output.vcf.gz", emit: phased_indel_snv_vcf
    path "${params.sampleid}_clair3_phased_merge_output.vcf.gz.tbi", emit: phased_indel_snv_vcf_index
    path "${params.sampleid}_clair3_phased_output.bam", emit: phased_bam
    path "${params.sampleid}_clair3_phased_output.bam.bai", emit: phased_bam_index
    path "run_clair3.log"

    // All about clair3 is described well in their github page https://github.com/HKU-BAL/Clair3
    // After clair3 finished, we rename the files (by adding clair3 to the name), just to distinguish them from deepvariant
    // Similar to deepvariant, the most interesting output here is the vcf files which contain the called variants.
    script:
    """
    run_clair3.sh \
        --bam_fn=$sorted_bam \
        --ref_fn=$genomeref \
        --threads=$task.cpus \
        --platform=${params.clair3_platform} \
        --model_path=/opt/models/${params.clair3_model} \
        --output=clair3_output \
        --sample_name=${params.sampleid} \
        --remove_intermediate_dir \
        --use_whatshap_for_final_output_haplotagging \
        --enable_phasing 
    echo "Finish properly, time to move results"
    cp clair3_output/merge_output.vcf.gz ${params.sampleid}_clair3_merge_output.vcf.gz
    cp clair3_output/merge_output.vcf.gz.tbi ${params.sampleid}_clair3_merge_output.vcf.gz.tbi
    cp clair3_output/phased_merge_output.vcf.gz ${params.sampleid}_clair3_phased_merge_output.vcf.gz
    cp clair3_output/phased_merge_output.vcf.gz.tbi ${params.sampleid}_clair3_phased_merge_output.vcf.gz.tbi
    cp clair3_output/phased_output.bam ${params.sampleid}_clair3_phased_output.bam
    cp clair3_output/phased_output.bam.bai ${params.sampleid}_clair3_phased_output.bam.bai
    cp clair3_output/run_clair3.log run_clair3.log
    echo "Finish properly, time to copy results"
    """
}
