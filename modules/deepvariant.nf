process deepvariant_snv_calling {
    label 'deepvariant'
    label 'time_high' // not used
    label ( workflow.profile.contains('qsub') ? null: 'cpu_low' )

    stageInMode 'copy'
    // Same as in samtools, store the output to ./data/Example_out_test_run/test_run/
    publishDir path: "${params.outdir}/${params.sampleid}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input: // The input include three things: bam file, bam index and genome reference
    path sorted_bam
    path bam_index
    path genomeref

    // We output three things: 1) indel_snv_vcf file with its index, phased_indel_snv_vcf with its index and haplotagged_bam with its index
    // We also output stats files.
    output:
    path "deepvar_out/${params.sampleid}_deepvariant.vcf.gz", emit: indel_snv_vcf
    path "deepvar_out/${params.sampleid}_deepvariant.vcf.gz.tbi", emit: indel_snv_vcf_index
    path "deepvar_out/${params.sampleid}_deepvariant_phased.vcf.gz", emit: phased_indel_snv_vcf
    path "deepvar_out/${params.sampleid}_deepvariant_phased.vcf.gz.tbi", emit: phased_indel_snv_vcf_index

    path "deepvar_out/*_deepvariant_haplotagged.bam", emit: haplotagged_bam
    path "deepvar_out/*_deepvariant_haplotagged.bam.bai", emit: haplotagged_bam_idx
    path "deepvar_out/*stats"

    script:

    // Just define the number of threads, here I just defined them in the config file to be 8 
    def localproc = ( workflow.profile.contains('qsub') ? 36: task.cpus )

    // The pepper method and the parameters are explained on their github page
    // https://github.com/kishwarshafin/pepper/blob/r0.8/docs/misc/pepper_methods.md
    // https://github.com/kishwarshafin/pepper/blob/r0.8/docs/usage/usage_and_parameters.md

    // After  call_variant finished, we rename  the output files (just adding the deepvariant word to the name, to distinguish from clair3) 
    // We recall the samtools on the _haplotagged.bam file
    // Generally, the most intersting file here is the vcf files which include the called variants.
    """
    run_pepper_margin_deepvariant call_variant \
        -b $sorted_bam \
        -f $genomeref \
        -o deepvar_out \
        -p ${params.sampleid} \
        -s ${params.sampleid} \
        -t ${localproc} \
        --${params.deepvariant_model} \
	    --phased_output
    cp ./deepvar_out/${params.sampleid}.vcf.gz ./deepvar_out/${params.sampleid}_deepvariant.vcf.gz
    cp ./deepvar_out/${params.sampleid}.vcf.gz.tbi ./deepvar_out/${params.sampleid}_deepvariant.vcf.gz.tbi
    cp ./deepvar_out/${params.sampleid}.phased.vcf.gz ./deepvar_out/${params.sampleid}_deepvariant_phased.vcf.gz
    cp ./deepvar_out/${params.sampleid}.phased.vcf.gz.tbi ./deepvar_out/${params.sampleid}_deepvariant_phased.vcf.gz.tbi
    cp ./deepvar_out/${params.sampleid}.haplotagged.bam ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam
    samtools index -b -@ 18 ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam.bai
    samtools flagstat ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam > ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam.flagstats
    samtools idxstats ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam > ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam.idxstats
    samtools stats ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam > ./deepvar_out/${params.sampleid}_deepvariant_haplotagged.bam.stats
    """
}
