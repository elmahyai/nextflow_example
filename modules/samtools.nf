// We can delete all non used code

process sam_to_sorted_bam {
    label 'cpu_mid' // not used
    label 'mem_high' // not used
    label 'time_mid' // not used
    label 'samtools' // not used

    // define main directory to store the output of the process
    publishDir path: "${params.outdir}/${params.sampleid}/${task.process.replaceAll(':', '_')}/", mode: 'copy'
    
    input:  // the input three things: sam file, genomeref and minimap2log (just to copy it to the output directory)
    path mapped_sam
    path genomeref
    path minimap2log

    output: // we output the bam file and its index and stats files
    path "${params.sampleid}_sorted.bam", emit: sorted_bam
    path "${params.sampleid}_sorted.bam.bai", emit: bam_index
    path "*stats"
    path minimap2log, emit: minimap2log


    // Sort: Sort alignments
    // flagstat:  Does a full pass through the input file to calculate and print statistics to stdout. 
    // idxstats: Retrieve and print stats in the index file corresponding to the input file.
    // stats: samtools stats collects statistics from BAM files and outputs in a text format. 
    script:
 
    """
    samtools sort  \
        --write-index \
        -o ${params.sampleid}_sorted.bam##idx##${params.sampleid}_sorted.bam.bai \
        --reference $genomeref \
        -T sorttmp_${params.sampleid}_sorted \
        $mapped_sam
    samtools flagstat ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.flagstats
    samtools idxstats -@ 4 ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.idxstats
    samtools stats -@ 4 ${params.sampleid}_sorted.bam > ${params.sampleid}_sorted.bam.stats
    """

}



// not used in the workflow
process get_haplotype_readids {

    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'samtools'

    publishDir path: "${params.outdir}/${params.sampleid}/${task.process.replaceAll(':', '_')}/", mode: 'copy'

    input:
    path haplotagged_bam

    output:
    path "hap1ids", emit: hap1ids
    path "hap2ids", emit: hap2ids
    path "hap0ids"

    script:
    """
    samtools view -@ $task.cpus -d "HP:0" $haplotagged_bam | cut -f 1 | shuf > hap0ids
    split -n l/2 hap0ids
    mv xaa hap1ids
    mv xab hap2ids
    samtools view -@ $task.cpus -d "HP:1" $haplotagged_bam | cut -f 1 >> hap1ids
    samtools view -@ $task.cpus -d "HP:2" $haplotagged_bam | cut -f 1 >> hap2ids
    """

}



// not used in the workflow
process index_bam {

    label 'cpu_low'
    label 'mem_low'
    label 'time_low'
    label 'samtools'

    input:
    path mapped_bam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bam_index

    """
    samtools index -b -@ $task.cpus $mapped_bam
    """
}






