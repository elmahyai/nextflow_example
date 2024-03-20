nextflow.enable.dsl=2

// Each of the four processes is described in a separate file, with the include statement, we load all of them
// Some processes were not used, so I didn't import them e.g. the  get_haplotype_readids process from samtools.nf
include { minimap_alignment as minimap } from '/home/ahmed/Documents/bioinf_exercise/Nextflow_EXERCISE/modules/minimap2'
include { sam_to_sorted_bam as samtobam } from '/home/ahmed/Documents/bioinf_exercise/Nextflow_EXERCISE/modules/samtools'

include { deepvariant_snv_calling as deepvariant } from '/home/ahmed/Documents/bioinf_exercise/Nextflow_EXERCISE/modules/deepvariant'
include { clair3_snv_calling as clair3 } from '/home/ahmed/Documents/bioinf_exercise/Nextflow_EXERCISE/modules/clair3'



// This process allign the sequences (fastq file) against the reference genome fastq > sam then from sam > bam

workflow minimap_align_bamout {
    take:
        genomeref
        fastq
    
    main:
        // First, use the minimap to allign the fastq (reads) against the reference genome (genomeref) this outputs a sam file (mapped.sam)
        minimap( genomeref, fastq )
        // Second, use samtools to sort the sam file and generate a bam file
        samtobam( minimap.out.mapped_sam, genomeref, minimap.out.runlog )

    emit:
        bam = samtobam.out.sorted_bam
        idx = samtobam.out.bam_index
}


/*
* Main workflow
*/
workflow {
    // Load the genome reference
    genomeref = Channel.fromPath( params.genomeref, checkIfExists: true )  
    genomeref_index = Channel.fromPath( params.genomeref_index, checkIfExists: true )
    
 
    if (params.fastqs != "") // ensure that the fastq directory id defined
    {
        process_fastq = "yes"
        // Load each fastq file
		fastq = Channel.fromPath( params.fastqs + "*.fastq.gz", followLinks: true, checkIfExists: true ).collect()

        if (params.mapping == "yes") 
        {
            // Part 1: allign the sequences (fastq file) against the reference genome fastq > sam them from sam > bam
            minimap_align_bamout(genomeref,fastq)
        }
        if (params.snv_calling == "yes" ){   // I enabled this to allow variant calling using both deepvariant and clair3
            bam = minimap_align_bamout.out.bam
            bam_index = minimap_align_bamout.out.idx
            deepvariant( bam, bam_index, genomeref ) 
            clair3( bam, bam_index, genomeref, genomeref_index )    
        }
    }
    

}


