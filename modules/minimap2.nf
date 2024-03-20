
process minimap_alignment {
    label 'cpu_low'  // I keep it at 8 threads
    label 'mem_high' // not used 
    label 'time_mid' // not used
    label 'minimap' // not used as I used a local version

    input:  // inputs are two things: genomeref and reads (fastq)
    path genomeref
    path reads

    output: // outputs are the sam file and the logs
    path "mapped.sam", emit: mapped_sam
    val "minimap2", emit: aligner
    path "minimap2.command.log", emit: runlog
    
    // As discussed in the config file, here I am using the local version of minimap2 (not the right thing to do, but I just don't  want to pull more docker images locally, space limited :D)
    script: //minimap2 will allign the fastq file (reads) to the reference genome (genomeref) and will store the output in sam format (mapped.sam)
    """
    /home/ahmed/Documents/workstation_bioinf/minimap2-2.26_x64-linux/minimap2 -ax ${params.minimap2_model} -k 17 -t $task.cpus -L --secondary=no --MD --cap-kalloc 1g -K 10g $genomeref $reads > mapped.sam
    cp .command.log minimap2.command.log
    """
}