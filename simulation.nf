#!/usr/bin/env nextflow

types = Channel.from(["DEL", "DUP", "INV", "TRA", "INS"])
bam_list = Channel.fromFilePairs("test_data/*{.bam,.bam.bai}", flat: true)
reference = Channel.fromPath("reference/WS245/WS245.fa.gz")
bam_crossed = bam_list.into { bam_cross1; bam_cross2 }

bam_input = Channel.fromPath("test_data/ECA257.bam")
bam_input_index = Channel.fromPath("test_data/ECA257.bam.bai")

process generate_simulated_reference {

    input:
        file("reference.fa.gz") from reference

    output:
        file("insert_sites.tsv") into insert_sites
        file("simulated.fa") into simulated_uncompressed

    '''
    #!/usr/bin/env Rscript --vanilla
    library(RSVSim)
    library(Biostrings)
    library(tidyverse)

    # Read reference genome
    genome <- readDNAStringSet("reference.fa.gz")

    # Simulate Genome
    simulated_genome <- simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, bpSeqSize=6,
    seed=246, verbose=FALSE)

    insertion_sites <- metadata(simulated_genome) %>% as.data.frame() 

    readr::write_tsv(insertion_sites, "insert_sites.tsv")

    writeXStringSet(simulated_genome, "simulated.fa", compress=FALSE, format="fasta")
    '''
}

process compress_simulated {


    input:
        file("simulated.fa") from simulated_uncompressed

    output:
        file("simulated.fa.gz") into simulated_compressed

    """
        bgzip simulated.fa 
    """
}

process index_simulated {

    input:
        file("simulated.fa.gz") from simulated_compressed

    output:
        set file('simulated.fa.gz.amb'), file('simulated.fa.gz.ann'), file('simulated.fa.gz.bwt'), file('simulated.fa.gz.pac'), file('simulated.fa.gz.sa')  into sim_set

    """
    bwa index simulated.fa.gz
    """

}


process bam_to_fastq {

    input:
        file("a.bam") from bam_input
        file("a.bam.bai") from bam_input_index

    output:
        set file("a1.fq.gz"), file("a2.fq.gz") into fq_out

    """
        samtools sort -n a.bam > rg_sorted.bam
        samtools fastq -1 a1.fq.gz -2 a2.fq.gz rg_sorted.bam
    """
}

process perform_alignment {

    publishDir "simulated_bam", mode: 'copy'

    input:
        set file('simulated.fa.gz.amb'), file('simulated.fa.gz.ann'), file('simulated.fa.gz.bwt'), file('simulated.fa.gz.pac'), file('simulated.fa.gz.sa') from sim_set
        set file("a1.fq.gz"), file("a2.fq.gz") from fq_out
    output:
        set file("out.bam"), file("out.bam.bai")

    
    """
        bwa mem simulated.fa.gz a1.fq.gz a2.fq.gz | \\
        sambamba view --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --show-progress --tmpdir=. --out=out.bam /dev/stdin
        sambamba index out.bam
    """
}

