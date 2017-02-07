#!/usr/bin/env nextflow

types = Channel.from(["DEL", "DUP", "INV", "TRA", "INS"])
bam_list = Channel.fromFilePairs("test_data/*{.bam,.bam.bai}", flat: true)
reference = Channel.fromPath("reference/WS245/WS245.fa")

bam_crossed = bam_list.spread(types).spread(reference).into { bam_cross1; bam_cross2 }

bam_cross1.println()

process run_delly {

    echo true

    input:
        set val("SM"), file("in.bam"), file("in.bam.bai"), val("type"), file("ref.fa.gz") from bam_cross2


    """
        echo ${SM}
        delly call -g ref.fa.gz -t ${type} in.bam
    """
}

