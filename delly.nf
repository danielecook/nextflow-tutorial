#!/usr/bin/env nextflow

types = Channel.from(["DEL", "DUP", "INV", "TRA", "INS"])
bam_list = Channel.fromFilePairs("test_data/*{.bam,.bam.bai}", flat: true)
reference = Channel.fromPath("reference/WS245/WS245.fa.gz")
bam_crossed = bam_list.into { bam_cross1; bam_cross2 }

bam_cross1.println()

process download_large_NIL {

    publishDir "test_data/", mode: 'copy'

    output:
    set val("ECA501"), file("ECA501.bam"), file("ECA501.bam.bai") into extra_bam
    
    """
        curl https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam > ECA501.bam
        curl https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam.bai > ECA501.bam.bai
    """

}


process extract_reference {

    input:
        file("ref.fa.gz") from reference

    output:
        file("reference.fa") into unzip_reference

    """
        gunzip -kfc ref.fa.gz > reference.fa
    """
}

delly_set = bam_cross2.concat(extra_bam).spread(types).spread(unzip_reference)

process run_delly {

    echo true

    tag { SM }

    input:
        set val("SM"), file("in.bam"), file("in.bam.bai"), val("type"), file("reference.fa") from delly_set

    output:
        val("SM") into sv_SM
        set file("${SM}.sv.bcf") into sv_bcf

    """
        delly call -g reference.fa -t ${type} in.bam > ${SM}.sv.bcf
    """
}


process combine {

    publishDir 'results/', mode: 'copy'

    input:
        val("SM") from sv_SM
        set file("${SM}.sv.bcf"), file("${SM}.sv.bcf.csi") from sv_bcf.toSortedList()

    output:
        file("sv.bcf")

    """
        delly merge `ls *.bcf` > sv.bcf
    """

}
