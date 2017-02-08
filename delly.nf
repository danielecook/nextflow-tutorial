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

/*
process extract_reference {

    input:
        file("ref.fa.gz") from reference

    output:
        file("reference.fa") into unzip_reference

    """
        gunzip -kfc ref.fa.gz > reference.fa
    """
}
*/

delly_set = bam_cross2.concat(extra_bam).spread(types).spread(reference)

process run_delly {

    echo true

    tag { "$SM - $type" }

    input:
        set val("SM"), file("in.bam"), file("in.bam.bai"), val("type"), file("reference.fa.gz") from delly_set

    output:
        file("${SM}.${type}.sv.bcf") into sv_bcf

    """
        touch ${SM}.${type}.sv.bcf
        delly call -g reference.fa.gz -t ${type} -o ${SM}.${type}.sv.bcf in.bam 
    """
}


process combine {

    publishDir 'results/', mode: 'copy'

    input:
        val "sv_set" from sv_bcf.toSortedList()

    output:
        file("sv.bcf")

    """
        merge_set=""
        for i in ${sv_set.join(" ")}; do
            if [ -s \${i} ]
            then
                    merge_set="\${merge_set} \${i}"
            fi
        done;

        echo \${merge_set}

        delly merge \${merge_set}
    """

}
