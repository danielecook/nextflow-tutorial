#!/usr/bin/env nextflow

github = Channel.from( ["Hello", "From", "GitHub"])

process sayHello {

    echo true

    input:
    val x from github

    """
        echo $x
    """
}

process download_large_NIL {

    publishDir "NIL/", mode: 'copy'

    output:
    file("ECA501.bam")
    file("ECA501.bam.bai")
    
    """
        curl https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam > ECA501.bam
        curl https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam.bai > ECA501.bam.bai
    """

}
