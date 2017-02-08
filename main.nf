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

    output:
    file("https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam")
    file("https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam.bai")
    
    """
        mkdir NIL
        curl https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam > NIL/ECA501.bam
        curl https://storage.googleapis.com/andersen/nextflow-tutorial/ECA501.bam.bai > NIL/ECA501.bam.bai
    """

}
