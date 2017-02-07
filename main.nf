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