#!/usr/bin/env nextflow

// Validate parameters
if (params.sra_accession == '') {
    error "You must provide an SRA accession number with --sra_accession"
}

accessions = params.sra_accession.split(',')
accessions_channel = Channel.from(accessions)

process fasterq_dump {
    container 'ncbi/sra-tools:aarch64-3.0.1'

    // Specify "/bin/sh" for this container
    shell '/bin/sh'

    input:
    val(sra_accession) from accessions_channel

    output:
    set val(sra_accession), file("${sra_accession}*.fastq") into reads

    """
    fasterq-dump ${sra_accession}
    """
}

process fastqc {
    tag "$name"
    container 'fastqc_trimmomatic:latest'
    
    input:
    set val(name), file(reads) from reads

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}