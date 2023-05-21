process fastq_dump {
    container 'ncbi/sra-tools:3.0.1'

    shell '/bin/sh'

    publishDir 'results/fastq_dump', mode: 'copy'

    input:
    val sra_accession

    output:
    tuple val(sra_accession), path("${sra_accession}_*.fastq.gz")

    script:
    """
    fastq-dump --split-3 --skip-technical --gzip ${sra_accession}
    """
}