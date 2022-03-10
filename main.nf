params.ncbi_api_key = '24e61d2d45cab93c6268b0d572d9fe82ac09'

Channel
    .from('SRR11462797')
    .view()
    .set { SRA_ch }

process fastq {

    label 'fastq'
    echo true

    input:
    val(id) from SRA_ch

    output:
    tuple val(id), file("*_{1,2}.fastq") into fastq_ch

    script:
    """
    fasterq-dump $id
    """

}

fastq_ch.view()