Channel
    .fromSRA(params.accession, apiKey: params.ncbi_api_key)
    .view()
    .set { SRA_ch }

Channel
    .of(3, 7)
    .view()
    .set { human_chrm_ch }

process chrm {

    label 'chrm'
    echo true

    input:
    val(chr) from human_chrm_ch

    output:
    file "human_genome.fa" into human_genome_ch

    script:
    """
    wget -O "$chr".fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome."$chr".fa.gz
    gunzip -c *.fa.gz > human_genome.fa
    """

}

process index {

    label 'index'
    echo true

    input:
    file "human_genome.fa" from human_genome_ch

    output:
    file 'SAindex' into index_ch

    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeFastaFiles 'human_genome.fa' --genomeSAindexNbases 12
    """
}