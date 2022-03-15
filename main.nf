c_green = "\033[0;32m";
c_blue='\033[0;34m'
c_purple = "\033[0;35m";
c_red =  "\033[0;31m";
b_green = "\033[1;32m";
c_reset = "\033[0m";

Channel
    .fromSRA(params.accession, apiKey: params.ncbi_api_key)
    .into { SRA_ch; SRA_ch2 }

Channel
    .fromPath(params.genome)
    .set { genome_ch }
process index {

    label 'index'
    echo true

    input:
    file 'genome' from genome_ch

    output:
    tuple path(genome), file("*.{amb,ann,bwt,pac,sa}") into bwa_index_ch

    script:
    """
    bwa index -a bwtsw $genome
    """
}

process align {

    label "align"
    echo true

    input :
    tuple val(sample_id), path(reads) from SRA_ch
    tuple path(genome), file("*.sa") from bwa_index_ch

    output:
    file("*.bam") into (align_out_ch, align_out_ch2, align_out_ch3)

    script:
    """
    bwa mem $genome $reads | samtools sort | samtools view -b - > ${sample_id}.bam
    """

}

process mapping {

    input:
    file(align) from align_out_ch

    output:
    file("*.bai") into mapping_ch

    script:
    """
    samtools index $align
    """

}

process annotation {

    output:
    file("*.gtf") into annotation_ch

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
    gunzip *.gff.gz
    gffread -E *.gff -T -o- | more > *.gtf
    """
}

Channel
    .of(0, 1, 2)
    .set { strand_ch }

process counting {

    input:
    file(bam) from align_out_ch2
    file(bai) from mapping_ch
    file(annot) from annotation_ch
    val(strand) from strand_ch

    output:
    file("*.counts") into counting_ch

    script:
    """
    featureCounts -g gene_id -s $strand -a $annot -o *.counts $bam
    """

}

process fastqc {

    input:
    file(bam) from align_out_ch3
    tuple val(sample_id), path(reads) from SRA_ch2

    output:
    file("*.fastqc.html") into fastqc_ch

    script:
    """
    fastqc $bam > *.fastqc.html
    fastqc $reads > *.fastqc.html
    """

}

process deseq {

    input:
    file(count) from counting_ch

    output:
    file("*.pdf") into pdf_ch
    file("*.txt") into txt_ch

    script:
    """
    "script_analyse_deseq.R" > *.pdf, *.txt
    """

}

// END processes
workflow.onComplete {
    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es)${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount}${c_reset}-"
        log.info "-${c_green}Number of successfully run process(es) : ${workflow.stats.succeedCount}${c_reset}-"
    }
    if (workflow.success) {
        log.info "-${c_purple}[${params.pipeline_name}]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        log.info "-${c_purple}[${params.pipeline_name}]${c_red} Pipeline completed with errors${c_reset}-"
    }
}