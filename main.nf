c_green = "\033[0;32m";
c_blue='\033[0;34m'
c_purple = "\033[0;35m";
c_red =  "\033[0;31m";
b_green = "\033[1;32m";
c_reset = "\033[0m";

Channel
    .fromSRA(params.accession, apiKey: params.ncbi_api_key)
    .set { SRA_ch }

Channel
    .fromPath(params.genome)
    .set { genome_ch }

process index {

    label 'index'
    echo true

    input:
    file 'genome' from genome_ch

    output:
    tuple path(genome), file("*.{amb,ann,bwt,pac,sa}") into (bwa_index_ch, bwa_index_ch2)

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
    file("*.bam") into align_out_ch

    script:
    """
    bwa mem $genome $reads | samtools sort | samtools view -b - > ${sample_id}.bam
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