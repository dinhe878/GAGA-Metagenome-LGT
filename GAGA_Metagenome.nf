/*
 * pipeline input parameters
 */

// ls -l | gawk 'match($0, /.* ([A-Z]+.[0-9]+).*/, a) {print a[1]}' will generate genome_id

params.genome =
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
         GAGA Metagenome - N F   P I P E L I N E
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


/*
 * define the `index` process that create a binary index
 * given the transcriptome file
 */
process index {

     input:
     path transcriptome from params.transcriptome

     output:
     path 'index' into index_ch

     script:
     """
     salmon index --threads $task.cpus -t $transcriptome -i index
     """
}

Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch }

process quantification {
    tag "$pair_id"

    input:
    path index from index_ch
    tuple pair_id, path(reads) from read_pairs_ch

    output:
    path pair_id into quant_ch

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"

    input:
    tuple sample_id, path(reads) from read_pairs2_ch

    output:
    path "fastqc_${sample_id}_logs" into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    path '*' from quant_ch.mix(fastqc_ch).collect()

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow.onComplete {
        log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
