#!/usr/bin/env nextflow


workflow {

    channel.fromPath(params.samplesheet)
            .splitCsv(header:true)
            .map { row -> 
                meta = [id:row.id, construct:file(row.construct)]
                [meta, file(row.r1), file(row.r2)]
            }
            | set {samples}

    read_stats(samples)


    // get the flanking sequences from the .dna file
    (flanking, cutadapt_bc) = get_flanks(samples)

    // Quality filtering and merging pairs
    (merged, fastp_json, fastp_html) = filter_and_merge(samples) 
    
    cleaned = rename_reads(merged)
    
    (barcodes, report) = extract_barcodes(cleaned.join(cutadapt_bc))

    (filtered_bc, bc_stats, bc_filt_stats) = filter_barcodes(barcodes)
    counts = barcode_counts(filtered_bc)

    if ( params.correct) {
    correcteed = barcode_correct(counts)
    }

    // report
    template = channel.fromPath("${projectDir}/assets/report_template.ipynb") 
    
    prepare_report(template)
}

process get_flanks {
    publishDir("$params.outdir/$meta.id")
    tag("$meta.id")

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("flanking.fasta")
    tuple val(meta), path("cutadapt_bc.fasta")

    script:
    """
    get_flanking.py $meta.construct
    bc_template.py flanking.fasta cutadapt > cutadapt_bc.fasta
    """
}


process read_stats {

    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    cpus params.cores 

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "read_stats.tsv"

    script:
    """
    seqkit stats -T -j $task.cpus $r1 $r2 > read_stats.tsv
    """

}

process filter_and_merge {

    cpus params.cores
    memory '2 GB'

    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(r1), path(r2) 
    
    output:
    tuple val(meta), path("merged_reads.fastq")
    path "fastp_report.json"
    path "fastp_report.html"


    script:
    """
    fastp -i $r1 -I $r2 --correction -m --merged_out merged_reads.fastq \
            -w $task.cpus \
            --json fastp_report.json \
            --html fastp_report.html
    """
}

process rename_reads {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    cpus params.cores 

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("cleaned_reads.fastq")

    script:
    """
    seqkit replace -j $task.cpus -p .+ -r "read_{nr}" $reads > cleaned_reads.fastq
    """

}


process extract_barcodes {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")


    cpus params.cores
    mem '16 GB'

    input:
    tuple val(meta), path(reads), path(flanking)

    output:
    tuple val(meta), path("barcodes.fastq")
    path "cutadapt_report.json"

    script:
    """
    cutadapt \
        -g  file:$flanking \
        --discard-untrimmed \
        --revcomp \
        --cores $task.cpus \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o barcodes.fastq \
        --json cutadapt_report.json \
        $reads
    """
}

process filter_barcodes {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path("barcodes_filtered.fasta")
    path "barcode_stats.tsv"
    path "barcodes_filtered_stats.tsv"

    script:
    """
    seqkit stats -T -j $task.cpus $barcodes > barcode_stats.tsv
    seqkit seq -j $task.cpus --min-len $params.min_bc_len --max-len $params.max_bc_len \
        $barcodes | seqkit fq2fa -j $task.cpus > barcodes_filtered.fasta
    seqkit stats -T -j $task.cpus barcodes_filtered.fasta > barcodes_filtered_stats.tsv
    """
}


process barcode_counts {

    publishDir("$params.outdir/$meta.id")
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path('barcode_counts.tsv')

    script:
    """
     seqkit fx2tab -j $task.cpus -i $barcodes | cut -f2 | sort | uniq -c | \
        awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}


process barcode_correct {
    
    publishDir("$params.outdir/$meta.id")
    tag("$meta.id")

    memory '32 GB'

    input:
    tuple val(meta), path(barcode_counts)

    output:
    path "barcodes_corrected.tsv"

    script:
    """
    barcodetool_correct.py $barcode_counts 
    """



}

process prepare_report {

    publishDir("$params.outdir")
    tag 'Preparing report'

    input:
    path report

    output:
    path 'report.ipynb'

    script:
    """
    cp $report 'report.ipynb'
    """
}
