#!/usr/bin/env nextflow

params.outdir = "results"
params.samplesheet = "samples.csv"


// barcode searching parameters
params.error_rate = 0.1
params.min_overlap = 3

// barcode searching parameters
params.min_bc_len = 20
params.max_bc_len = 60



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
    flanking = get_flanks(samples)

    // Quality filtering and merging pairs
    reads = filter_and_merge(samples) 
        | rename_reads
    
    (barcodes, report) = extract_barcodes(reads, flanking) 
    counts = barcodes    
        | filter_barcodes
        | barcode_counts

    // report
    template = channel.fromPath("${projectDir}/assets/report_template.ipynb") 
    
    preparereport(template, counts)
}

process get_flanks {
    publishDir("$params.outdir/$meta.id")
    tag("$meta.id")

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("flanking.fasta")

    script:
    """
    get_flanking.py $meta.construct
    """
}


process read_stats {

    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "read_stats.tsv"

    script:
    """
    seqkit stats -T $r1 $r2 > read_stats.tsv
    """

}

process filter_and_merge {

    cpus 4
    
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(r1), path(r2) 
    
    output:
    tuple val(meta), path("merged_reads.fastq")


    script:
    """
    fastp -i $r1 -I $r2 --correction -m --merged_out merged_reads.fastq \
            --include_unmerged -w $task.cpus \
            --json fastp_report.json
    """
}

process rename_reads {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("cleaned_reads.fastq")

    script:
    """
    seqkit replace -p .+ -r "read_{nr}" $reads > cleaned_reads.fastq
    """

}


process extract_barcodes {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(flanking)

    output:
    tuple val(meta), path("barcodes.fasta")
    path "cutadapt_report.json"

    script:
    """
    cutadapt \
        -g \$(bc_template.py $flanking cutadapt) \
        --discard-untrimmed \
        --revcomp \
        -e $params.error_rate \
        -O $params.min_overlap \
        -o barcodes.fasta \
        --json cutadapt_report.json \
        $reads
    """
}

process filter_barcodes {
    publishDir "$params.outdir/$meta.id"
    tag("$meta.id")

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path("barcodes_filtered.fasta")

    script:
    """
    seqkit seq --min-len $params.min_bc_len --max-len $params.max_bc_len \
        $barcodes > barcodes_filtered.fasta
    """
}


process barcode_counts {

    publishDir("$params.outdir/$meta.id")
    tag("$meta.id")

    input:
    tuple val(meta), path(barcodes)

    output:
    path 'barcode_counts.tsv'

    script:
    """
     seqkit fx2tab -i $barcodes | cut -f2 | sort | uniq -c | \
        awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}


process preparereport {

    publishDir("$params.outdir")
    tag 'Preparing report'

    input:
    path report
    path barcode_counts

    output:
    path 'report.ipynb'

    script:
    """
    cp $report 'report.ipynb'
    """
}
