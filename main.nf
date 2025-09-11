#!/usr/bin/env nextflow

// Define the usage string
def helpMessage() { 
        log.info """
Usage: nextflow run Pioneer-Research-Labs/short-read-pipeline  [options]

Options:
---------

General:
--outdir <path>                Output directory (default: "results")
--samplesheet <path>           Path to the samplesheet CSV file (default: "samples.csv")
--barcode_cutoff <list>        List of barcode count cutoffs. Produces one barcode output 
                                file per cutoff (default: [0, 5])

Barcode searching:
--error_rate <float>           Error rate for barcode searching (default: 0.1)
--min_overlap <int>            Minimum overlap for barcode searching (default: 3)
--min_bc_len <int>             Minimum barcode length for filtering (default: 20)
--max_bc_len <int>             Maximum barcode length for filtering (default: 60)

Barcode correction:
--correct                      Enable barcode correction (default: false)
--min_centroid <int>           Minimum centroid for barcode correction (default: 2)
--correct_error_rate <float>   Error rate for barcode correction (default: 0.1)
--max_edits <int>              Maximum edits for barcode correction (default: 3)

Resources:
--cores <int>                  Number of CPU cores to use (default: 4)
--big_mem <string>             Memory allocation for big memory processes (default: "128 GB")
--correct_mem <string>         Memory allocation for barcode correction processes (default: "128 GB")

Profiles:
-profile standard              Run pipeline locally with Docker
-profile awsbatch              Run pipeline on AWS Batch

    """
}



workflow {

    log.info """
▗▄▄▖▗▄▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄▖▗▄▄▖     ▗▄▄▖▗▄▄▄▖▗▄▄▖ ▗▄▄▄▖▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖ ▗▄▄▖
▐▌ ▐▌ █  ▐▌ ▐▌▐▛▚▖▐▌▐▌   ▐▌   ▐▌ ▐▌    ▐▌ ▐▌ █  ▐▌ ▐▌▐▌   ▐▌     █  ▐▛▚▖▐▌▐▌   ▐▌   
▐▛▀▘  █  ▐▌ ▐▌▐▌ ▝▜▌▐▛▀▀▘▐▛▀▀▘▐▛▀▚▖    ▐▛▀▘  █  ▐▛▀▘ ▐▛▀▀▘▐▌     █  ▐▌ ▝▜▌▐▛▀▀▘ ▝▀▚▖
▐▌  ▗▄█▄▖▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▙▄▄▖▐▌ ▐▌    ▐▌  ▗▄█▄▖▐▌   ▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▐▙▄▄▖▗▄▄▞▘

Short Read Processing Pipeline          
    """           

    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    samples = channel.fromPath(params.samplesheet)
            .splitCsv(header:true)
            .map { row -> 
                meta = [id:row.id]
                [meta, file(row.r1), file(row.r2)]
            }


    constructs = channel.fromPath(params.samplesheet)
            .splitCsv(header:true)
            .map { row -> 
                meta = [id:row.id,]
                [meta, file(row.construct)]
            }

    r_stats = read_stats(samples)

    // get the flanking sequences from the .dna file
    flanks = get_flanks(constructs)

    // Quality filtering and merging pairs
    filtered = filter_and_merge(samples) 

    // extract barcodes
    barcodes = extract_barcodes(filtered.merged.join(flanks.cutadapt_bc))

    // size filter for barcodes
    filtered_bc = filter_barcodes(barcodes.barcodes)

    // get stats for barcodes
    bc_stats = barcode_stats(barcodes.barcodes.join(filtered_bc))

    // get barcode counts
    counts = barcode_counts(filtered_bc)

    // combine stats
    stats_ch = bc_stats.barcodes \
        .join(bc_stats.barcodes_filtered) \
        .join(r_stats)


    cutoffs = Channel.fromList([[0,5]])

    if ( params.correct) {
        bc_correct = barcode_correct(counts)
        stats_ch = stats_ch.join(bc_correct.corrected_stats)
        freqs = bc_correct.corrected \
            | combine(cutoffs) \
            | add_freq
    } else {
        stats_ch = stats_ch.map { it + [ [] ] }
        freqs = counts \
            | combine(cutoffs) \
            | add_freq
    }

    // Aggregate stats
    // stats_ch | combine_stats
    //     | collect \
    //     | agg_stats


    // Aggregate barcode counts and unique counts
    // freqs.freq \
    //     | groupTuple \
    //     | join(freqs.uniq | groupTuple) \
    //     | agg_barcode_counts


    // report
    template = channel.fromPath("${projectDir}/assets/report_template.ipynb")
    prepare_report(template)
}

process get_flanks {
    
    tag("$meta.id")

    input:
    tuple val(meta), path(construct)

    output:
    tuple val(meta), path("flanking.gb"), emit: flanking
    tuple val(meta), path("cutadapt_bc.fasta"), emit: cutadapt_bc

    script:
    """
    get_flanking.py $construct
    bc_template.py flanking.gb cutadapt > cutadapt_bc.fasta
    """
}


process read_stats {

    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores
    memory params.big_mem

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("read_stats.tsv")

    script:
    """
    seqkit stats -T -j $task.cpus $r1 $r2 > read_stats.tsv
    """

}

process filter_and_merge {

    cpus params.cores
    memory params.big_mem

    //publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(r1), path(r2) 
    
    output:
    tuple val(meta), path("merged_reads.fastq"), emit: merged
    path "fastp_report.json", emit: fastp_json
    path "fastp_report.html", emit: fastp_html


    script:
    """
    fastp -i $r1 -I $r2 --correction -m --merged_out merged_reads.fastq \
            -w $task.cpus \
            --json fastp_report.json \
            --html fastp_report.html
    """
}

process rename_reads {
    //publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores 
    memory params.big_mem

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

    tag("$meta.id")

    cpus params.cores
    memory params.big_mem

    input:
    tuple val(meta), path(reads), path(flanking)    

    output:
    tuple val(meta), path("barcodes.fastq"), emit: barcodes
    path "cutadapt_report.json", emit: report

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
    //publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    cpus params.cores

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path("barcodes_filtered.fasta")
    
    script:
    """
    seqkit seq -j $task.cpus --min-len $params.min_bc_len --max-len $params.max_bc_len \
        $barcodes | seqkit fq2fa -j $task.cpus > barcodes_filtered.fasta
    """
}

process barcode_stats {
    
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")
    cpus params.cores

    input: 
    tuple val(meta), path(barcodes), path(barcodes_filtered)

    output:
    tuple val(meta), path("barcode_stats.tsv"), emit: barcodes
    tuple val(meta), path("barcodes_filtered_stats.tsv"), emit: barcodes_filtered

    script:
    """
    seqkit stats -T -j $task.cpus $barcodes > barcode_stats.tsv
    seqkit stats -T -j $task.cpus $barcodes_filtered > barcodes_filtered_stats.tsv
    """

}

process combine_stats {
    publishDir "$params.outdir/$meta.id",  mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(barcodes), path(barcodes_filtered), path(read_stats), path(correct_stats)

    output:
    path "${meta.id}_stats.csv"

    script:
    """
    combine_stats.py $barcodes $barcodes_filtered $read_stats $correct_stats --sample_name $meta.id
    """

}

process agg_stats {
    
    publishDir "$params.outdir",  mode: 'copy'

    input:
    path stats_files

    output:
    path "all_stats.csv"

    script:
    """ 
    csvtk concat -k $stats_files > all_stats.csv
    """

}

process barcode_counts {

    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag("$meta.id")

    cpus params.cores
    memory params.big_mem

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
    
    //publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag("$meta.id")

    memory params.correct_mem

    input:
    tuple val(meta), path(barcode_counts)

    output:
    tuple val(meta), path("barcodes_corrected.tsv"), emit: corrected
    tuple val(meta), path("barcodes_corrected_low_count.tsv"),  emit: corrected_low_count
    tuple val(meta), path("correct_stats.csv"), emit: corrected_stats

    script:
    """
    barcodetool_correct.py --path $barcode_counts --min-centroid $params.min_centroid \
        --error-rate $params.correct_error_rate --max-edits $params.max_edits
    """

}

process add_freq {
    publishDir("$params.outdir/$meta.id"),  mode: 'copy'
    tag("$meta.id")

    input:
    tuple val(meta), path(barcode_counts), val(cutoff)

    output:
    tuple val(cutoff), path("${meta.id}_${cutoff}_barcodes_freq.csv"), emit: freq
    tuple val(cutoff), path("${meta.id}_${cutoff}_uniq_barcodes.csv"), emit: uniq

    script:
    """
    add_freq.py $barcode_counts $meta.id --cutoff $cutoff
    """

}

process agg_barcode_counts {
    publishDir("$params.outdir"),  mode: 'copy'

    input:
    tuple val(cutoff), path(freq_files), path(uniq_files)

    output:
    path "all_barcodes_freq_${cutoff}.csv"
    path "all_uniq_barcodes_${cutoff}.csv"

    script:
    """
    csvtk concat -k $freq_files > all_barcodes_freq_${cutoff}.csv
    csvtk concat -k $uniq_files > all_uniq_barcodes_${cutoff}.csv
    """
}


process prepare_report {

    publishDir("$params.outdir"),  mode: 'copy'
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
