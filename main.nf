
params.outdir = "results"
params.samplesheet = "samples.csv"
params.mid = "CAGA"


// mapping paramaters for finding barcodes and inserts
params.kmer_size = 6
params.chaining_score = 5 
params.mismatch = 4


workflow {

    input_ch = channel.fromPath(params.samplesheet)
            .splitCsv(header:true)
            .map { row -> 
                meta = [id:row.id, construct:file(row.construct)]
                [meta, file(row.r1), file(row.r2)]
            }

    read_stats(input_ch)

    // get the flanking sequences from the .dna file
    flanking = getflanks(input_ch)

    // Quality filtering and merging pairs
    merged_reads = qc_reads(input_ch)


    



    // extract barcodes with bartender
    bartender_extract(flanking, merged_reads) 
        | barcode_counts
    
    // count unique barcodes
    //barcodecounts(tab)
    
    /*
    | bartender_cluster

    */
}

process getflanks {
    publishDir("$params.outdir/$meta.id")
    tag ("Extracting flanks for $meta.id")

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
    tag("Getting stats for sample $meta.id")

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    path "read_stats.tsv"

    script:
    """
    seqkit stats -T $r1 $r2 > read_stats.tsv
    """

}

process qc_reads {

    cpus 4
    
    publishDir "$params.outdir/$meta.id"
    tag("Getting stats for sample $meta.id")

    input:
    tuple val(meta), path(r1), path(r2) 
    
    output:
    tuple val(meta), path("merged_reads.fastq")

    script:
    """
    fastp -i $r1 -I $r2 --correction -m --merged_out merged_reads.fastq --include_unmerged -w $task.cpus
    """
}





process bartender_extract {

    publishDir "$params.outdir/$meta.id"
    tag("Bartender barcode extraction for sample $meta.id")

    input:
    tuple val(meta), path(flanking)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("bartender_extracted_barcode.txt")
    
    script:
    """
    bartender_extractor_com -f $reads -o bartender_extracted -p \$(bc_template.py $flanking $params.mid) -m 2
    """

}

process bartender_cluster {
    
    cpus 4

    publishDir "$params.outdir/$meta.id"
    tag("Bartender barcode clustering for sample $meta.id")

    input:
    tuple val(meta), path(barcodes)

    output:
    tuple val(meta), path("bt_clustered_barcode.csv"), path("bt_clustered_cluster.csv"), path("bt_clustered_quality.csv")

    script:
    """
    bartender_single_com -f $barcodes -o bt_clustered -t $task.cpus
    """


}



process barcode_counts {

    publishDir("$params.outdir/$meta.id")
    tag 'Counting unique barcodes'

    input:
    tuple val(meta), path(bc_seqs_tab)

    output:
    path 'barcode_counts.tsv'

    script:
    """
    cut -f1 -d, $bc_seqs_tab | sort | uniq -c | awk '{print \$2"\t"\$1}' > barcode_counts.tsv
    """
}



/*

process barcodetools_cluster {

    publishDir "$params.outdir/$meta.id"
    tag("Barcodetools barcode clustering for sample $meta.id") 

    input:
    tuple val(meta), path(barcodes) 

    script:
    """
    
    """

}
*/