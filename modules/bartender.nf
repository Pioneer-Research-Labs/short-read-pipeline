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
    bartender_extractor_com -f $reads -o bartender_extracted \
        -p \$(bc_template.py $flanking bartender $params.mid) -m 2
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