## Short read analysis pipeline

This is the nextflow pipeline designed to extract barcodes from short reads, optionally cluster them, and do some basic analysis.

### Usage

Run this pipeline with default parameters:
```
nextflow run Pioneer-Research-Labs/short-read-pipeline
```

Provide a sample sheet:
```
nextflow run Pioneer-Research-Labs/short-read-pipeline --samplesheet samples.csv
```

The sample sheet is a three column csv file as follows:
Column 1: sample id
Column 2: construct .dna file from SnapGene with barcode flanking regions annotated (BARCODEUP, BARCODEDN)
Column 3: S3 bucket location of the sequencing file

### Outputs

Output is the `results` folder, or whichever name was specified when the pipeline was run.  There will be a subfolder for each sample that has the output files.  
The results directory will also contain a `report.ipynb` file which is a Jupyter Notebook with some basic analysis and plots.

