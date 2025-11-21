## Short read analysis pipeline

This is the nextflow pipeline designed to extract barcodes from short reads, optionally cluster them, and do some basic analysis.
The pipeline expects paired read data by default, but can also be configured to run on single-end sequencing data.

### Usage

Run :runner: this pipeline with default parameters:
```
nextflow run Pioneer-Research-Labs/short-read-pipeline
```

Provide a sample sheet:
```
nextflow run Pioneer-Research-Labs/short-read-pipeline --samplesheet samples.csv
```

The sample sheet is a three column csv file as follows:

* id: sample id that will be attached to output data
* r1: S3 bucket location of the sequencing file for read 1
* r2 (optional): S3 bucket location of the sequencing file for read 2 if paired end reads
* construct: construct .dna file from SnapGene with barcode flanking regions annotated (BARCODEUP, BARCODEDN)

### Notes On Parameter Selection
For most analyses, you will want to run barcode error correction and define a barcode cutoff. Both of these
values default to 2, which removes all singletons and is currently used for selection experiment analysis.  This is desirable to clean up noise when you are looking
at a high depth sample, but note this may not be appropriate for analyses where you are not sequencing most 
barcodes very deeply (e.g. library evenness checks or low depth sequencing). In this case, you may wish to retain all barcodes by setting
the cutoff to 1 (no cutoff). **If performing error correction it is important to also set the min_centroid to 1 or else no singletons will pass error correction**

### Outputs

Outputs are sorted into a folder per {id}, or whichever name was specified when the pipeline was run and will have at least the following files:
  - `read_stats.csv`: per sample statistics about number of reads
  - `{id}_barcodes_freq_{cutoff}.csv`: Barcode counts with frequencies and filtered for barcode count >= cutoff.
  - `{id}_uniq_barcodes_{cutoff}.csv`: Lists number of unique barcodes per sample per cutoff selected.

The results directory will also contain a `report.ipynb` file which is a Jupyter Notebook with some basic analysis and plots.

## All options

```
    ▗▄▄▖▗▄▄▄▖ ▗▄▖ ▗▖  ▗▖▗▄▄▄▖▗▄▄▄▖▗▄▄▖     ▗▄▄▖▗▄▄▄▖▗▄▄▖ ▗▄▄▄▖▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖ ▗▄▄▖
    ▐▌ ▐▌ █  ▐▌ ▐▌▐▛▚▖▐▌▐▌   ▐▌   ▐▌ ▐▌    ▐▌ ▐▌ █  ▐▌ ▐▌▐▌   ▐▌     █  ▐▛▚▖▐▌▐▌   ▐▌   
    ▐▛▀▘  █  ▐▌ ▐▌▐▌ ▝▜▌▐▛▀▀▘▐▛▀▀▘▐▛▀▚▖    ▐▛▀▘  █  ▐▛▀▘ ▐▛▀▀▘▐▌     █  ▐▌ ▝▜▌▐▛▀▀▘ ▝▀▚▖
    ▐▌  ▗▄█▄▖▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▙▄▄▖▐▌ ▐▌    ▐▌  ▗▄█▄▖▐▌   ▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▐▙▄▄▖▗▄▄▞▘

    Short Read Processing Pipeline          
    

    Usage: nextflow run Pioneer-Research-Labs/short-read-pipeline  [options]

    Options:
    ---------

    General:
    --outdir <path>                Output directory (default: "timestamp_samplesheet")
    --samplesheet <path>           Path to the samplesheet CSV file (default: "samples.csv")
    --barcode_cutoff <integer>     Barcode threshold where count >= cutoff for frequency calculations (default: 2)
    --paired_reads <boolean>       Whether the input reads are paired. (default: true). If false only column r1 will be used from sample sheet

    Barcode searching:
    --error_rate <float>           Error rate for barcode searching (default: 0.1)
    --min_overlap <int>            Minimum overlap for barcode searching (default: 3)
    --min_bc_len <int>             Minimum barcode length for filtering (default: 20)
    --max_bc_len <int>             Maximum barcode length for filtering (default: 60)

    Barcode correction:
    --correct                      Enable barcode correction (default: true)
    --min_centroid <int>           Minimum centroid for barcode correction (default: 2)
    --correct_error_rate <float>   Error rate for barcode correction (default: 0.1)
    --max_edits <int>              Maximum edits for barcode correction (default: 3)

    Resources:
    --cores <int>                  Number of CPU cores to use (default: 32)
    --big_mem <string>             Memory allocation for big memory processes (default: "128 GB")
    --correct_mem <string>         Memory allocation for barcode correction processes (default: "128 GB")

    Profiles:
    -profile standard              Run pipeline locally with Docker
    -profile awsbatch              Run pipeline on AWS Batch

```

## Developing

### Updating the container

If you need to modify or update the Docker container that is used to run the pipeline then you'll need to rebuild the container and push to the repo.

#### Authenticate Docker with Pioneer's ECR (local development)

After setting up aws cli authentication (use `aws sso` to do this) you can run the following:

```
aws ecr get-login-password --region us-west-1 | docker login --username AWS --password-stdin 116981805068.dkr.ecr.us-west-1.amazonaws.com/pioneer/pipelines
```

#### Building

Assumes that you have ssh keys setup for Github.  If it's not added already make sure your github ssh key is added to the ssh-agent.  This is required to install Pioneer's private python package with internal functions used by the pipeline.

```
eval `ssh-agent`
ssh-add ~/.ssh/my_key
```

Build with secure ssh mount

```
docker build --ssh default=$SSH_AUTH_SOCK . -t 116981805068.dkr.ecr.us-west-1.amazonaws.com/pioneer/pipelines:short-read-processing --platform linux/amd64
```

Push to private repo

```
docker push 116981805068.dkr.ecr.us-west-1.amazonaws.com/pioneer/pipelines:short-read-processing
```
