## Short read analysis pipeline

This is the nextflow pipeline designed to extract barcodes from short reads, optionally cluster them, and do some basic analysis.

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

* Column 1: sample id
* Column 2: construct .dna file from SnapGene with barcode flanking regions annotated (BARCODEUP, BARCODEDN)
* Column 3: S3 bucket location of the sequencing file

### Outputs

Output is the `results` folder by default, or whichever name was specified when the pipeline was run and will have the following files:
  - `all_stats.csv`: per sample statistics about number of reads and barcodes
  - `all_barcodes_freq_{cutoff}.csv`: Barcode counts with frequencies and filtered for barcode count > cutoff.  There will be one file per cutoff selected in the config.
  - `all_uniq_barcodes_{cutoff}.csv`: Lists number of unique barcodes per sample per cutoff selected.

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
    --big_mem <string>             Memory allocation for big memory processes (default: "16 GB")
    --correct_mem <string>         Memory allocation for barcode correction processes (default: "16 GB")

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
