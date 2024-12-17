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
