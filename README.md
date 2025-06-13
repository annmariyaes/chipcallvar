# chipcallvar

**chipcallvar** is a reproducible [Nextflow](https://www.nextflow.io/) workflow to call variants (SNVs and INDELs) on ChIP-seq data.
There are 3 different variant callers in this pipeline: [macs3 callvar](https://macs3-project.github.io/MACS/docs/callvar.html), [GATK mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), [freebayes](https://github.com/freebayes/freebayes)


## 🧬 Workflow Overview

This pipeline performs the following steps:

1. **Read mapping** – `bwa-mem2`
2. **Merging technical duplicates** – `samtools`
3. **Create intervals** - `bedtools`
4. **Peak calling** – `macs3 callpeak`
5. **Variant calling** – `macs3 callvar, GATK mutect2, freebayes`
6. **Variant annotation** – Ensembl `vep`
7. **Variant filtering** – `bcftools`
8. **Downstream analysis** – `vcf2maf, maftools`
8. **MultiQC** - `fastqc`, `samtools`, `mosdepth`, `bcftools`, `ensembl vep`

<img width="1340" alt="workflow" src="https://github.com/user-attachments/assets/a1821c20-c71e-4d9f-ba12-5c5abc14fe74" />


## Usage

### 📂 Input: Sample Sheet

Each row represents a pair of fastq files (paired end). Prepare a CSV file (`samplesheet.csv`) with the following format:

```csv
patient,sample,replicate,fastq_1,fastq_2,control,control_replicate
OCI-AML3,OCI-AML3_input,1,test/hs_ChIP_OCI-AML3_rep1_Input_R1_001.fastq.gztest/hs_ChIP_OCI-AML3_rep1_Input_R2_001.fastq.gz,,
OCI-AML3,OCI-AML3_H3K27ac,1,test/hs_ChIP_OCI-AML3_rep1_H3K27ac_R1_001.fastq.gz,test/hs_ChIP_OCI-AML3_rep1_H3K27ac_R2_001.fastq.gz,OCI-AML3_input,1
OCI-AML3,OCI-AML3_H3K27ac,2,test/hs_ChIP_OCI-AML3_rep2_H3K27ac_R1_001.fastq.gz,test/hs_ChIP_OCI-AML3_rep2_H3K27ac_R2_001.fastq.gz,OCI-AML3_input,1
```

### ⚙️ Parameters File (params.yaml)
```yaml
samplesheet: './samplesheet_example.csv'
outdir: "./nf-macs3"
fasta: './reference/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta'
fai: './reference/resources-broad-hg38-v0-Homo_sapiens_assembly38.fasta.fai'
dict: './reference/resources-broad-hg38-v0-Homo_sapiens_assembly38.dict'
assembly: 'GRCh38'
genome_size: 'hs'
step: 'mapping'
tools: 'macs3,mutect2,freebayes'
email: 'example@gmail.com'
```

### 🚀 Running the Pipeline

```bash
nextflow run main.nf \
   -params-file params.yaml \
   -resume
```


## 👩‍💻 Author

Ann Mariya
[GitHub](https://github.com/annmariyaes)
[Email](annmariya.elayani@gmail.com)
