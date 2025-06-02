# chipcallvar

**chipcallvar** is a reproducible [Nextflow](https://www.nextflow.io/) workflow to call variants (SNVs and INDELs) on ChIP-seq data primarily taking advantage of `macs3 callvar` variant caller.


## 🧬 Workflow Overview

This pipeline performs the following steps:

1. **Sequencing quality control** – FastQC  
2. **Read mapping** – `bwa-mem2`
3. **Merging technical duplicates** – `samtools`
4. **Peak calling** – MACS3 [callpeak](https://macs3-project.github.io/MACS/docs/callpeak.html)  
5. **Variant calling** – MACS3 [callvar](https://macs3-project.github.io/MACS/docs/callvar.html)  
6. **Variant annotation** – Ensembl `vep`
7. **Variant filtering** – BCFtools, VCF2MAF

<img width="1340" alt="workflow" src="https://github.com/user-attachments/assets/a1821c20-c71e-4d9f-ba12-5c5abc14fe74" />


## Usage

### 📂 Input: Sample Sheet

Each row represents a pair of fastq files (paired end). Prepare a CSV file (`samplesheet.csv`) with the following format:

```csv
patient,sample,replicate,fastq_1,fastq_2,control,control_replicate
OCI-AML3,OCI-AML3_input,1,ChIP_OCI-AML3_rep1_Input_R1_001.fastq.gz, ChIP_OCI-AML3_rep1_Input_R2_001.fastq.gz,,
OCI-AML3,OCI-AML3_H3K27ac,1,ChIP_OCI-AML3_rep1_H3K27ac_R1_001.fastq.gz,ChIP_OCI-AML3_rep1_H3K27ac_R2_001.fastq.gz,OCI-AML3_input,1
OCI-AML3,OCI-AML3_H3K27ac,2,ChIP_OCI-AML3_rep2_H3K27ac_R1_001.fastq.gz,ChIP_OCI-AML3_rep2_H3K27ac_R2_001.fastq.gz,OCI-AML3_input,1
```

### ⚙️ Parameters File (params.yaml)
```yaml
OUTDIR: "./nf-macs3"
email: 'example@gmail.com'
GENOME: 'WholeGenomeFasta/genome.fa'
GENOME_SIZE: "hs"
SAMPLESHEET: "./samplesheet_example.csv"
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
