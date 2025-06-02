# chipcallvar

**chipcallvar** is a reproducible [Nextflow](https://www.nextflow.io/) workflow for variant calling on ChIP-seq data.


## 🖼️ Workflow Diagram

<img width="1340" alt="workflow" src="https://github.com/user-attachments/assets/a1821c20-c71e-4d9f-ba12-5c5abc14fe74" />


## 🧬 Workflow Overview

This pipeline performs the following steps:

1. **Sequencing quality control** – FastQC  
2. **Read mapping** – BWA-mem2  
3. **Merging technical duplicates** – SAMtools merge  
4. **Peak calling** – MACS3 `callpeak`   
5. **Variant calling** – MACS3 `callvar`  
6. **Variant annotation** – Ensembl `vep`  
7. **Variant filtering** – BCFtools, VCF2MAF


# Usage

## 📂 Input: Sample Sheet

Each row represents a pair of fastq files (paired end). Prepare a CSV file (`samplesheet.csv`) with the following format:

```csv
patient,sample,replicate,fastq_1,fastq_2,control,control_replicate
OCI-AML3,OCI-AML3_input,1,hs_ChIP_OCI-AML3_rep1_Input_60PE_JB_S13_R1_001.fastq.gz, hs_ChIP_OCI-AML3_rep1_Input_60PE_JB_S13_R2_001.fastq.gz,,
OCI-AML3,OCI-AML3_H3K27ac,1,hs_ChIP_OCI-AML3_rep1_H3K27ac_60PE_JB_S9_R1_001.fastq.gz,hs_ChIP_OCI-AML3_rep1_H3K27ac_60PE_JB_S9_R2_001.fastq.gz,OCI-AML3_input,1
OCI-AML3,OCI-AML3_H3K27ac,2,hs_ChIP_OCI-AML3_rep2_H3K27ac_60PE_JB_S10_R1_001.fastq.gz,hs_ChIP_OCI-AML3_rep2_H3K27ac_60PE_JB_S10_R2_001.fastq.gz,OCI-AML3_input,1
```


## 🚀 Running the Pipeline

```bash
nextflow run main.nf \
   -params-file params.yaml \
   -resume
```

Adjust the `-profile` to match your environment (e.g., `singularity`, `conda`, or institutional profile).


## 👩‍🔬 Author

Annmariya 
