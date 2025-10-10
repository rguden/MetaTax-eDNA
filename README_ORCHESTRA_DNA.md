# ORCHESTRA-DNA
**Author:** Rodgee Mae Guden  
**Marine Biology Research Group, Ghent University**  
**Contact:** rodgeemae.guden@ugent.be  

---

## Overview
**ORCHESTRA-DNA** is a modular Bash workflow for metagenomic DNA sequencing analysis, integrating quality control, read trimming, taxonomic classification, and downstream summary and OTU table generation for Bacteria, Archaea, and Nematoda.  

This pipeline is designed for metagenomic datasets derived from complex benthic or environmental communities, enabling high-throughput taxonomic profiling through **Kraken2** and streamlined post-processing for domain-specific community summaries and abundance tables suitable for downstream ecological and functional analyses.

---

## Key Features
- ✔️ Quality Control with **FastQC** & **MultiQC**  
- ✔️ Adapter trimming using **Fastp**  
- ✔️ Taxonomic classification with **Kraken2** (paired-end reads)  
- ✔️ Summarization of **Bacteria**, **Archaea**, **Nematoda**, and **Unclassified** domains  
- ✔️ Domain-specific filtering of Kraken2 reports (taxid-based)  
- ✔️ Automated OTU table generation (base R; genus and species levels)  
- ✔️ Pre-sorted domain summary table integrating metadata  
- ✔️ Email notifications and detailed timestamped logging  
- ✔️ Modular design: run all or selected steps independently  

---

## Dependencies

| Tool      | Purpose                                 |
|------------|------------------------------------------|
| fastqc     | Read quality control                     |
| multiqc    | Aggregate QC reports                     |
| fastp      | Adapter trimming                         |
| kraken2    | Taxonomic classification                 |
| R          | OTU table processing and summarization   |

No additional R packages are required — the pipeline uses **base R** only.  

---

## Usage

### 1) Make the script executable
```bash
chmod +x orchestra_DNA_pipeline_v2_kraken.sh
```

### 2) Run the workflow
To execute all steps sequentially:
```bash
bash orchestra_DNA_pipeline_v2_kraken.sh
```

To run specific steps:
```bash
bash orchestra_DNA_pipeline_v2_kraken.sh -s qc,trim,postqc,kraken,summarize_domains,filter_bacteria,generate_bacteria_otu_R,filter_archaea,generate_archaea_otu_R,filter_nematoda,generate_nematoda_otu_R
```

Each step can be executed independently or in combination.

---

## Output
The workflow generates the following structured outputs:

- **Quality reports** before and after trimming (`FastQC`, `MultiQC`)  
- **Trimmed FASTQ files** (`fastp`)  
- **Kraken2 taxonomic classification reports** for each sample  
- **Domain summary tables** integrating metadata (`Bacteria`, `Archaea`, `Nematoda`, `Unclassified`)  
- **Domain-specific filtered Kraken2 reports** (per taxid)  
- **OTU tables** (genus and species levels) for each domain, generated via R  
- **Timestamped logs** and optional email notifications  

### Output Directory Structure
```
DNA_output/
├── quality_control/
│   ├── *.zip
│   ├── multiqc_report/
├── trimming/
│   ├── fastp/
│   └── fastqc/fastp/
├── community_analyses/
│   └── kraken_DNA_hybrid/
│       ├── *_report.txt
│       ├── domain_summary_with_metadata.tsv
│       ├── filtered_bacteria/
│       ├── filtered_archaea/
│       └── filtered_nematoda/
└── logs/
    ├── ORCHESTRA_DNA_pipeline_kraken_YYYYMMDD_HHMMSS.log
    ├── ORCHESTRA_DNA_pipeline_kraken_mem_YYYYMMDD_HHMMSS.log
    ├── latest.log → symlink to most recent
    └── latest_mem.log → symlink to most recent
```

---

## Notes
- Paths to input, output, and database directories are configurable at the top of the script.  
- The workflow assumes paired-end FASTQ files ending in `_R1.fastq.gz` and `_R2.fastq.gz`.  
- The Kraken2 database (`KRAKEN_DB`) should be pre-built and accessible.  
- Metadata file (`sample_condition.tsv`) must include the following columns:
  ```
  sample    condition    replicate
  ```
- The script includes **automatic timestamped logging** and **optional email alerts** for pipeline start, step completion, and failures.  

---

## Citation & License
© 2025 **Rodgee Mae Guden**  
**Marine Biology Research Group, Ghent University**  
_All rights reserved._

For permission to reuse or adapt this workflow, please contact the author.

---

## Acknowledgments
This workflow was developed to support **metagenomic community profiling** for the **ORCHESTRA** project, integrating bacterial, archaeal, and nematode analyses under a unified DNA-based taxonomic framework.

---

## Questions & Issues
If you encounter any issues, please open an **Issue** on the GitHub repository or contact the author directly:  
📧 **rodgeemae.guden@ugent.be**
