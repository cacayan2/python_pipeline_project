# HCMV Transcriptome Analysis Pipeline
## Author: Emil Cacayan

### Table of Contents

- **[Introduction](#introduction)**
- **[Data](#data)**
- **[Usage](#usage)**
    + **[Dependencies](#dependencies)**
    + **[The `input_srr` File](#the-input_srr-file)**
    + **[Command Arguments](#command-arguments)**
- **[RNA-Seq Transcriptomic Pipeline Workflow](#rna-seq-transcriptomic-pipeline-workflow)**
- **[Disclaimer](#disclaimer)**
- **[References](#references)**

### Introduction

Transcriptomics is the study of the sum set of all RNA transcripts in a cell, organism, or other entity that manipulates genetic information. It is distinct from genomics as it provides a real-time view of the nexus between the information stored in genetic material and how it may be expressed. 

Researchers have characterized transcriptomes across many different organisms, but as they can provide a synchronous view of an organism's response to stimuli, it is useful to do transcriptomic analysis preceding and subsequent some event or as a time-series analysis. 

One of the vectors of genetic information of great interest to researchers are viruses - both in their clinical repercussions and clinical/nonclinical applications. Transcriptomic analysis of viruses produces data that allows researchers to more comprehensively study gene expression patterns, which gives the following potentially actionable insights:

1. **Characterization of viral life cycles.**
    - Over the course of infection, viruses undergo a number of changes in its genetic expression that can be tracked using transcriptomic analysis (Sudhagar et al., 2018).
    - Tracking these changes can provide researchers valuable information about potential drug targets (Mo et al., 2016).
2. **Identification of host response pathways.**
    - While this is not one of the objectives of this particular analysis, establishing the genetic expression profile of a host cell over the course of infection can help outline important infection response genes (Zhou et al., 2021).
3. **Annotation/Characterization of viral genes.**
    - Transcriptomic analysis allows for the classification of viral genes that may cause/escalate pathogenesis (Ivanov et al., 2023).
4. **Classification/prediction of disease severity.**
    - One of the hallmark results of bioinformatic analysis is identifications of biomarkers.
    - Transcriptomic analysis of host/viral transcriptomes allows researchers to identify biomarkers to elucidate disease progression/severity (Arriaga-Canon et al., 2022).
5. **Drug development and discovery.**
    - In understanding the information above, researchers can take advantage of the molecular mechanisms of viral infection to develop drugs to disrupt virus processes or interactions between the virus and its host (Zhou et al., 2021).

### Data

In this analysis, we wish to compare transcriptomes 2-days and 6-days post-infection (dpi) of patients infected with human herpesvirus 5, also known as human cytomegalovirus (HCMV), data which has been collected and published in a study by Cheng et al., 2017.

The sequence read archive data used in this analysis is used here:

- **Donor 1 (2 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896360
- **Donor 1 (6 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896363
- **Donor 3 (2 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896374
- **Donor 3 (6 dpi):** https://www.ncbi.nlm.nih.gov/sra/SRX2896375

Because this analysis was performed in a UNIX environment, the `wget` command was used to download the data, e.g.:

More information about how to correctly obtain these files to can be found in the [Usage](#usage) section under [The `input_srr` File](#the-input_srr-file).

These transcriptomes will be indexed against the HCMV transcriptome (NCBI accession [NC_006273.2](https://www.ncbi.nlm.nih.gov/sra/SRX2896375)) (Gatherer et al., 2011).

### Usage
The goal of the design of this pipeline is to make its use as user-friendly and oriented as possible.

There are only two things that you have to change:
1. The `input_srr` file.
2. The command itself.

#### Dependencies
Use of this repository requires the following software:
- [Linux/Unix](https://www.linux.org/pages/download/)
- [Python3](https://www.python.org/downloads/)
- [Biopython](https://biopython.org/wiki/Download)
- [Kallisto](https://pachterlab.github.io/kallisto/download)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [SPAdes](https://github.com/ablab/spades)
- [`R` (Software)](https://cran.r-project.org/bin/windows/base/)
- [BLAST+ suite](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [SRA Toolkit (for fasterq-dump)](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) 

#### The `input_srr` File
In the root directory, you will find a file called `input_srr`. In this file, four SRR accession numbers will already be present. In this file, just put in any number of SRR accession numbers and the pipeline will execute on all the SRR numbers (for now, this implementation works on only paired-end reads).

The SRR accession numbers already present, which are:

```
SRR5660030
SRR5660033
SRR5660044
SRR5660045
```

Are the accession numbers associated with the HCMV study outlined in the beginning of the [Data](#data) section. If you would like to use different sequence read archive numbers, simply open the SRA entry and copy-and-paste the run number into the `input_srr` file (keep in mind run numbers and SRA accession numbers are distinct from one another - the SRR number associated with the SRA accession number SRX2896360 has an associated run number of SRR5660030 - in other words, what goes into the `input_srr` file isn't the SRA accession number but run numbers).

This is important because to obtain the raw transcriptomic reads, the RNA-Seq data was obtained from the sequence read archive database using `wget` after setting the current working directory as the root directory (`/python_pipeline_project` - the links provided below are the direct download links, found by opening the SRA run browser and navigating to the data access tab, the data used here is from AWS in the SRA normalized format):

**Please note the below command is shown for the sake of demonstration - by entering your SRR numbers into the `input_srr` file, these links are generated automatically**. 
```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030 
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033 
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044 
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045 
```

The `input_srr` file will download the data associated with the sequence read archive run in its entirety. This means each time that you run the pipeline, the SRR is entirely downloaded each time. If you would like to run the pipeline on a truncated sample dataset, please see the [Command Arguments](#command-arguments) section.

#### Command Arguments
To run the script, please clone the repository:

```
git clone https://github.com/cacayan2/python_pipeline_project
```

The repository contains a test dataset (a subset of 40,000 lines corresponding to 10,000 reads) from the input fastq files (used the `head` command to do this) after manually downloading with `wget` and unpacking with `fasterq-dump`. 

Use of the test data set will allow you to test a smaller dataset without having to wait for the SRR to download (since it is already present in the repository) to allow for troubleshooting.

The command to run the script is as follows:

```
python3 wrapper.py -r <reference_accession> -c <t/f>
```

Here is a brief explanation of the arguments:
```
-r/--reference: 
The accession number to be used as a reference genome/transcriptome. For this particular analysis, the accession NC_006273.2 is used, but in theory any accession can be used as an argument.

-c/--condensed:
This is where whether the pipeline will be run on the test dataset or pull the full dataset from NCBI bnefore running. Valid values are true/True/t/T or false/False/f/F. The false option allows the pipeline to run with the full dataset that it pulls from NCBI, while the true option allows the pipeline to run with the test dataset. 
```

So for this analysis, to run the pipeline on the full dataset, use the command:
```
python3 wrapper.py -r "NC_006273.2" -c F
```
This command iterates through the `input_srr` file and downloads/unpacks the appropriate sequence read archive runs from NCBI. 

To run on the sample dataset, use the command:
```
python3 wrapper.py -r "NC_006273.2" -c T
```
This command skips downloading and unpacking the appropriate sequence read archive runs from NCBI and instead conducts the rest of the analysis on the files contained in the `sample_data` directory. 

Upon running either of these commands for the first time, you will notice a directory called `processed_data` is created - please do not move/delete any files/directories in this directory while the pipeline is running - the pipeline relies on the structure of this folder in order to complete its analysis. Also, please keep in mind that every time that this pipeline is started, if the `processed_data` directory is present in the root directory, it and its contents will be deleted and replaced with empty versions to prevent data pollution in subsequent uses of the pipeline.  

The full output will be updated as the pipeline moves along in the file `processed_data/output_log/output.txt`, and minimal updates will be sent to the terminal. 

For this assignment, one of the requirements was the creation of a log file that outputted certain results from the analysis. The original log file can be found in `PipelineProject_Emil_Cacayan/PipelineProject.log`. To prevent any changes from going to this log file, a new directory with a new log file was made and the current implementation of this pipeline writes to this new log file (`PipelineProject/PipelineProject.log`). As a last note, for the BLAST alignment, creation of the database did not only include RefSeq accessions as that was not explicitly required. 

### RNA-Seq Transcriptomic Pipeline Workflow

This pipeline (source code contained in `wrapper.py`) takes RNA-seq viral transcriptomic information and performs the following analyses (many of the results will be printed to the `PipelineProject_Emil_Cacayan/PipelineProject.log` file - the directory and file will be created if not already present). For the sake of demonstration, this pipeline will be explained using the data listed in the [data](#data) section. 

1. **Quantification of transcripts per million (TPM) using [kallisto](https://pachterlab.github.io/kallisto/about) (Bray et al., 2016).**
    - First, a transcriptome index for HCMV (see [data](#data)) is built.
    - Alongside the creation of the transcriptome index, the coding sequence (CDS) features will be extracted using [Biopython tools](https://biopython.org) and a `cds.fasta` file will be created which contains the coding sequence (CDS) features from the genome using the RefSeq `protein_id` as the header .
    - The `cds.fasta` file will be used as the input for the kallisto index command. To the `PipelineProject.log` file will be written the number of coding sequences in the genome in the following fashion:
    ```
    The HCMV genome (NC_006273.2) has # CDS.
    ```
    where `#` is replaced with the number of CDS features. 
    - For each sample and condition, the minimum, median, mean, and maximum TPM will be calculated from the results in the `abundance.tsv` kallisto output file. The output will be written as a tab-delimited table in `PipelineProject.log` in the following fashion:
    ```
    sample  condition   min_tpm   med_tpm   mean_tpm   max_tpm
    ```
2. **Detection of differentially expressed genes using `R` package [sleuth](https://pachterlab.github.io/sleuth/about) (Pimentel et al., 2017).**
    - The output from kallisto will be used as input for sleuth to find differentially expressed genes between the two conditions. 
    - A transcript will be considered significant if the calculated false discovery rate (FDR) falls below 0.05 ($\text{FDR} < 0.05$).
    - The following details for each significant transcript will be written as a tab-delimited table to `PipelineProject.log` in the following fashion:
    ```
    target_id   test_stat   pval   qval
    ```
6. **Filtering of transcriptome reads using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (Langmead & Salzberg, 2012).**
    - To compare strain similarity between patient samples, the transcriptome reads will be assembled (we don't expect the entire genome to be produced - just enough useful information to be used in alignment tools).
    - To ensure that only viral reads are included (as sometimes viral transcriptomic data contains reads from the host), we will index the reads against the genome index. Using Bowtie2, we will create a genome index using the reference genome for the target virus (HCMV).
    - For the subsequent assembly, only reads that map to the index will be used.
    - The number of reads in each transcriptome before and after the Bowtie2 mapping will be written to `PipelineProject.log` in the following fashion:
    ```
    Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    Donor 1 (6dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    Donor 3 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    Donor 3 (6dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
    ```
6. **Assembly using [SPAdes](https://github.com/ablab/spades) (Bankevich et al., 2012).**
    - Using the output from Bowtie2, SPAdes will be utilized to generate two assemblies - one for each patient/donor with a k-mer size of 77.
    - The commands sent to the SPAdes assembler can be found in `PipelineProject.log`. 
7. **Determination of the strain of each assembly using [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (Camacho et al., 2009).** 
    - The pipeline will then retrieve the longest contig from each SPAdes assembly. 
    - This contig will be used as a blast+ input to query the nr nucleotide database limited to members of the *Betaherpesvirinae* subfamily. This is done by creating a local database of just sequences from the *Betaherpesvirinae* subfamily. 
    - The blast+ runs only keep the best alignment (HSP) for any single query-subject pair of sequences.
    - For the top hit runs for each assembly, the following are written to `PipelineProject.log`: subject accession, percent identity, alignment length, start of alignment in query, end of alignment in query, start of alignment in subject, end of alignment in subject, bit score, E-value, and subject title in a tab-delimited table in the following fashion:
    ```
    Donor1:
    sacc   pident    length   qstart   qend   sstart   send   bitscore   evalue   stitle

    Donor3:
    saccc   pident   length   qstart   qend   sstart   send   bitscore   evalue   stitle
    ```

## Disclaimer

This is a project for the computational biology class (COMP 483) taught by Dr. Heather Wheeler at Loyola University Chicago. I credit much of the structure of the pipeline to the objectives outlined in her class. I strongly discourage viewing the source code if you are currently in this class or similar classes and suggest you find and implement solutions to bioinformatic analysis using your own resources. I take no responsibility for any plagiarism that occurs from this repository as a result of my posting of this repository. 

## References

Arriaga-Canon, C., Contreras-Espinosa, L., Rebollar-Vega, R., Montiel-Manríquez, R., Cedro-Tanda, A., García-Gordillo, J. A., Álvarez-Gómez, R. M., Jiménez-Trejo, F., Castro-Hernández, C., & Herrera, L. A. (2022). Transcriptomics and RNA-Based Therapeutics as Potential Approaches to Manage SARS-CoV-2 Infection. International journal of molecular sciences, 23(19), 11058. https://doi.org/10.3390/ijms231911058

Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology : a journal of computational molecular cell biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. BMC bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

Cheng, S., Caviness, K., Buehler, J., Smithey, M., Nikolich-Žugich, J., & Goodrum, F. (2017). Transcriptome-wide characterization of human cytomegalovirus in natural infection and experimental latency. Proceedings of the National Academy of Sciences of the United States of America, 114(49), E10586–E10595. https://doi.org/10.1073/pnas.1710522114

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics (Oxford, England), 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163

Gatherer, D., Seirafian, S., Cunningham, C., Holton, M., Dargan, D. J., Baluchova, K., Hector, R. D., Galbraith, J., Herzyk, P., Wilkinson, G. W., & Davison, A. J. (2011). High-resolution human cytomegalovirus transcriptome. Proceedings of the National Academy of Sciences of the United States of America, 108(49), 19755–19760. https://doi.org/10.1073/pnas.1115861108

Ivanov, S. M., Tarasova, O. A., & Poroikov, V. V. (2023). Transcriptome-based analysis of human peripheral blood reveals regulators of immune response in different viral infections. Frontiers in immunology, 14, 1199482. https://doi.org/10.3389/fimmu.2023.1199482

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923

Mo, Z. Q., Li, Y. W., Wang, H. Q., Wang, J. L., Ni, L. Y., Yang, M., Lao, G. F., Luo, X. C., Li, A. X., & Dan, X. M. (2016). Comparative transcriptional profile of the fish parasite Cryptocaryon irritans. Parasites & vectors, 9(1), 630. https://doi.org/10.1186/s13071-016-1919-1

Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. Nature methods, 14(7), 687–690. https://doi.org/10.1038/nmeth.4324

Sudhagar, A., Kumar, G., & El-Matbouli, M. (2018). Transcriptome Analysis Based on RNA-Seq in Understanding Pathogenic Mechanisms of Diseases and the Immune System of Fish: A Comprehensive Review. International journal of molecular sciences, 19(1), 245. https://doi.org/10.3390/ijms19010245

Zhou, A., Dong, X., Liu, M., & Tang, B. (2021). Comprehensive Transcriptomic Analysis Identifies Novel Antiviral Factors Against Influenza A Virus Infection. Frontiers in immunology, 12, 632798. https://doi.org/10.3389/fimmu.2021.632798