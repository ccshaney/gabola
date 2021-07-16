# #RUN ALL in one command

Use -h show description: **`python /opt/10x_program/runAll.py -h`**

## Usage
```
usage: runStep1to3.py [-h] -f FASTQS -g GENOME --id PROJECTID -o OUTPUT_FOLDER
                      [-t THREADS] [-a {bwa_mem,kart}] [--check {False,True}]
                      [-q QUALITY] [-l LENGTH] [--trimr2 {False,True}]
                      [-d NONDUP] [-c CIGAR_MAP_QUALITY] [-m MAPQ] [-r RANGE]
                      [-gap WITHGAP] [--detail {False,True}]
                      [--min_rp THRESHOLD_OF_BX_RP_PER_END]
                      [--rpN ESTABLISH_COMBINATION_OF_READ_PAIR]
                      [--BXN ESTABLISH_NUMBER_OF_BX]

Combine THREE part of process: ...

```

# #1 Filter SAM file into selected CIGAR match% and MAPQ score:

**Process BWA / Kart output sam file into specific qualify reads**

Use -h show description: ** `python /opt/10x_program/step1_preprocessFastq.py -h`**

### Dependencies:

1. Python3+

## Usage
```
usage: step1_preprocessFastq.py [-h] -f FASTQS --id PROJECTID -o OUTPUT_FOLDER
                                [-c {False,True}] [-q QUALITY] [-l LENGTH]
                                [--trimr2 {False,True}] [-d NONDUP]
                                [-t THREADS]

Wraps 10X genomic's longranger basic to takes FASTQ files, which is created by longranger mkfastq, to performs basic barcode processing including error
correction, barcode white-listing, and attaching barcodes to reads. 
Then, use TrimGalore to apply adapter and quality trimming to FastQ files.
Then, optionally you can choose to trim 23mer of R2 prefix or not. 
Finally, process all fastq files to remove duplicate reads.

...

```

# #2 Mapping and Alignment

Use -h show description: ** `python /opt/10x_program/step2_alignment_andFilter.py -h`**

### Dependencies:

1. Python3+

2. BWA MEM / KART

### Process step:

1. Mapping pre-processed reads on to template scaffold by BWA MEM or KART.

2. Filter out un-qualified reads by CIGAR & MAPQ score.

## Usage:
```
usage: step2_alignment_andFilter.py [-h] -a {bwa_mem,kart} -g GENOME -f1
                                    FASTQ_R1 -f2 FASTQ_R2 -o OUTPUTPREFIX
                                    [-c CIGAR_MAP_QUALITY] [-m MAPQ]
                                    [-t THREADS]

To begin with align reads to genome scaffolds by BWA MEM or Kart. 
Then, filter proper pair by SAMBAMBA (v0.6.9). 
Finally, filter alignment by CIGAR match quality and MAPQ score. 
# proper pair define: not (unmapped or mate_is_unmapped or secondary_alignment or failed_quality_control or duplicate or supplementary or chimeric)

...

```

# #3 Process Genome Scaffold Head/Tail Information

Use -h show description: ** `python /opt/10x_program/step3_process_samfile.py -h`**

### Dependencies:

1. Python3+

### Process step:

1. Process sam file and genome fasta file to generate scaffold information

2. Mark scaffold information (Head/Tail) for each fastq reads

```
usage: step3_process_samfile.py [-h] -f FASTA [-r RANGE] -s SAM
                                [--max_number_of_scafendcnt SCAFENDCNT]
                                [--min_rp THRESHOLD_OF_BX_RP_PER_END]
                                [--rpN ESTABLISH_COMBINATION_OF_READ_PAIR]
                                [--BXN ESTABLISH_NUMBER_OF_BX] [-g WITHGAP]
                                [-t THREADS]

Process genome fasta file to get Head/Tail position for each scaffold. And use
this scaffold Head/Tail info with input sam file to process each reads'
position information.

...

```
