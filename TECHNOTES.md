## GABOLA Technical Notes
### § Preprocess Module:

#### *Step I. Produce non-duplicate split fastqs sorted by barcodes from raw linked reads*

**Input:**

- Path to the directory storing raw linked reads (see: https://github.com/10XGenomics/longranger for read file format)

**Output:**

- Non-duplicate, trimmed and filtered FASTQ files for Read 1 and Read 2 separately (NonDupR1.fq and NonDupR2.fq)
- Directories for every barcode containing their respective reads that are trimmed and filtered (OUTPUT_FOLDER/nonDupFq/split). 
> For example, Barcode name: AAAACCCCGTGTGTGT
> Format: OUTPUT_FOLDER/nonDupFq/split/AAAACCCCGTGTGTGT_1
```
usage: /opt/10x_program/step1_preprocessFastq.py [-h] -f FASTQS --id PROJECTID -o OUTPUT_FOLDER [-q QUALITY] 
                                                 [-l LENGTH] [--trimr2 {False,True}] [-d NONDUP] [-t THREADS]
positional arguments:
  -f FASTQS, --fastqs FASTQS
                        Input the folder path of FASTQ files which were
                        produced from longranger mkfastq. Files in the folder
                        must be named like: [Sample Name]_S1_L00[Lane
                        Number]_[Read Type]_001.fastq.gz. (see:
                        https://github.com/10XGenomics/longranger) If your
                        fastq files are in a different folder, please move
                        them to the same folder. (input file must be in fastq
                        format, such as .fq or .fastq)
  --id PROJECTID        Input project name. It will be your output prefix.
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder path. Your output files will be in this
                        folder. ps. longranger will output your folder at pwd,
                        so you must cd to your output folder.
                        
optional arguments:
  -q QUALITY, --quality QUALITY
                        Input the minimum quality score for quality trimming
                        of FASTQ files. Same as the option for trim_galore -q
                        [default:20]
  -l LENGTH, --length LENGTH
                        Input the minimum length for quality trimming of FASTQ
                        files. Same as the option for trim_galore --length
                        [default:50]
  --trimr2 {False,True}
                        (Optional) R2 trimming. For better mapping quality,
                        choose to trim barcode (16 mer) and adapter (6+1 mer),
                        overall 23 mer, from R2 prefix. If you choose not to
                        trim R2, set this parameter to False or 0.
                        [default:True]
  -d NONDUP, --deduplicate NONDUP
                        Get non-duplicate read pairs for each barcode and keep
                        BX size larger than 3 reads. If you do not want to de-
                        duplicate, please set "-d 0". [default:3]
  -t THREADS, --threads THREADS
                        Number of threads. Input 0 will use all
                        threads.(default: 1)
  -h, --help            show this help message and exit

```
#### *Step II. Filter alignment by CIGAR match quality and MAPQ score*

**Input:**

- Draft assembly

**Output:**

- A SAM file of filtered read-to-scaffold alignment with CIGAR match quality >= 70 and Mapping Quality = 60 
> File name ends with *C70M60.sam*
```
usage: /opt/10x_program/step2_alignment_andFilter.py [-h] -a {bwa_mem,kart} -g GENOME -f1 FASTQ_R1 -f2 FASTQ_R2 
                                                     -o OUTPUTPREFIX [-c CIGAR_MAP_QUALITY] [-m MAPQ] [-t THREADS]
positional arguments:
  -a {bwa_mem,kart}, --aligner {bwa_mem,kart}
                        Align reads to scaffolds. Choose one of the following
                        aligners: BWA MEM or Kart.
  -g GENOME, --genome GENOME
                        Path of genome fasta.
  -f1 FASTQ_R1, --fastqR1 FASTQ_R1
                        Path of qualified and reformated Fastq R1.
  -f2 FASTQ_R2, --fastqR2 FASTQ_R2
                        Path of qualified and reformated Fastq R2.
  -o OUTPUTPREFIX, --outputPrefix OUTPUTPREFIX
                        Prefix of output samfile. Do not use "." in this
                        prefix. ie, -o PROJ_NAME will get PROJ_NAME.sam &
                        PROJ_NAME_ProperPair.sam & PROJ_NAME_C70M60.sam
optional arguments:
  -h, --help            show this help message and exit
  -c CIGAR_MAP_QUALITY, --CIGAR CIGAR_MAP_QUALITY
                        CIGAR match quality: count([M=X])/len(seq)*100. From 1
                        to 100. [default: 70]
  -m MAPQ, --MAPQ MAPQ  MAPQ score. From 0 to 60. [default: 60]
  -t THREADS, --threads THREADS
                        Number of threads. Input 0 will use all
                        threads.[default: 1]
```
#### *Step III. Generate barcode information on each scaffold*

**Input:**
- Read-to-assembly filtered SAM file from Step II

**Output:**
- A TSV file for possible scaffold pairs and shared barcode counts 
> File name ends with *C70M60_ScafA_ScafB_BXCnt.tsv*
- A TSV file for barcodes on every scaffold end
> File name ends with *C70M60_ScafHeadTail_BX_pairSum.tsv*

```
usage: /opt/10x_program/step3_process_samfile.py [-h] -f FASTA [-r RANGE] -s SAM
                                                 [--max_number_of_scafendcnt SCAFENDCNT]
                                                 [--min_rp THRESHOLD_OF_BX_RP_PER_END] 
                                                 [--rpN ESTABLISH_COMBINATION_OF_READ_PAIR]
                                                 [--BXN ESTABLISH_NUMBER_OF_BX] [-g WITHGAP] [-t THREADS]

positional arguments:
  -f FASTA, --fasta FASTA
                        Reference genome. (fasta format, named as .fa or .fasta)
  -s SAM, --sam SAM     Path of the SAM file you want to process.

optional arguments:
  -h, --help            show this help message and exit
  -r RANGE, --range RANGE
                        Range length of each scaffold Head/Tail.
                        [default:20000]
  --max_number_of_scafendcnt SCAFENDCNT
                        Max number of Scaffold ends sharing the same BX.
                        [default:0]
  --min_rp THRESHOLD_OF_BX_RP_PER_END
                        Threshold of BX read pair/per end. The certified
                        barcode control value. If the value is set to 5 means
                        barcode mapping to genome must be supported by 6 read
                        pair. [default:0]
  --rpN ESTABLISH_COMBINATION_OF_READ_PAIR
                        Establish Number of barcode read pair at each scaffold
                        end. If value is set to 10 means this program will
                        generate combination from (1,1) to (10*,10*).
                        [default:10]
  --BXN ESTABLISH_NUMBER_OF_BX
                        Determines the number of barcode supporting pairs of
                        scaffold end. If the value is set to 20 means this
                        program will generate BX1 to BX20*. [default:20]
  -g WITHGAP, --withgap WITHGAP
                        Determines whether the Head/Tail range is with or
                        without gap. True means that the Head/Tail range will
                        be counted with gap. If you input [-r 20000], then the
                        Head/Tail range will be calculated from both ends
                        (Head/Tail) to 20000 base, and there may be some gaps
                        in said range(N base). False means that the Head/Tail
                        range will count without gap (N base). (default: True)
  -t THREADS, --threads THREADS
                        Number of threads you want to use. Input 0 will use
                        all threads.(default: 1)

```
#### *Run Step I to III altogether*
This program is the combination of the three parts of the preprocess module:

```
usage: /opt/10x_program/runStep1to3.py [-h] -f FASTQS -g GENOME --id PROJECTID -o OUTPUT_FOLDER
                      [-t THREADS] [-a {bwa_mem,kart}] [-q QUALITY]
                      [-l LENGTH] [--trimr2 {False,True}] [-d NONDUP]
                      [-c CIGAR_MAP_QUALITY] [-m MAPQ] [-r RANGE]
                      [-gap WITHGAP] [--detail {False,True}]
                      [--min_rp THRESHOLD_OF_BX_RP_PER_END]
                      [--rpN ESTABLISH_COMBINATION_OF_READ_PAIR]
                      [--BXN ESTABLISH_NUMBER_OF_BX]
positional arguments:
  -f FASTQS, --fastqs FASTQS
                        Input the folder path of FASTQ files which were
                        produced from longranger mkfastq. Files in the folder
                        must be named like: [Sample Name]_S1_L00[Lane
                        Number]_[Read Type]_001.fastq.gz. (see:
                        https://github.com/10XGenomics/longranger). If your
                        fastq files are in a different folder, please move to
                        the same folder. (input file must be in fastq format,
                        such as .fq or .fastq)
  -g GENOME, --genome GENOME
                        Path of genome fasta.
  --id PROJECTID        Input project name. It will be your output prefix.
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Output folder path. Your output files will be in this
                        folder.
                        
optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads. Input 0 will use all
                        threads.[default: 1]
  -a {bwa_mem,kart}, --aligner {bwa_mem,kart}
                        Choose one of the following aligners: BWA MEM or Kart,
                        to align reads. BWA MEM: https://github.com/lh3/bwa,
                        KART: https://github.com/hsinnan75/Kart.
                        [default:bwa_mem]
  -q QUALITY, --quality QUALITY
                        Input the minimum quality score for quality trimming
                        of FASTQ files. Same as the option for trim_galore -q
                        [default:20]
  -l LENGTH, --length LENGTH
                        Input the minimum length for quality trimming of FASTQ
                        files. Same as the option for trim_galore --length
                        [default:50]
  --trimr2 {False,True}
                        (Optional) R2 trimming. For better mapping quality,
                        choose to trim barcode (16 mer) and adapter (6+1 mer),
                        overall 23 mer, from R2 prefix. If you choose not to
                        trim R2, set this parameter to False or 0.
                        [default:True]
  -d NONDUP, --deduplicate NONDUP
                        (Optional) Get non-duplicate read pairs for each
                        barcode and keep BX size larger than 3 reads. If you
                        do not want to de-duplicate, please set "-d 0".
                        [default:3]
  -c CIGAR_MAP_QUALITY, --CIGAR CIGAR_MAP_QUALITY
                        CIGAR match quality: count([M=X])/len(seq)*100. From 1
                        to 100. [default: 70]
  -m MAPQ, --MAPQ MAPQ  MAPQ score. From 0 to 60. [default: 60]
  -r RANGE, --range RANGE
                        The range length of each scaffold Head/Tail.
                        [default:20000]
  -gap WITHGAP, --withgap WITHGAP
                        Determines whether the Head/Tail range is with or
                        without gap. True means the Head/Tail range will be
                        counted with gap. If you input [-r 20000], then the
                        Head/Tail range will be calculated from both ends
                        (Head/Tail) to 20000 base, and there may be some gaps
                        in that range (N base). False means the Head/Tail
                        range will be counted without gap (N base). [default:
                        True]
  --detail {False,True}
                        Print out the details of process log. [default: False]
  --min_rp THRESHOLD_OF_BX_RP_PER_END
                        Threshold of BX read pair/per end. The certified
                        barcode control value. If the value is set as 5, this
                        means the number of barcodes mapping to genome must be
                        supported by 6 read pairs. [default:0]
  --rpN ESTABLISH_COMBINATION_OF_READ_PAIR
                        Establish Number of barcode read pairs at each end of
                        the scaffolds. If the value is 10, this means this
                        program will generate a combination ranging from (1,1)
                        to (10*,10*). [default:10]
  --BXN ESTABLISH_NUMBER_OF_BX
                        Determine the number of barcodes supporting scaffold
                        ends. If the value is set as 20, this means this
                        program will generate BX1 to BX20*. [default:20]

```
### § Main Module 1:
