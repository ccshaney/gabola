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
