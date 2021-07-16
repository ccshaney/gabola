import argparse
import os
from src.gapfill_lib.funlib import proper_threads, run_command_with_popen_communicate
# from src.gapfill_lib.funlib import run_command_with_popen_communicate_onlyReturn

# ###############
# argparse part:
# ###############

version_number = '0.0.1'
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=
'''This program is the combination of the three parts of the preprocess module:

   During the first part, we perform FASTQ preprocessing:
   It takes FASTQ files created by longranger mkfastq
   and perform basic barcode processing including error correction,barcode white-listing,
   and attaching barcodes to reads. 
   Then, we use TrimGalore to apply adapter and quality trimming to FastQ files.
   After that, you can choose to trim 23mer of R2 prefix or not. 
   Finally, we process all fastq files to remove duplicate reads. 
   (See help: python /opt/10x_program/step1_preprocessFastq.py -h for more detail)

   During the second part, we map reads to genome and filter them:
   We begin with aligning reads to genome scaffolds by BWA MEM or Kart.
   Then, we filter proper pairs by SAMBAMBA (v0.6.9). 
       #definition of proper pair:
       not (unmapped or mate_is_unmapped or secondary_alignment or 
       failed_quality_control or duplicate or supplementary or chimeric)
   Finally, we filter the alignment by match quality calculated from CIGAR and MAPQ score.
   (See help: python /opt/10x_program/step2_alignment_andFilter.py -h for more detail)

   During the third part, we process SAM file to get each read's position information.
   (See help: python /opt/10x_program/step3_process_samfile.py -h for more detail) ''',
                                 epilog='''If you run into any questions or have any suggestions,please contact Fu Po-Ying <spashleyfu@gmail.com> 
process_sam_to_HeadTailInfo.py version: {}'''.format(version_number))

### required:
parser.add_argument('-f', '--fastqs', dest="fastqs", required=True, type=str,
                    help='''Input the folder path of FASTQ files which were produced from longranger mkfastq.
                    Files in the folder must be named like:
                    [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz.
                    (see: https://github.com/10XGenomics/longranger). 
                    If your fastq files are in a different folder, please move to the same folder.
                    (input file must be in fastq format, such as .fq or .fastq)''')
                    
parser.add_argument('-g', '--genome', dest="genome", required=True, type=str,
                    help='Path of genome fasta.')
parser.add_argument('--id', dest="projectID", type=str, required=True,
                    help='''Input project name. It will be your output prefix.''')
parser.add_argument('-o', '--output_folder', dest="output_folder", type=str, required=True,
                    help='''Output folder path. Your output files will be in this folder.''')

# optional:
parser.add_argument('-t', '--threads', dest="threads", type=int, default=1,
                    help='Number of threads. Input 0 will use all threads.[default: 1]')
parser.add_argument('-a', '--aligner', dest="aligner", default='bwa_mem', type=str,
                    choices=['bwa_mem', 'kart'],
                    help='''Choose one of the following aligners: BWA MEM or Kart, to align reads.
                    BWA MEM: https://github.com/lh3/bwa, 
                    KART: https://github.com/hsinnan75/Kart.
                    [default:bwa_mem]''')

parser.add_argument('-q', '--quality', dest="quality", type=int, default=20,
                    help='''Input the minimum quality score for quality trimming of FASTQ files.
                    Same as the option for trim_galore -q [default:20]''')
parser.add_argument('-l', '--length', dest="length", type=int, default=50,
                    help='''Input the minimum length for quality trimming of FASTQ files. 
                    Same as the option for trim_galore --length [default:50]''')
parser.add_argument('--trimr2', dest="trim_r2",
                    type=lambda x: False if str(x).lower() == 'false' else True,
                    default=True, choices=[False, True],
                    help='''(Optional) R2 trimming. For better mapping quality, choose to trim barcode (16 mer) and adapter (6+1 mer), 
                    overall 23 mer, from R2 prefix. If you choose not to trim R2, set this parameter to False or 0. [default:True]''')
parser.add_argument('-d', '--deduplicate', dest="nonDup", type=int, default=3,
                    help='''(Optional) Get non-duplicate read pairs for each barcode and keep BX size larger than 3 reads.
                    If you do not want to de-duplicate, please set "-d 0". [default:3]''')
parser.add_argument('-c', '--CIGAR', dest="cigar_map_quality", default=70, type=int,
                    help='''CIGAR match quality: count([M=X])/len(seq)*100. 
                    From 1 to 100. [default: 70]''')
parser.add_argument('-m', '--MAPQ', dest="mapq", type=int, default=60,
                    help='''MAPQ score. From 0 to 60. [default: 60]''')
parser.add_argument('-r', '--range', dest="range", type=int, default=20000,
                    help='The range length of each scaffold Head/Tail. [default:20000]')
parser.add_argument('-gap', '--withgap', dest="withgap", type=bool, default=True,
                    help='''Determines whether the Head/Tail range is with or without gap.
                    True means the Head/Tail range will be counted with gap. 
                    If you input [-r 20000], then the Head/Tail range will be calculated from both ends (Head/Tail) to 20000 
                    base, and there may be some gaps in that range (N base).
                    False means the Head/Tail range will be counted without gap (N base). [default: True]''')
parser.add_argument('--detail', dest="detail",
                    type=lambda x: False if str(x).lower()=='false' else True,
                    default=False, choices=[False, True],
                    help='''Print out the details of process log. [default: False]''')
# For calculate combination:
parser.add_argument('--min_rp', dest="threshold_of_BX_rp_per_end", type=int, default=0,
                    help='''Threshold of BX read pair/per end. The certified barcode control value. 
                    If the value is set as 5, this means the number of barcodes mapping to genome must be supported by 6 read pairs. [default:0]''')
parser.add_argument('--rpN', dest="establish_combination_of_read_pair", type=int, default=10,
                    help='''Establish Number of barcode read pairs at each end of the scaffolds.
                    If the value is 10, this means this program will generate a combination ranging from (1,1) to (10*,10*). [default:10]''')
parser.add_argument('--BXN', dest="establish_number_of_bx", type=int, default=20,
                    help='''Determine the number of barcodes supporting scaffold ends.
                    If the value is set as 20, this means this program will generate BX1 to BX20*. [default:20]''')

options = parser.parse_args()
# --detail:
detail = options.detail
# --output_folder:
output_folder = options.output_folder
output_folder = os.path.abspath(output_folder) + '/'

### 1. step1_preprocessFastq.py:
# --fastqs:
input_fastqs_folder = options.fastqs # ie. /data/PROJ/
input_fastqs_folder = os.path.abspath(input_fastqs_folder)
# --id:
project_id = options.projectID # ie. PROJ_NAME
# --check:
input_file_check = options.check_input_file
# --output_galore: Same as output_folder.
# --quality:
fastq_quality = int(options.quality)
# --length:
fastq_length = int(options.length)
# --trimr2:
trimR2 = options.trim_r2
# --deduplicate:
nonDup_int = int(options.nonDup)
# --threads:
threads = int(options.threads)
threads = proper_threads(threads)

### 2. step2_alignment_andFilter.py:
# --aligner:
aligner = options.aligner
# --genome:
input_genomeFasta = options.genome # ie. /data/PROJ/genome.fasta
input_genomeFasta = os.path.abspath(input_genomeFasta)
# --fastqR1:
input_R1Fastq = output_folder + 'nonDupFq/NonDupR1.fq'
# --fastqR2:
input_R2Fastq = output_folder + 'nonDupFq/NonDupR2.fq'
# --outputPrefix: # PROJ_NAME
output_prefixofSamFile = output_folder + project_id
# --CIGAR:
cigar_in = int(options.cigar_map_quality)
# --MAPQ:
mapQ_in = int(options.mapq)
# --threads: same as threads

### 3. step3_process_samfile.py:
# --fasta: same as input_genomeFasta
# --range:
head_tail_range = options.range
# --sam:
input_samfile_path = output_prefixofSamFile + '_{}_C{}M{}.sam'.format(aligner, cigar_in, mapQ_in)
# --min_rp:
dropOut = options.threshold_of_BX_rp_per_end
# --rpN:
rpN = options.establish_combination_of_read_pair
# --BXN:
BXNumb = options.establish_number_of_bx
# --withgap:
with_or_without = options.withgap
# --threads: same as threads

################
# Process Step:
################
### CHECK INPUT FILE:
# check FASTQ folder:
assert os.path.exists(input_fastqs_folder), "Your Input folder path does not exists! Please Check it."
assert os.path.isdir(input_fastqs_folder), "Your Input path is not a folder! Please Check it."

### STEP1. python /opt/10x_program/step1_preprocessFastq.py -h
print('\n[PROCESS STEP 1]...')
Step1_commandStr = 'python /opt/10x_program/step1_preprocessFastq.py ' \
                   '-f {} --id {} -o {} -q {} -l {} --trimr2 {} -d {} -t {}'.format(input_fastqs_folder,
                                                                                    project_id,
                                                                                    input_file_check,
                                                                                    output_folder,
                                                                                    fastq_quality,
                                                                                    fastq_length,
                                                                                    trimR2,
                                                                                    nonDup_int,
                                                                                    threads)

split_fastq_folder = output_folder + 'nonDupFq/split/'
split_fastq_folder = os.path.abspath(split_fastq_folder)

if os.path.exists(split_fastq_folder):
    print('\nThe directory for Non-duplicate Split Fastqs already exists! Skipping PROCESS STEP 1')
elif detail:
    print(run_command_with_popen_communicate(Step1_commandStr))
else:
    run_command_with_popen_communicate(Step1_commandStr)


### STEP2. python /opt/10x_program/step2_alignment_andFilter.py -h
# Check Fastq R1 / R2:
assert os.path.exists(input_R1Fastq), "Your Input R1 FASTQ file does not exists! Please Check it."
assert os.path.exists(input_R2Fastq), "Your Input R2 FASTQ file does not exists! Please Check it."

print('\n[PROCESS STEP 2]...')
Step2_commandStr = 'python /opt/10x_program/step2_alignment_andFilter.py ' \
                   '-a {} -g {} -f1 {} -f2 {} -o {} -c {} -m {} -t {}'.format(aligner,
                                                                              input_genomeFasta,
                                                                              input_R1Fastq,
                                                                              input_R2Fastq,
                                                                              output_prefixofSamFile,
                                                                              cigar_in,
                                                                              mapQ_in,
                                                                              threads)
if detail:
    print(run_command_with_popen_communicate(Step2_commandStr))
else:
    run_command_with_popen_communicate(Step2_commandStr)

### STEP3. python /opt/10x_program/step3_process_samfile.py -h
# Check SAM file:
assert os.path.exists(input_samfile_path), "Your Input sam file (CIGAR {} MAPQ {}) does not exists! Please Check it.".format(cigar_in, mapQ_in)

print('\n[PROCESS STEP 3]...')
Step3_commandStr = 'python /opt/10x_program/step3_process_samfile.py ' \
                   '-f {} -r {} -s {} -g {} -t {} --min_rp {} --rpN {} --BXN {}'.format(input_genomeFasta,
                                                                                        head_tail_range,
                                                                                        input_samfile_path,
                                                                                        with_or_without,
                                                                                        threads,
                                                                                        dropOut, rpN, BXNumb)
if detail:
    print(run_command_with_popen_communicate(Step3_commandStr))
else:
    run_command_with_popen_communicate(Step3_commandStr)

print("[PROCESS ALL DONE]... Finished!")
