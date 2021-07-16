import argparse
import subprocess
import os
import sys

# #############
# ERROR class:
# #############

class InputPathError(Exception):
    pass

class FileFormatError(Exception):
    pass

class ValuesError(Exception):
    pass

# ##########
# Function:
# ##########

def find_proper_fastq_byname(file_list):
    import re
    I1 = 0
    R1 = 0
    R2 = 0
    # pattern: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
    I1_prog = re.compile(r"^[\w\-_\.]+_S\d+_L\d+_I1_001\.fastq(\.gz)?$")
    R1_prog = re.compile(r"^[\w\-_\.]+_S\d+_L\d+_R1_001\.fastq(\.gz)?$")
    R2_prog = re.compile(r"^[\w\-_\.]+_S\d+_L\d+_R2_001\.fastq(\.gz)?$")
    file_dict = {'index': [], 'r1': [], 'r2': []}
    for file in file_list:
        if I1_prog.match(file):
            I1 += 1
            file_dict['index'].append(file)
        elif R1_prog.match(file):
            R1 += 1
            file_dict['r1'].append(file)
        elif R2_prog.match(file):
            R2 += 1
            file_dict['r2'].append(file)
    return [I1, R1, R2, file_dict]

def check_input_fastq_name(file_list):
    # R1/R2 file: must has same number of R1 and R2 file.
    # Index file: must has same number of index file with R1/R2 file or no index file!
    I1, R1, R2, file_dict = find_proper_fastq_byname(file_list)
    #print("I1 = {}, R1 = {}, R2 = {}, file_dict = {}".format(I1, R1, R2, file_dict))
    if len(file_list) == 0:
        raise InputPathError("Your Input folder is empty! Please check your input folder.")
    if R1 > 0 and R1 == R2:
        if I1 == R1:
            return 'hasIndex'
        elif I1 == 0:
            return 'noIndex' # this condition contain: no index file and some of index file.
        else:
            raise InputPathError("""Number of index file is not match number of R1/R2 file.
            For more information please see: 
            https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/fastq-input""")
    elif R1 == 0 or R2 == 0:
        raise InputPathError("Cannot find R1 or R2 file in input folder.")
    else:
        raise InputPathError("""Cannot find proper files in folder. Please Check your input folder.
        The files must named like: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz.
        [Read Type] must have R1 and R2, and I1 can generate by this program. 
        For more information please see: 
        https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/fastq-input""")
    
def create_index(fastqs_folder, prog_path):
    # it must: (R1==R2 and R1,R2>0) && (I1==0)
    import subprocess
    import os
    import re
    files_list = os.listdir(fastqs_folder)
    I1, R1, R2, file_dict = find_proper_fastq_byname(files_list)
    prog = re.compile(r"^[\w\-_\.]+_S\d+_L\d+")
    # Use R1 file to generate index file:
    print("\nInput folder dosen't have INDEX file.\n[CREATE INDEX]...")
    for r1file in file_dict['r1']:
        m = prog.search(r1file)
        projectName = r1file[m.start():m.end()]
        print("Start to create {}'s index file...".format(projectName))
        indexName = fastqs_folder + projectName + '_I1_001.fastq.gz'
        r1fileFullPath = fastqs_folder + r1file
        # Check is a gz file or not: call different script.
        gz_prog = re.compile(r"\.fastq.gz$")
        if gz_prog.search(r1file):
            # Use gz's shell script: create_10xIndex.sh
            print("sh {}create_10xIndex.sh {} {}".format(prog_path, r1fileFullPath, indexName))
            createIndex_gz = subprocess.Popen(
                "sh {}create_10xIndex.sh {} {}".format(prog_path, r1fileFullPath, indexName),
                shell=True, stdout=subprocess.PIPE)
            createIndex_gz.wait()
            print("{}".format(createIndex_gz.stdout.read().decode()))
        else:
            # Use other script: create_10xIndex_noZip.sh
            print("sh {}create_10xIndex_noZip.sh {} {}".format(prog_path, r1fileFullPath, indexName))
            createIndex = subprocess.Popen(
                "sh {}create_10xIndex_noZip.sh {} {}".format(prog_path, r1fileFullPath, indexName),
                shell=True, stdout=subprocess.PIPE)
            createIndex.wait()
            print("{}".format(createIndex.stdout.read().decode()))

def proper_threads(input_num):
    import os
    # Get total CPU number:
    cpu = os.cpu_count()
    if input_num < cpu and input_num != 0:
        cpu = input_num
    print("\n[RUN THREADS]... {} cores".format(cpu))
    return cpu

def run_command_with_popen(commandStr):
    import subprocess
    print(commandStr)
    command = subprocess.Popen(commandStr, shell=True, stdout=subprocess.PIPE)
    command.wait()
    print("{}".format(command.stdout.read().decode()))

def run_popen_with_return(commandStr):
    import subprocess
    command = subprocess.Popen(commandStr, shell=True, stdout=subprocess.PIPE)
    command.wait()
    return command.stdout.read().decode()

# ###############
# argparse part:
# ###############

version_number = '0.0.1'
parser = argparse.ArgumentParser(description='''Wraps 10X genomic's longranger basic to takes FASTQ files, 
which is created by longranger mkfastq, to performs basic barcode processing including error correction, 
barcode white-listing, and attaching barcodes to reads. Then, use TrimGalore to apply adapter and quality 
trimming to FastQ files. Then, optionally you can choose to trim 23mer of R2 prefix or not. 
Finally, process all fastq files to remove duplicate reads.''',
                                 epilog='''Any questions or comments, please contact Fu Po-Ying <spashleyfu@gmail.com> 
process_sam_to_HeadTailInfo.py version: {}'''.format(version_number))

parser.add_argument('-f', '--fastqs', dest="fastqs", required=True, type=str,
                    help='''Input the folder path of FASTQ files which produced from longranger mkfastq. 
                    Same as the option of longranger basic --fastqs, and files in the folder must named like:
                    [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
                    (see: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/fastq-input)
                    If your fastq files are in different folder, please move to same folder.
                    (input file must be fastq format, named as .fq or .fastq)''')
parser.add_argument('-id', '--id', dest="projectID", type=str, required=True,
                    help='''Input project name. Same as the option of longranger basic --id.''')
parser.add_argument('-c', '--check', dest="check_input_file", 
                    type=lambda x: False if (str(x).lower() in ['false', '0']) else True, 
                    default=True,
                    help='''Check fastqs input file naming format first, then send input file folder (-f) to longranger.
                    If your file naming format is not like: [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz, 
                    but yous sure the longranger basic can process, you can set this parameter to False or 0, 
                    then the program will do the bypass function to run longranger basic directly. 
                    [default:True]''')
parser.add_argument('-o', '--output', dest="outputFolder", type=str, default='notSet',
                    help='''Output folder path. Same as the option of trim_galore -o [default: same as -f input path]''')
parser.add_argument('-q', '--quality', dest="quality", type=int, default=20,
                    help='''Input the minimum quality score for quality trimming to FASTQ files.
                    Same as the option of trim_galore -q [default:20]''')
parser.add_argument('-l', '--length', dest="length", type=int, default=50,
                    help='''Input the minimum length for quality trimming to FASTQ files. 
                    Same as the option of trim_galore --length [default:50]''')
parser.add_argument('-r2', '--trimr2', dest="trim_r2", 
                    type=lambda x: False if (str(x).lower() in ['false', '0']) else True, 
                    default=True,
                    help='''(Optional) R2 trimming. For better mapping quality, choose to trim barcode (16 mer) and adapter (6+1 mer), 
                    overall 23 mer, from R2 prefix. If choose not to trim R2, set this parameter to False or 0. [default:True]''')
parser.add_argument('-d', '--deduplicate', dest="nonDup", type=int, default=3,
                    help='''Get non-duplicate read pairs for each barcode and keep BX size large than 3 reads.
                    If you do not want to de-duplicate, please set "-d 0". [default:3]''')
parser.add_argument('-t', '--threads', dest="threads", type=int, default=1,
                    help='Number of threads. Input 0 will use all threads.(default: 1)')

options = parser.parse_args()
input_fastqs_folder = options.fastqs
output_folder = options.outputFolder
project_id = options.projectID
input_file_check = options.check_input_file
trimR2 = options.trim_r2
fastq_quality = int(options.quality)
fastq_length = int(options.length)
nonDup_int = int(options.nonDup)
threads = int(options.threads)

# ######
# Main:
# ######

# check input threads number:
cores = proper_threads(threads)

# ie. prog_path = '/opt/10x_program/src/'
prog_path = os.path.dirname(os.path.abspath(sys.argv[0])) + '/src/'
print("prog_path = {}".format(prog_path))

# check -o: If not input any, set the path same as -f.
if output_folder == 'notSet':
    output_folder = input_fastqs_folder

# get absolute path:
input_fastqs_folder = os.path.abspath(input_fastqs_folder) + '/'
output_folder = os.path.abspath(output_folder) + '/'

# get pwd:
pwd = os.getcwd()

if os.path.exists(input_fastqs_folder):
    if os.path.isdir(input_fastqs_folder):
        files_list = os.listdir(input_fastqs_folder)
        # Check input file:
        if input_file_check:
            # 1. Check input folder has index files or not?
            #     If not has index file, create index files!
            #     If has index, do nothing.
            print("\n[CHECK INPUT FOLDER]...")
            if check_input_fastq_name(files_list) == 'noIndex':
                create_index(input_fastqs_folder, prog_path)
            print("DONE.")
        else:
            print("\n[CHOOSE NOT TO CHECK INPUT FILE]...")
        # 2. RUN longranger basic:
        print("\n[CALL LONGRANGER BASIC] it may take hours...".format())
        runLongranger = "longranger basic --id={} --fastqs={} --localcores={}".format(project_id,
                                                                                      input_fastqs_folder, 
                                                                                      cores)
        run_command_with_popen(runLongranger)
        longranger_outFQ_path = pwd + '/' + project_id + '/outs/barcoded.fastq.gz' # Output will generate at pwd
        print("longranger OUTPUT file at {}\nDONE.".format(longranger_outFQ_path))

        # 3. Separate FQ to R1 & R2: split_longrangerFQ_toR1R2.sh
        raw_R1_path = input_fastqs_folder + project_id + '_barcoded_R1.fastq'
        raw_R2_path = input_fastqs_folder + project_id + '_barcoded_R2.fastq'
        print("\n[SPLIT LONGRANGER OUTPUT INTO R1/R2]...".format())
        runSplitFQ = "sh {}split_longrangerFQ_toR1R2.sh {} {} {}".format(prog_path,
                                                                  longranger_outFQ_path, 
                                                                  raw_R1_path, 
                                                                  raw_R2_path)
        outResult = run_popen_with_return(runSplitFQ)
        # process output info:
        emp1, r1_out, r2_out, emp2 = outResult.split("\n")
        r1_readCnt, r1_readCntWithBX = r1_out.strip().split(", ")
        r2_readCnt, r2_readCntWithBX = r2_out.strip().split(", ")
        # Get each number:
        r1_readCnt = int(r1_readCnt.replace('TOTAL READ: ',''))
        r1_readCntWithBX = int(r1_readCntWithBX.replace('TOTAL READ WITH BX: ',''))
        r2_readCnt = int(r2_readCnt.replace('TOTAL READ: ',''))
        r2_readCntWithBX = int(r2_readCntWithBX.replace('TOTAL READ WITH BX: ',''))

        if r1_readCnt == 0 or r2_readCnt == 0:
            raise ValuesError('''[ERROR] Something went wrong! 
            The fastq file may be EMPTY, or /opt/10x_program/src/split_longrangerFQ_toR1R2_v1.sh dosen't count correctly.''')
        elif r1_readCntWithBX == 0 or r2_readCntWithBX == 0:
            raise ValuesError('''[ERROR] Your FASTQ file dosen't have any barcode...
            It may cause by file format uncorrect.
            This Program search for barcode by /BX:Z:[ATGC]{16}-1/''')
        else:
            BXPercentageofR1 = str(round((r1_readCntWithBX/r1_readCnt)*100, 2)) + '%'
            BXPercentageofR2 = str(round((r2_readCntWithBX/r2_readCnt)*100, 2)) + '%'

        if r1_readCnt != r2_readCnt or r1_readCntWithBX != r2_readCntWithBX:
            raise FileFormatError("""[ERROR] Something WRONG with FASTQ files, the read is not paired!
            Please check your INPUT file, which should be the out file from LONGRANGER!""")
        else:
            print("\nR1 - Number of read = {} (100%), Number of read with BX = {} ({})".format(r1_readCnt, r1_readCntWithBX, BXPercentageofR1))
            print("R2 - Number of read = {} (100%), Number of read with BX = {} ({})".format(r2_readCnt, r2_readCntWithBX, BXPercentageofR2))
        # 4. trim_galore:
        print("\n[RUN TRIM GALORE]...")
        runTrimGalore = "trim_galore -q {} --phred33 -a AGATCGGAAGAGC --length {} --paired {} {} -o {}".format(fastq_quality,
                                                                                                               fastq_length, 
                                                                                                               raw_R1_path, 
                                                                                                               raw_R2_path, 
                                                                                                               input_fastqs_folder)
        run_command_with_popen(runTrimGalore)
        # Print out output file info:
        val_trimgalore_R1 = input_fastqs_folder + project_id + '_barcoded_R1_val_1.fq'
        val_trimgalore_R2 = input_fastqs_folder + project_id + '_barcoded_R2_val_2.fq'
        checkTrimGaloreOut = "ls -lh {} {}".format(val_trimgalore_R1, val_trimgalore_R2)
        run_command_with_popen(checkTrimGaloreOut)
        
        # 5. Trim R2 23mer: (Optional)
        trimR2_out = input_fastqs_folder + project_id + '_barcoded_R2_val_trim23mer.fq'
        if trimR2:
            print("\n[START TO TRIM R2 23mer]...")
            runTrimR2 = "sh {}trimR2_23mer.sh {} {}".format(prog_path, val_trimgalore_R2, trimR2_out)
            run_command_with_popen(runTrimR2)
        
        # 6. Filter out the read without BX and reFormat readID:
        print("\n[START TO FILTER OUT READ WITHOUT BX AND REFORMAT]...")
        # R1:
        keepBX_R1 = input_fastqs_folder + project_id + '_barcoded_R1_val_reformat.fq'
        runKeepBXR1 = "sh {}keepFQwithBXandReformat.sh {} {}".format(prog_path, val_trimgalore_R1, keepBX_R1)
        run_command_with_popen(runKeepBXR1)
        # R2:
        if trimR2:
            keepBX_R2 = input_fastqs_folder + project_id + '_barcoded_R2_val_trim23mer_reformat.fq'
            runKeepBXR2Trim = "sh {}keepFQwithBXandReformat.sh {} {}".format(prog_path, trimR2_out, keepBX_R2)
            run_command_with_popen(runKeepBXR2Trim)
        else:
            keepBX_R2 = input_fastqs_folder + project_id + '_barcoded_R2_val_reformat.fq'
            runKeepBXR2 = "sh {}keepFQwithBXandReformat.sh {} {}".format(prog_path, val_trimgalore_R2, keepBX_R2)
            run_command_with_popen(runKeepBXR2)
        # 7. non-dup: 
        if nonDup_int > 0:
            print("\n[START TO GET NON-DUPLICATE READ PAIRS FOR EACH BARCODE]...")
            runNonDup = "/opt/bcGenNonDupFQ/bcGenNonDupFastq {} {} nonDupFq {}".format(keepBX_R1, keepBX_R2, nonDup_int)
            run_command_with_popen(runNonDup)

        print("DONE.\n")
        
    else:
        raise InputPathError("Your Input path is not a folder! Please Check it.")
else:
    raise InputPathError("Your Input folder path is not exists! Please Check it.")
