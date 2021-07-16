import argparse
import subprocess
import os
import re
import sys
import time
from src.gapfill_preprocess_10xdata.process_path import getFileNameFromAbsPath, splitFileSuffixThenAddNewSuffix
from src.gapfill_preprocess_10xdata.funlib import proper_threads, run_command_with_popen, run_popen_with_return
from src.gapfill_preprocess_10xdata.funlib import run_command_with_subprocessRUN_forParallel

# #############
# ERROR class:
# #############

class InputValueError(Exception):
    pass


# #########
# Fuction:
# #########



# #########
# Program:
# #########

version_number = '0.0.1'
parser = argparse.ArgumentParser(
    description='''Select qualify alignment from sam file which is output from BWA or Kart.
    First, filter proper pair by SAMBAMBA (v0.6.9). Then, filter by CIGAR match percentage and MAPQ score.\n
    # proper pair define:
     not (unmapped or mate_is_unmapped or secondary_alignment or failed_quality_control or duplicate or 
     supplementary or chimeric)''',
    epilog='''Any questions or comments, please contact Fu Po-Ying <spashleyfu@gmail.com> 
    filter_samfile.py version: {}'''.format(version_number))
parser.add_argument('-s', '--sam', dest="sam", required=True, type=str,
                    help='un-filter SAM file path.')
parser.add_argument('-c', '--CIGAR', dest="cigar_map_quality", default=70, type=int,
                    help='''CIGAR match quality: count([M=X])/len(seq)*100. 
                    From 0 to 100. [default: 70]''')
parser.add_argument('-m', '--MAPQ', dest="mapq", type=int, default=60,
                    help='''MAPQ score. From 0 to 60. [default: 60]''')
parser.add_argument('-t', '--nthreads', dest="nthreads", type=int, default=1,
                    help='Number of threads. Input 0 will use all threads.(default: 1)')

options = parser.parse_args()
input_samfile = options.sam
cigar_in = int(options.cigar_map_quality)
mapQ_in = int(options.mapq)
cores = int(options.nthreads)

# Check Input parameter:
# -t:
cpu = proper_threads(cores)
# -c:
if cigar_in <= 0 or cigar_in > 100:
    raise InputValueError("""
    Input CIGAR match percentage ERROR! 
    Value range from 1 to 100. MIN = 1, MAX = 100, unsigned integer.""")
# -m:
if mapQ_in < 0 or mapQ_in > 60:
    raise InputValueError("""
    Input MAPQ score ERROR! 
    Value range from 0 to 60. MIN = 0, MAX = 60, unsigned integer.""")
# -s, sam file path:
if os.path.islink(input_samfile) == True:
    # get folder path of link file: (output file HERE!)
    soft_input_sam = os.path.abspath(input_samfile)
    samfile_folder = os.path.dirname(soft_input_sam) + '/'
    # get actual file path:
    input_samfile = os.path.realpath(input_samfile)
    print("INPUT SAM FILE IS A LINK FILE...\nACTURE FILE at {}, BUT OUTPUT FILE will at {}\n".format(input_samfile,
                                                                                                   samfile_folder))
else:
    input_samfile = os.path.abspath(input_samfile)
    samfile_folder = os.path.dirname(input_samfile) + '/'
    print("INPUT SAM FILE at {}, AND OUTPUT FILE will at = {}\n".format(input_samfile, samfile_folder))

# Generate OUTPUT file path: (at input_samfile's folder)
# 1) Proper pair:
out_proper_pair_sam_path = samfile_folder + \
                           splitFileSuffixThenAddNewSuffix(getFileNameFromAbsPath(input_samfile), '_ProperPair.sam')
print("[PROPER PAIR SAM FILE]... {}".format(out_proper_pair_sam_path))

# 2) filter CIGAR + MAPQ:
out_CIGAR_MAPQ_sam_path = samfile_folder + \
                          splitFileSuffixThenAddNewSuffix(getFileNameFromAbsPath(input_samfile),
                                                          '_C{}M{}.sam'.format(cigar_in, mapQ_in))
print("[QUALIFY CIGAR AND MAPQ]... = {}".format(out_CIGAR_MAPQ_sam_path))

# Get program absolute path:
prog_path = sys.argv[0]
prog_path = os.path.abspath(prog_path)
prog_folder = os.path.dirname(prog_path) + '/'
prog_folder_src = prog_folder + 'src/'
print("prog_folder_src = {}".format(prog_folder_src))

# ###### Process Proper pair: Use sambamba ######
print("\n[PROCESS PROPER PAIR]...")
runProperPair = "sambamba view -F 'not (unmapped or mate_is_unmapped or " \
                "secondary_alignment or failed_quality_control or duplicate " \
                "or supplementary or chimeric)' -S {} -o {} -t {}".format(input_samfile, out_proper_pair_sam_path, cpu)
print(runProperPair)
print(run_popen_with_return(runProperPair))

# print output file info:
getFileinfo = "ls -lh {}".format(out_proper_pair_sam_path)
run_command_with_popen(getFileinfo)

# ###### split file dependent on CPU number: ######
if cpu == 1:
    print("\n[PROCESS CIGAR{} MAPQ{}]...".format(cigar_in, mapQ_in))
    runCIGAR_MAPQ_1cpu = "python {}filter_CIGAR_MAPQ.py -s {} -o {} -c {} -m {}".format(prog_folder_src,
                                                                                        out_proper_pair_sam_path,
                                                                                        out_CIGAR_MAPQ_sam_path,
                                                                                        cigar_in, mapQ_in)
    run_command_with_popen(runCIGAR_MAPQ_1cpu)
    # print output file info:
    getFileinfo = "ls -lh {}".format(out_CIGAR_MAPQ_sam_path)
    run_command_with_popen(getFileinfo)
else:
    print("\n[PROCESS CIGAR{} MAPQ{}]...".format(cigar_in, mapQ_in))
    outSplitPrefix = samfile_folder + 'tmpSplit_'
    splitNPart = cpu if cpu % 2 == 0 else (cpu-1)
    print("[FOR MULTI-THREADS] SPLIT PROPER SAM FILE INTO {} PART...".format(splitNPart))
    runSplit = "sh {}splitFile_bylines.sh {} {} {}".format(prog_folder_src,
                                                           out_proper_pair_sam_path,
                                                           splitNPart,
                                                           outSplitPrefix)
    run_command_with_popen(runSplit)
    # print output file info:
    getFileinfo = "ls -lh {}* | head".format(outSplitPrefix)
    run_command_with_popen(getFileinfo)

    # ###### RUN split file separately: ######
    check = 1
    CIGAR_MAPQ_PAET_suffix = "_C{}M{}.sam".format(cigar_in, mapQ_in)
    proper_sam_split_match = "^tmpSplit_[0-9]{2,4}$"
    proper_sam_split_prog = re.compile(proper_sam_split_match)
    for entry in os.scandir(samfile_folder):
        if not entry.name.startswith('.') and entry.is_file() and proper_sam_split_prog.match(entry.name):
            print("++++++> [PROCESS] ", entry.name)
            # RUN in parallel:
            runFilterParallel = "nohup python {}filter_CIGAR_MAPQ.py " \
                                "-s {} -o {} -c {} -m {} > {}.log 2> {}.err &".format(prog_folder_src,
                                                                                      samfile_folder + entry.name,
                                                                                      samfile_folder + entry.name + CIGAR_MAPQ_PAET_suffix,
                                                                                      cigar_in, mapQ_in,
                                                                                      samfile_folder + entry.name,
                                                                                      samfile_folder + entry.name)
            run_command_with_subprocessRUN_forParallel(runFilterParallel)
    # Check all process, need to wait until all done:
    print("process in parallel, please wait...\n")
    while check > 0:
        run_command_with_popen("ps aux | grep -G 'python .*src/filter_CIGAR_MAPQ\.py'")
        check = int(run_popen_with_return("ps aux | grep -G 'python {}filter_CIGAR_MAPQ\.py' | grep -v 'grep' | wc -l".format(prog_folder_src)))
        print("{} program still running...sleep 2 sec".format(check))
        time.sleep(2)
    print(run_popen_with_return("ls {}*_C{}M{}.sam".format(outSplitPrefix, cigar_in, mapQ_in)))

    # ###### Combind files into ONE file: ######
    print("\n[COMBIND ALL SPLIT FILES INTO ONE FILE]...")
    if os.path.isfile(out_CIGAR_MAPQ_sam_path):
        print("\t{} already exist, remove it!\n".format(out_CIGAR_MAPQ_sam_path))
        os.remove(out_CIGAR_MAPQ_sam_path)
    CIGARMAPQ_sam_split_match = "^tmpSplit_[0-9]{2,4}" + CIGAR_MAPQ_PAET_suffix #"_C{}M{}.sam".format(cigar_in, mapQ_in)
    print("CIGARMAPQ_sam_split_match = ", CIGARMAPQ_sam_split_match)
    CIGARMAPQ_split_prog = re.compile(CIGARMAPQ_sam_split_match)
    # list split files:
    for entry in os.scandir(samfile_folder):
        print("entry = ", entry)
        print("CIGARMAPQ_split_prog.match(entry.name) = ", CIGARMAPQ_split_prog.match(entry.name))
        if not entry.name.startswith('.') and entry.is_file() and CIGARMAPQ_split_prog.match(entry.name):
            print("entry.name = ", entry.name)
            runCatFile = "cat {} >> {}".format(samfile_folder + entry.name, out_CIGAR_MAPQ_sam_path)
            run_command_with_popen(runCatFile)
    run_command_with_popen("ls -lh {}".format(out_CIGAR_MAPQ_sam_path))
    # remove files:
    print("\n[REMOVE TMP FILE]...\n\t{}".format(subprocess.run("rm -rf {}".format(samfile_folder+"tmpSplit_*"), shell=True)))

print("\nCOMPLETED SUCCESSFULLY!")

