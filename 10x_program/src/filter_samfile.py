import argparse
import subprocess
import os
import re
import sys
from gapfill_lib.process_path import getFileNameFromAbsPath, splitFileSuffixThenAddNewSuffix
from gapfill_lib.funlib import proper_threads, run_command_with_popen_communicate
from gapfill_lib.funlib import run_command_with_subprocessRUN_forParallel, print_sameline
from gapfill_lib.funlib import run_command_with_popen_communicate_onlyReturn

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
assert cigar_in > 0, "CIGAR MAPPING QUALITY cannot small than zero! ( -c > 0 )"
assert cigar_in <= 100, "CIGAR MAPPING QUALITY cannot large than 100! ( -c <= 100 )"

# -m:
assert mapQ_in >= 0, "MAPQ cannot small than zero! ( -m >= 0 )"
assert mapQ_in <= 60, "MAPQ cannot large than 60! ( -m <= 60 )"

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
prog_folder = os.path.dirname(prog_path) + '/' # NOW at /opt/10x_program/src/
# print("prog_folder = {}".format(prog_folder))

# ###### Process Proper pair: Use sambamba ######
print("\n[PROCESS PROPER PAIR]...")
runProperPair = "sambamba view -F 'not (unmapped or mate_is_unmapped or " \
                "secondary_alignment or failed_quality_control or duplicate " \
                "or supplementary or chimeric)' -S {} -o {} -t {}".format(input_samfile, out_proper_pair_sam_path, cpu)
print(runProperPair)
run_command_with_popen_communicate(runProperPair)

# print output file info:
getFileinfo = "ls -lh {}".format(out_proper_pair_sam_path)
print(run_command_with_popen_communicate(getFileinfo))

# ###### split file dependent on CPU number: ######
if cpu == 1:
    print("\n[PROCESS CIGAR{} MAPQ{}]...".format(cigar_in, mapQ_in))
    runCIGAR_MAPQ_1cpu = "python {}filter_CIGAR_MAPQ.py -s {} -o {} -c {} -m {}".format(prog_folder,
                                                                                        out_proper_pair_sam_path,
                                                                                        out_CIGAR_MAPQ_sam_path,
                                                                                        cigar_in, mapQ_in)
    run_command_with_popen_communicate(runCIGAR_MAPQ_1cpu)
    # print output file info:
    getFileinfo = "ls -lh {}".format(out_CIGAR_MAPQ_sam_path)
    print(run_command_with_popen_communicate(getFileinfo))
else:
    print("\n[PROCESS CIGAR{} MAPQ{}]...".format(cigar_in, mapQ_in))
    outSplitPrefix = samfile_folder + 'tmpSplit_'
    splitNPart = cpu if cpu % 2 == 0 else (cpu-1)
    print("[FOR MULTI-THREADS] SPLIT PROPER SAM FILE INTO {} PART...".format(splitNPart))
    runSplit = "sh {}splitFile_bylines.sh {} {} {}".format(prog_folder,
                                                           out_proper_pair_sam_path,
                                                           splitNPart,
                                                           outSplitPrefix)
    run_command_with_popen_communicate(runSplit)
    # print output file info:
    getFileinfo = "ls -lh {}* | head".format(outSplitPrefix)
    print(run_command_with_popen_communicate(getFileinfo))

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
                                "-s {} -o {} -c {} -m {} > {}.log 2> {}.err &".format(prog_folder,
                                                                                      samfile_folder + entry.name,
                                                                                      samfile_folder + entry.name + CIGAR_MAPQ_PAET_suffix,
                                                                                      cigar_in, mapQ_in,
                                                                                      samfile_folder + entry.name,
                                                                                      samfile_folder + entry.name)
            run_command_with_subprocessRUN_forParallel(runFilterParallel)
    # Check all process, need to wait until all done:
    print("process in parallel, please wait...\n")
    while check > 0:
        # run_command_with_popen("ps aux | grep -G 'python .*src/filter_CIGAR_MAPQ\.py'")
        check = int(run_command_with_popen_communicate_onlyReturn("ps aux | grep -G 'python {}filter_CIGAR_MAPQ\.py' | grep -v 'grep' | wc -l".format(prog_folder)))
        strOut = "{} program still running...".format(check)
        print_sameline(strOut, 0.01)
    print(run_command_with_popen_communicate("ls -lh {}*_C{}M{}.sam".format(outSplitPrefix, cigar_in, mapQ_in)))

    # ###### Combind files into ONE file: ######
    print("\n[COMBIND ALL SPLIT FILES INTO ONE FILE]...")
    if os.path.isfile(out_CIGAR_MAPQ_sam_path):
        print("\t{} already exist, remove it!\n".format(out_CIGAR_MAPQ_sam_path))
        os.remove(out_CIGAR_MAPQ_sam_path)
    CIGARMAPQ_sam_split_match = "^tmpSplit_[0-9]{2,4}" + CIGAR_MAPQ_PAET_suffix #"_C{}M{}.sam".format(cigar_in, mapQ_in)
    # print("CIGARMAPQ_sam_split_match = ", CIGARMAPQ_sam_split_match)
    CIGARMAPQ_split_prog = re.compile(CIGARMAPQ_sam_split_match)
    # list split files:
    for entry in os.scandir(samfile_folder):
        # print("entry = ", entry)
        # print("CIGARMAPQ_split_prog.match(entry.name) = ", CIGARMAPQ_split_prog.match(entry.name))
        if not entry.name.startswith('.') and entry.is_file() and CIGARMAPQ_split_prog.match(entry.name):
            # print("entry.name = ", entry.name)
            runCatFile = "cat {} >> {}".format(samfile_folder + entry.name, out_CIGAR_MAPQ_sam_path)
            run_command_with_popen_communicate(runCatFile)
    print(run_command_with_popen_communicate("ls -lh {}".format(out_CIGAR_MAPQ_sam_path)))
    # remove files:
    print("\n[REMOVE TMP FILE]...\n\t{}".format(subprocess.run("rm -rf {}".format(samfile_folder+"tmpSplit_*"), shell=True)))

print("\nCOMPLETED SUCCESSFULLY!")
