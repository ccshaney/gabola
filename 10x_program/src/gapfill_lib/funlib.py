def proper_threads(input_num):
    import os
    # Get total CPU number:
    cpu = os.cpu_count()
    if input_num > cpu:
        print("[WARNING] Your input CPU is {}, but you only have {} CPUs.".format(input_num, cpu))
    if input_num < cpu and input_num != 0:
        cpu = input_num
    print("[RUN CPUs]... {}".format(cpu))
    return cpu

def run_command_with_popen_communicate(commandStr):
    import subprocess
    print('\n{}\n'.format(commandStr))
    p = subprocess.Popen(commandStr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, error) = p.communicate()
    assert p.poll() == 0, "[ERROR] at step: {}\nreturn code is {}, {}".format(commandStr, p.poll(), error.decode("utf-8"))
    return stdout.decode("utf-8")

def run_command_with_popen_communicate_onlyReturn(commandStr):
    import subprocess
    p = subprocess.Popen(commandStr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, error) = p.communicate()
    assert p.poll() == 0, "[ERROR] at step: {}\nreturn code is {}, {}".format(commandStr, p.poll(), error.decode("utf-8"))
    return stdout.decode("utf-8")

# def run_command_with_popen_readline_withReturn(commandStr):
#     import subprocess
#     import sys
#     p = subprocess.Popen(commandStr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     while p.poll() == None:
#         buff = p.stdout.readline()
#         sys.stdout.write(buff.decode('utf-8'))
#         sys.stdout.flush()
#     (stdout, error) = p.communicate()
#     assert p.poll() == 0, "[ERROR] at step: {}\nReturn code is {}, stderr = {}".format(commandStr, p.poll(), error.decode("utf-8"))
#     return stdout.decode("utf-8")
#
# def run_command_with_popen_readline(commandStr):
#     import subprocess
#     import sys
#     p = subprocess.Popen(commandStr, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     while p.poll() == None:
#         buff = p.stdout.readline()
#         sys.stdout.write(buff.decode('utf-8'))
#         sys.stdout.flush()
#     (stdout, error) = p.communicate()
#     assert p.poll() == 0, "[ERROR] at step: {}\nReturn code is {}, stderr = {}".format(commandStr, p.poll(), error.decode("utf-8"))

def run_command_with_subprocessRUN_forParallel(commStr):
    import subprocess
    print(commStr)
    subprocess.run(commStr, shell=True)

def print_sameline(string, flushTime):
    import sys
    import time
    sys.stdout.write(string)
    print('\r', end='')
    sys.stdout.flush()
    time.sleep(flushTime)

