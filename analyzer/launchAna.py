from subprocess import call
import sys

file_path = sys.argv[1]
file_name = sys.argv[2]

call(["root", "-l", 'run.C("' + file_path + '", "' + file_name + '")'], shell=False)
