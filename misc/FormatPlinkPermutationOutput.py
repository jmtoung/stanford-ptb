#!/usr/bin/env python

"""
    Format the output from plink permutation runs

    Command Line Usage:
        FormatPlinkPermutationOutput.py [<options>...] 

    Options:
        --perm-dir 	Directories containing the output of plink permutation runs
	--best-or-all	Indicate whether plink was run in '--mperm-save' (specify 'best') or '--mperm-all' (specify 'all')
	--max-perm	Indicate the maximum number of permutations you want to analyze
	--min-perm	Indicate the minimum number of permutations you want to analyze
	--output	Output file that contains all plink permutation outputs
"""

import pdb
import os
import sys
import argparse
import types
import re
import logging
import subprocess
from collections import defaultdict
import random

#---- global data
logging.basicConfig(level=logging.DEBUG, format='%(levelname)s | %(asctime)s | %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

#---- internal support stuff

#---- module API

def formatPlinkPermutationOutput(directories,best_or_all,max_perm,min_perm,output):
    
    # loop through directories and find files
    perm_files = []
    abs_dir_regex = re.compile("^\/")
    for directory in directories:
	if(not abs_dir_regex.match(directory)):
	    raise Exception("please give full path for directory -> " + directory)
    	
	# look up all the files
	for dirname, dirnames, filenames in os.walk(directory):
	    for subdirname in dirnames:
		if (re.match("^[0-9]+-[0-9]+$", subdirname)):
		    logging.info("\t".join(["found sub perm dir",subdirname,"in",dirname]))

		for dirname2, dirnames2, filenames2 in os.walk(os.path.join(dirname,subdirname)):
		    for filename2 in filenames2:
			filename2_info = os.stat(os.path.join(dirname2, filename2))
			if (filename2_info.st_size == 0):
			    continue
			if ((best_or_all == "best") & (filename2 == "plink.mperm.dump.best")) | ((best_or_all == "all") & (filename2 =="plink.mperm.dump.all")):
				perm_files.append(os.path.join(dirname2,filename2))
				logging.info("\t".join(["adding mperm file",os.path.join(dirname2,filename2)]))
    
    # trim off files based on the number of permutations we need
    num_permutations_count_total = 0
    perm_files_trim = []
    for perm_file in perm_files:
	num_permutations_count = 0

	# if the total number of permutations we've seen thus far is greater than max_perm, then break
	if num_permutations_count_total >= max_perm:
	    break
	perm_file_fh = open(perm_file)
	# subtract 1 for the header file
	num_permutations_count_total += -1
	num_permutations_count += -1
	for line in perm_file_fh:
	    num_permutations_count += 1
	    num_permutations_count_total += 1
	
	perm_files_trim.append(perm_file)
	logging.info("saw " + str(num_permutations_count) + " permutations for " + perm_file)

    logging.info("trimmed " + str((len(perm_files) - len(perm_files_trim))) + " files")
    perm_files = perm_files_trim
	
    # check that a minimum number of permutations have been found
    if num_permutations_count_total < min_perm:
	raise Exception("we only saw " + str(num_permutations_count_total) + " but we require a minimum of " + str(min_perm))
    
    # cat all files into one file
    max_cat_at_once = 20
    
    files_to_cat = []
    while(perm_files):
	idx_last = 0
	if len(perm_files) >= max_cat_at_once:
	    idx_last = max_cat_at_once
        else:
	    idx_last = len(perm_files)
	files_to_cat.extend(perm_files[0:idx_last])
	del(perm_files[0:idx_last])
	_cat_files(files_to_cat,output)
	files_to_cat = [output]
    
    # open output file and get rid of repetitive 0 line
    output_temp = output + ".temp"
    output_temp_fh = open(output_temp,'w')
    output_fh = open(output)
    line_0 = []
    perm_count = 0
    for line in output_fh:
	line = line.rstrip('\n\r').split()
	if line[0] == "0":
	    if not line_0:
		line_0 = line
		output_temp_fh.write("\t".join(line) + "\n")
	    else:
		if len(line) != len(line_0):
		    raise Exception("line 0 don't match")
		for i in range(0,len(line)):
		    if line_0[i] != line[i]:
			raise Exception("value idx " + i + "doesn't match")
	else:
	    perm_count += 1
	    if (perm_count > max_perm):
	        break
	    output_temp_fh.write("\t".join(line) + "\n")
    
    _runCommand("mv " + output_temp + " " + output)
    
        
	    
def _cat_files(files_to_cat,output):
    
    temp = output + ".temp"
    
    command = "cat " + " ".join(files_to_cat) + " > " + temp
    
    _runCommand(command)
    
    command = "mv " + temp + " " + output
    
    _runCommand(command)
    
def _runCommand(command):
    logging.info(command)
    try:
	subprocess.call(command, shell=True)
    except subprocess.CalledProcessError:
	logging.error("error: " + command)

#---- mainline

def main(argv):

    parser = argparse.ArgumentParser(description='Reformats output from `plink --recode` command for plotting in R')
    
    parser.add_argument('--perm-dir', required=True, nargs='+', help='directories containing the plink permutation output')
    parser.add_argument('--best-or-all', required=True, help='whether plink was run in --mperm-save or --mperm-save-all mode')
    parser.add_argument('--max-perm', required=True, type=int, help='maximum number of permutations to be analyzed')
    parser.add_argument('--min-perm', required=True, type=int, help='minimum number of permutations to be analyzed')
    parser.add_argument('--output', required=True, help='minimum number of permutations to be analyzed')

    args = parser.parse_args()
    
    if args.best_or_all == "best":
	pass
    elif args.best_or_all == "all":
	pass
    else:
	raise Exception("--best-or-all must equal 'best' or 'all'")
    
    logging.info("formatPlinkPermutationOuput(%r,%r,%r,%r,%r)", args.perm_dir, args.best_or_all, args.max_perm, args.min_perm, args.output)
    formatPlinkPermutationOutput(args.perm_dir, args.best_or_all, args.max_perm, args.min_perm, args.output)
        
if __name__ == "__main__":
    sys.exit(main(sys.argv))



