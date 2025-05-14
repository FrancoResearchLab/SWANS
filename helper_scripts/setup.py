'''
Author:	E. Reichenberger
Date:	4.30.24

Purpose: under project folder (data/endpoints/project_name/), create sample sub-folders & sub-folders for 10X or cellranger based on user-input in config file
'''

import os
import sys
from sample_list import get_samples

project = sys.argv[1] #project name passed as argument
starting_data = sys.argv[2] #starting data (cellranger or matrix) passed as argument
run_cellranger = sys.argv[3] #starting data (cellranger or matrix) passed as argument

SAMPLE_LIST = []
SAMPLE_LIST = get_samples(project) #get samples from samples.sample_list


for s in SAMPLE_LIST:
	sample_path = 'data'
	path_list = []

	if starting_data == 'fastq' and run_cellranger == 'y':
		path_list = ['endpoints', project, s, 'fastq']
	if starting_data == 'cellranger':
		path_list = ['endpoints', project, s, 'cellranger']
	if starting_data == 'matrix':
		path_list = ['endpoints', project, s, 'matrix']
	
	for p in path_list:
		sample_path = sample_path + '/' + p
		if not os.path.exists(sample_path):
			os.mkdir(sample_path)
