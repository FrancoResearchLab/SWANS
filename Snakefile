'''
Author:	E. Reichenberger
Date:		4.30.2024

Notes: Reference script to call all rules
'''

import sys
import glob
import os
import os.path
import subprocess
from helper_scripts.sample_list import get_samples
from pathlib import Path 
import shutil

configfile: 'configs/prelim_configs.yaml'

# set config parameters
# -------------------------------------------------------------------------
contact = config['contact']
print(contact)
PROJECT = config['PROJECT']
ORGANISM = config['ORGANISM']
STARTING_DATA = config['STARTING_DATA']
RUN_CELLRANGER = config['RUN_CELLRANGER']
RPATH = config['RPATH']
RUN_SOUPX = config['RUN_SOUPX']
RUN_DOUBLETFINDER = config['RUN_DOUBLETFINDER']
MITO = config['MITO']
RIBO = config['RIBO']
MIN_FEATURE_THRESHOLD = config['MIN_FEATURE_THRESHOLD']
MAX_FEATURE_THRESHOLD = config['MAX_FEATURE_THRESHOLD']
SPLIT_LAYERS_BY = config['SPLIT_LAYERS_BY']
COMPONENTS = config['COMPONENTS']
NUM_VARIABLE_FEATURES = config['NUM_VARIABLE_FEATURES']
SCALE_DATA_FEATURES = config['SCALE_DATA_FEATURES']
MITO_REGRESSION = config['MITO_REGRESSION']
RIBO_REGRESSION = config['RIBO_REGRESSION']
CELL_CYCLE_REGRESSION = config['CELL_CYCLE_REGRESSION']
CELL_CYCLE_METHOD = config['CELL_CYCLE_METHOD']
SEURAT_NORMALIZATION_METHOD = config['SEURAT_NORMALIZATION_METHOD']
SEURAT_INTEGRATION_METHOD = config['SEURAT_INTEGRATION_METHOD']
REFERENCE_BASED_INTEGRATION = config['REFERENCE_BASED_INTEGRATION']
RUN_AZIMUTH = config['RUN_AZIMUTH']
RUN_TRANSFERDATA = config['RUN_TRANSFERDATA']
TRANSFERDATA_REF_FILE = config['TRANSFERDATA_REF_FILE']
TRANSFERDATA_REDUCTION = config['TRANSFERDATA_REDUCTION']
TRANSFERDATA_ANNOCOL = config['TRANSFERDATA_ANNOCOL']
RESOLUTION = config['RESOLUTION']
TSNE = config['TSNE']
CONSERVED_GENES = config['CONSERVED_GENES']
STORAGE = config['STORAGE']
THREADS = config['THREADS']
MEMORY = config['MEMORY']

CELLRANGER = config['CELLRANGER']
CELLRANGER_REFERENCE = config['CELLRANGER_REFERENCE']
RUN_MULTIQC = config['RUN_MULTIQC']
MULTIQC = config['MULTIQC']
OUTPUT_BAM = config['OUTPUT_BAM']
RPATH = config['RPATH']
SOUPX_START = config['SOUPX_START']
REFERENCE_SAMPLES = config['REFERENCE_SAMPLES']
AZIMUTH_REFERENCE = config['AZIMUTH_REFERENCE']
USER_GENE_FILE = config['USER_GENE_FILE']
REGRESSION_FILE = config['REGRESSION_FILE']
VISUALIZATION = config['VISUALIZATION']

#-------------------------------------------------------------------------------------

# Set storage type to qs if none provided
#-------------------------------------------------------------------------------------
if STORAGE is None:
	STORAGE = 'qs'
#-------------------------------------------------------------------------------------

# Get samples
#-------------------------------------------------------------------------------------
sample_file = 'samples.sample_list'
if os.path.exists(sample_file) == False:
	print('You must supply a list of files')
	print('Create samples.sample_list first, then run this script again')
	sys.exit()

if os.path.exists(sample_file) == True:
	SAMPLE_LIST = []
	SAMPLE_LIST = get_samples(PROJECT) #get sample list from samples.sample_list 
#-------------------------------------------------------------------------------------

# DEFINITIONS/FUNCTIONS ----------------------------------------------------
# for params with potential for > 1 selection ---------------
def replace_spaces_commas(config_param_name, config_param):
	print(config_param_name)
	normalization_options = ['standard', 'sct']
	integration_options = ['cca', 'harmony', 'rpca']
	storage_options = ['qs', 'rds']
	visualization_options = ['feature', 'violin', 'ridge', 'dot']
	flag = 1

	if type(config_param) == float:
		config_param = str(config_param)

	while ' ' in config_param:
		config_param = config_param.replace(' ', ',')
	config_param = config_param.replace(',,', ',')

	if config_param.endswith(','):
		config_param = config_param.rsplit(',', 1)[0] #remove last comma

	options = []
	if ',' in config_param:
		options = config_param.split(',')
	if ',' not in config_param:
		options = [config_param]
		print(type(options))

	# added 8.26-----
	if config_param_name == 'RESOLUTION' and ',' not in config_param:
		return(config_param)
	#	flag = 0
	# added 8.26-----

	for o in options:
		if config_param_name == 'SEURAT_NORMALIZATION_METHOD':
			if o not in normalization_options:
				print('Please write \'standard\' or \'sct\' in ' + config_param_name + ' the configs/prelim_configs.yaml file, then run this script again')
				sys.exit()
			else:
				flag = 0
				#return(config_param)

		if config_param_name == 'SEURAT_INTEGRATION_METHOD':
			if o not in integration_options:
				print('Please write \'cca\', \'rpca\',  or \'harmony\' for ' + config_param_name + ' in the configs/prelim_configs.yaml file, then run this script again')
				sys.exit()
			else:
				flag = 0

		if config_param_name == 'STORAGE':
			if o not in storage_options:
				print('If you would like to save output as RDS, enter \'rds\' for ' + config_param_name + ' in the configs/prelim_configs.yaml file, then run this script again')
				sys.exit()
			else:
				flag = 0

		if config_param_name == 'VISUALIZATION':
			if USER_GENE_FILE is not None:
				if os.path.isfile(USER_GENE_FILE):
					if o not in visualization_options:
						print('You have selected NO feature, dot, ridge, or violin plots')
						print('and e.g., you typed \'viola\' instead of \'violin\', CTRL+c, and correct your error...')
						#flag = 0 #i think this should be an error instead of passing on the flag
						sys.exit()
					else:
						flag = 0

	if flag == 0:
		return(config_param)

# VISUALIZATION can be blank
# -----------------------------------------------------------------------------------------
if VISUALIZATION is not None and type(VISUALIZATION) != 'NoneType' and VISUALIZATION != ' ':
	VISUALIZATION = replace_spaces_commas('VISUALIZATION', VISUALIZATION)
	print(VISUALIZATION)
# -----------------------------------------------------------------------------------------

# check if null ----------------
def not_null_check_no_yes(config_param, config_param_name):
	numeric_values = ['MITO', 'MIN_FEATURE_THRESHOLD', 'MAX_FEATURE_THRESHOLD', 'COMPONENTS', 'NUM_VARIABLE_FEATURES', 'RESOLUTION', 'THREADS', 'MEMORY']
	do_not_lower = ['RPATH', 'MITO', 'MIN_FEATURE_THRESHOLD', 'MAX_FEATURE_THRESHOLD', 'COMPONENTS', 'NUM_VARIABLE_FEATURES', 'TRANSFERDATA_REF_FILE', 'TRANSFERDATA_REDUCTION', 'TRANSFERDATA_ANNOCOL', 'RESOLUTION', 'THREADS', 'MEMORY']
	multi_options = ['RESOLUTION', 'SEURAT_NORMALIZATION_METHOD', 'SEURAT_INTEGRATION_METHOD', 'STORAGE']
	yes_no = ['RUN_SOUPX', 'RUN_DOUBLETFINDER', 'MITO_REGRESSION', 'CELL_CYCLE_REGRESSION', 'REFERENCE_BASED_INTEGRATION', 'RUN_AZIMUTH', 'RUN_TRANSFERDATA', 'TSNE', 'CONSERVED_GENES']

	#print(config_param, config_param_name)

	if config_param is None or type(config_param) == 'NoneType' or config_param == ' ':
		print(config_param_name + ' cannot be empty.')
		print('Please write a value for ' + config_param_name + ' in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

	if config_param is not None and config_param != ' ':
		#print(config_param, config_param_name)

		if config_param_name not in numeric_values:
			config_param = ''.join(config_param.split()) # remove all spaces

		if config_param_name not in do_not_lower and config_param_name not in numeric_values:
			config_param = config_param.lower()

		if config_param_name in yes_no:
			yes_no_options = ['n', 'y']
			if config_param == 'yes': 
				config_param = 'y'
			if config_param == 'no': 
				config_param = 'n'
			if config_param not in yes_no_options:
				print('You must provide a \'y\' or a \'n\' for', config_param_name)
				print('Please write \'y\' or \'n\' in the configs/prelim_configs.yaml file, then run this script again')
				sys.exit()
		
		if config_param_name in numeric_values:
			if config_param_name != 'RESOLUTION':
				config_param = str(config_param)

				if config_param.isdigit() == False:
					print(config_param_name)
					print('You must provide a numeric value for ', config_param_name)
					print('Please write a numeric value (e.g., 12) for this parameter in the configs/prelim_configs.yaml file, then run this script again')
					sys.exit()

			if config_param_name == 'RESOLUTION':
				if str(config_param).isdecimal() == True:
					config_param = config_param
					print('single', config_param)
				if ',' in str(config_param):
					temp_res = str(config_param).replace(',', '').strip()
					temp_res = temp_res.replace('.', '')

					if temp_res.isdigit() == False:
						print('You must provide a decimal value for ', config_param_name)
						print('Please write decimal value(s) (e.g., 1.2 or 0.1,0.2) for this parameter in the configs/prelim_configs.yaml file, then run this script again')
						sys.exit()

		if config_param_name in multi_options:
			print(config_param, config_param_name)
			config_param = replace_spaces_commas(config_param_name, config_param)
			
		return(config_param)
# -------------------------------------------------------------------------

'''
# -- double check resolution
if str(RESOLUTION).isnumeric() == True: 
	RESOLUTION = RESOLUTION 
elif ',' in str(RESOLUTION): 
	replace_spaces_commas('RESOLUTION', RESOLUTION)
	#temp_res = str(RESOLUTION).replace(',', '').strip() 
	#temp_res = temp_res.replace('.', '') 
	#temp_res = check_numeric('RESOLUTION', temp_res) 
 
#RESOLUTION = RESOLUTION.strip() 
#replace_spaces_commas('RESOLUTION', RESOLUTION)
'''
# -------------------------------------------------------------------------

no_nulls_strings = ['contact', 'PROJECT', 'ORGANISM', 'STARTING_DATA', 'RUN_CELLRANGER', 'RPATH', 'RUN_SOUPX', 'RUN_DOUBLETFINDER', 'MITO', 'MIN_FEATURE_THRESHOLD', 'MAX_FEATURE_THRESHOLD', 'SPLIT_LAYERS_BY', 'COMPONENTS', 'NUM_VARIABLE_FEATURES', 'SCALE_DATA_FEATURES', 'MITO_REGRESSION', 'RIBO_REGRESSION', 'CELL_CYCLE_REGRESSION', 'SEURAT_NORMALIZATION_METHOD', 'SEURAT_INTEGRATION_METHOD', 'REFERENCE_BASED_INTEGRATION', 'RUN_AZIMUTH', 'RUN_TRANSFERDATA', 'RESOLUTION', 'TSNE', 'CONSERVED_GENES', 'STORAGE', 'THREADS']

no_nulls = [contact, PROJECT, ORGANISM, STARTING_DATA, RUN_CELLRANGER, RPATH, RUN_SOUPX, RUN_DOUBLETFINDER, MITO, MIN_FEATURE_THRESHOLD, MAX_FEATURE_THRESHOLD, SPLIT_LAYERS_BY, COMPONENTS, NUM_VARIABLE_FEATURES, SCALE_DATA_FEATURES, MITO_REGRESSION, RIBO_REGRESSION, CELL_CYCLE_REGRESSION, SEURAT_NORMALIZATION_METHOD, SEURAT_INTEGRATION_METHOD, REFERENCE_BASED_INTEGRATION, RUN_AZIMUTH, RUN_TRANSFERDATA, RESOLUTION, TSNE, CONSERVED_GENES, STORAGE, THREADS]

for index,n in enumerate(no_nulls):
	value = no_nulls_strings[index]
	n = not_null_check_no_yes(n, value)

# -- special cases

if RUN_CELLRANGER == 'y':

	if CELLRANGER_REFERENCE is None:
		print('You have selected \'y\' for RUN_CELLRANGER; the path for CELLRANGER_REFERENCE cannot be empty.')
		print('Please list the path to the genome to be used as reference for Cell Ranger in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()
	
	# Re-naming
	REFERENCE = CELLRANGER_REFERENCE

	if CELLRANGER is None:
		print('You have selected \'y\' for RUN_CELLRANGER; the path for CELLRANGER cannot be empty.')
		print('Please list the path to the Cell Ranger executable in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

	if STARTING_DATA != 'fastq':
		print('You have selected \'y\' for RUN_CELLRANGER; the STARTING_DATA for CELLRANGER should be FASTQ files.')
		print('Please make sure data are FASTQ files and that STARTING_DATA is \'fastq\' in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

	if OUTPUT_BAM is None:
		print('You have selected \'y\' for RUN_CELLRANGER; OUTPUT_BAM cannot be empty.')
		print('You must provide a \'y\' or a \'n\' for OUTPUT_BAM')
		sys.exit()

	if RUN_SOUPX == 'y':
		SOUPX_START = 'outs'

if RUN_CELLRANGER == 'y' and (STARTING_DATA == 'fastq' or STARTING_DATA == 'cellranger'):
	if RUN_MULTIQC is None:
		print('You have selected \'y\' for RUN_CELLRANGER or \'cellranger\' for STARTING_DATA: RUN_MULTIQC cannot be empty.')
		print('You must provide a \'y\' or a \'n\' for RUN_MULTIQC')
		sys.exit()

if RUN_MULTIQC == 'y':
	if MULTIQC is None:
		print('You have selected \'y\' for RUN_MULTIQC; MULTIQC cannot be empty.')
		print('Please list the path to the MULTIQC executable in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

	if STARTING_DATA != 'fastq' and STARTING_DATA != 'cellranger':
		print('You have selected \'y\' for RUN_MULTIQC; STARTING_DATA must be \'cellranger\' or RUN_CELLRANGER must \'y\' and STARTING_DATA must be \'fastq\'.')
		print('Please write \'cellranger\' for STARTING_DATA or \'y\' for RUN_CELLRANGER and \'fastq\' for STARTING_DATA in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

if RUN_CELLRANGER == 'n' and RUN_SOUPX == 'y' and STARTING_DATA != 'cellranger':
		print('SoupX requires cellranger output files, please provide cellranger output files and write \'cellranger\' for STARTING_DATA or change RUN_SOUPX to \'n\' in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

if RUN_SOUPX == 'y':
	if (SOUPX_START != 'outs' and SOUPX_START != 'no_clusters' and SOUPX_START != 'h5') or SOUPX_START is None:
		print('Please write \'outs\', \'no_clusters\', or \'h5\' for SOUPX_START in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

if (SPLIT_LAYERS_BY != 'Sample' and SPLIT_LAYERS_BY != 'Experiment'):
	print('Please write \'Sample\' or \'Experiment\' for SPLIT_LAYERS_BY in the configs/prelim_configs.yaml file for SPLIT_LAYERS_BY, then run this script again')
	sys.exit()

if (SCALE_DATA_FEATURES != 'all' and SCALE_DATA_FEATURES != 'variable'):
	print('Please write \'variable\' or \'all\' for SCALE_DATA_FEATURES in the configs/prelim_configs.yaml file, then run this script again')
	sys.exit()

if CELL_CYCLE_REGRESSION == 'y':
	if (CELL_CYCLE_METHOD != 'standard' and CELL_CYCLE_METHOD != 'alternative') or CELL_CYCLE_METHOD is None:
		print('Please write \'standard\' or \'alternative\' for CELL_CYCLE_METHOD in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

# if specific seurat options are being run, check that the additional parameters are present and correct
# if not run and if addtl params are null, configs set to 'does_not_exist' to preserve argument order in R
# -- reference-based integration config checks
if REFERENCE_BASED_INTEGRATION == 'y':
	ref_flag = 0

	if REFERENCE_SAMPLES is None:
		print('You have selected \'y\' for REFERENCE_BASED_INTEGRATION; REFERENCE_SAMPLES cannot be empty.')
		sys.exit()

	# check for sample names in SAMPLE_LIST
	ref_flag = ''
	if ',' not in REFERENCE_SAMPLES:
		if REFERENCE_SAMPLES not in SAMPLE_LIST:
			ref_flag = 1

	if ',' in REFERENCE_SAMPLES:
		ref_samples = REFERENCE_SAMPLES.split(',')
		for r in ref_samples:
			if r not in SAMPLE_LIST:
				ref_flag = 1

	if ref_flag == 1:
		print('You have selected to run REFERENCE_BASED_INTEGRATION but the reference sample(s) you provided does not exist in the sample list:')
		print('Here is the sample list:', SAMPLE_LIST)
		print('Enter a new sample name or names for \'REFERENCE_SAMPLES\' in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

	if ref_flag == 0:
		replace_spaces_commas('REFERENCE_SAMPLES', REFERENCE_SAMPLES)

elif REFERENCE_BASED_INTEGRATION == 'n':
	REFERENCE_SAMPLES = 'does_not_exist'
# --

# -- run azimuth config checks
if RUN_AZIMUTH == 'y':
	if ORGANISM == 'human':
		azimuth_ref_list = ['adiposeref', 'bonemarrowref', 'fetusref', 'heartref', 'humancortexref', 'kidneyref', 'lungref', 'pancreasref', 'pbmcref', 'tonsilref']

	if ORGANISM == 'mouse':
		azimuth_ref_list = ['mousecortexref']

	if (AZIMUTH_REFERENCE not in azimuth_ref_list) or AZIMUTH_REFERENCE is None:
		print('You have select \'y\' for RUN_AZIMUTH; AZIMUTH_REFERENCE must be selected.')
		print('Please write one of the following: (' + ", ".join(str(x) for x in azimuth_ref_list)  + ') for AZIMUTH_REFERENCE in the configs/prelim_configs.yaml file, then run this script again')
		sys.exit()

elif RUN_AZIMUTH == 'n':
	REFERENCE_SAMPLES = 'does_not_exist'
	AZIMUTH_REFERENCE = 'does_not_exist'
# --

# -- run transfer data config checks
if RUN_TRANSFERDATA == 'y':
	if TRANSFERDATA_REF_FILE is None:
		print('You have select \'y\' for RUN_TRANSFERDATA; TRANSFERDATA_REF_FILE must be selected.')
		print('Please provide a path (absolute or relative) to the Seurat object file (.rds or .qs) to be used as the TRANSFERDATA_REF_FILE in the configs/prelim_configs.yaml file, then run this script again.')
		sys.exit()

	if os.path.exists(TRANSFERDATA_REF_FILE) == False:
		print('The file you have provided for TRANSFERDATA_REF_FILE does not exist.')
		print('Please provide a proper name/path for the TRANSFERDATA_REF_FILE file.')
		sys.exit()

	if TRANSFERDATA_REDUCTION is None:
		print('You have select \'y\' for RUN_TRANSFERDATA; TRANSFERDATA_REDUCTION must be selected.')
		print('Please provide the name of the reduction (such as pca) to be used in the TransferData process in the configs/prelim_configs.yaml file, then run this script again.')
		sys.exit()

	if TRANSFERDATA_ANNOCOL is None:
		print('You have select \'y\' for RUN_TRANSFERDATA; TRANSFERDATA_ANNOCOL must be selected.')
		print('Please provide the name of a column in the Seurat reference object metadata to be used for TransferData annotation of the project dataset in the configs/prelim_configs.yaml file, then run this script again.')
		sys.exit()

elif RUN_TRANSFERDATA == 'n':
	TRANSFERDATA_REF_FILE = 'does_not_exist'
	TRANSFERDATA_REDUCTION = 'does_not_exist'
	TRANSFERDATA_ANNOCOL = 'does_not_exist'
# --


# -- handling of null user files
def user_files(user_file, user_file_name):
	if user_file is not None:
		if os.path.exists(user_file) == False:
			out_line = 'The file you have listed in configs/prelim_configs.yaml : ' + user_file_name + ' does not exist.'
			print(out_line)
			print('Please provide a proper name/path for the listed file')
			sys.exit()

	if user_file is None:
		user_file = 'does_not_exist'

	print(user_file)
	return(user_file)

REGRESSION_FILE =  user_files(REGRESSION_FILE, 'REGRESSION_FILE')
USER_GENE_FILE =  user_files(USER_GENE_FILE, 'USER_GENE_FILE')
print(USER_GENE_FILE)

if (ORGANISM != 'mouse' and ORGANISM != 'human') or ORGANISM is None:
	print('This pipeline only works with mouse or human data')
	print('Please write \'human\' or \'mouse\' in the configs/prelim_configs.yaml file, then run this script again')
	sys.exit()

if (STARTING_DATA != 'fastq' and STARTING_DATA != 'cellranger' and STARTING_DATA != 'matrix'):
	print('Please write \'fastq\', \'cellranger\' or \'matrix\' in the configs/prelim_configs.yaml file, then run this script again')
	sys.exit()

'''
# if you want to see all the params at once...
params = [contact, PROJECT, ORGANISM, STARTING_DATA, RUN_CELLRANGER, CELLRANGER, RUN_MULTIQC, MULTIQC, OUTPUT_BAM, RPATH, RUN_SOUPX, SOUPX_START, RUN_DOUBLETFINDER, MITO, RIBO, MIN_FEATURE_THRESHOLD, MAX_FEATURE_THRESHOLD,
SPLIT_LAYERS_BY, COMPONENTS, NUM_VARIABLE_FEATURES, SCALE_DATA_FEATURES, MITO_REGRESSION, RIBO_REGRESSION, CELL_CYCLE_REGRESSION, CELL_CYCLE_METHOD, SEURAT_NORMALIZATION_METHOD, SEURAT_INTEGRATION_METHOD, REFERENCE_BASED_INTEGRATION, REFERENCE_SAMPLES, RUN_AZIMUTH, AZIMUTH_REFERENCE, RESOLUTION, TSNE, CONSERVED_GENES, STORAGE, THREADS, MEMORY, USER_GENE_FILE]

for p in params:
	print(p)
'''


# REDEFINE SEURAT_INTEGRATION_METHOD -------------------------------------------------
if len(SAMPLE_LIST) == 1:
	SEURAT_INTEGRATION_METHOD == 'pca'
#-------------------------------------------------------------------------------------

# DEFINE LOG & BENCHMARK DIRS --------------------------------------------------------
log_directory = 'logs/' + PROJECT + '/'
cellranger_log = log_directory + 'cellranger/'
multiqc_log = log_directory + 'multiqc/'
doubletFinder_log = log_directory + 'doubletFinder/'
soupX_log = log_directory + 'soupX/'
create_initial_seurat_log = log_directory + 'create_initial_seurat/'
qc_report_log = log_directory + 'qc_report/'
memory_log = log_directory + 'memory_log/'
seurat_analysis_log = log_directory + 'seurat_analysis/'
create_images_DGE_log = log_directory + 'create_images_DGE/'

benchmark_log = 'benchmarks/' + PROJECT + '/'
#-------------------------------------------------------------------------------------

# --------------------------------- BIGGIE INPUT -------------------------------------
qc_files = []
cellranger_report = ''
memory_file = 'data/endpoints/' + PROJECT + '/analysis/' + PROJECT + '_memory.txt'
# ------------------------------------------------------------------------------------

# ------------------------------------- CONDITIONS -----------------------------------

# -- CELLRANGER ---------------------------------------------------
qc_cellranger = []
cellranger_output = [] #this will be input for soupx (if y), then doubletFinder (if y)

if RUN_CELLRANGER == 'y':
	qc_cellranger = expand(cellranger_log + PROJECT + '_{s}.stamp', s=SAMPLE_LIST)

# loop through samples, if fastq (new column) is in sample, cellranger can be run, append to list
	for q in qc_cellranger: 
		qc_files.append(q)
		cellranger_output.append(q)

# sometimes unexpected things happen in config w/ true/false
if OUTPUT_BAM == 'y':
	OUTPUT_BAM = 'true'
if OUTPUT_BAM == 'n':
	OUTPUT_BAM = 'false'

# multiqc report
if RUN_MULTIQC == 'y':
	# output
	cellranger_report = 'data/endpoints/' + PROJECT + '/analysis/report/cellranger/' + PROJECT + '_multiqc_cellranger_summary.html'
	cellranger_output.append(cellranger_report)
	qc_files.append(cellranger_report)

# change STARTING_DATA to cellranger, which will be used for everything after this
if RUN_CELLRANGER == 'y':
	STARTING_DATA = 'cellranger'
# -----------------------------------------------------------------

# -- DOUBLETFINDER ------------------------------------------------
doubletfinder_output = []

if RUN_DOUBLETFINDER == 'y':
	# doubletFinder output
	doubletfinder_output = expand('data/endpoints/' + PROJECT + '/' + '{s}/doubletFinder/tables/' + PROJECT + '_{s}_' + 'doublet_ids.txt', s=SAMPLE_LIST)
	for d in doubletfinder_output:
		qc_files.append(d)
# -----------------------------------------------------------------

# ---------- SOUPX ------------------------------------------------
soupX_path = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX', s=SAMPLE_LIST) 
soupX_output = []

if RUN_SOUPX == 'y':
	# soupX output
	barcodes_path_soupX = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX/barcodes.tsv.gz', s=SAMPLE_LIST)
	features_path_soupX = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX/features.tsv.gz', s=SAMPLE_LIST)
	matrix_path_soupX = expand('data/endpoints/' + PROJECT + '/' + '{s}/soupX/matrix.mtx.gz', s=SAMPLE_LIST)

	soupX_list = [barcodes_path_soupX, features_path_soupX, matrix_path_soupX]

	for s in soupX_list: #add each item individually
		for item in s:
			qc_files.append(item)
			soupX_output.append(item)
			#print(item)
# -----------------------------------------------------------------

# -----------------------------------------------------------------

# ---------- MAKE DIRECTORIES PRE-ANALYSIS OUTPUT ------------------------------------
out_parts = ['analysis', 'RDS', 'figures', 'tables', 'report']
string_path = 'data/endpoints/' + PROJECT + '/'

for o in out_parts:
	temp_path = string_path + o
	if o != 'analysis': #want subdirectories under analysis 
		temp_path = string_path + 'analysis/' + o
	if not os.path.exists(temp_path):
		os.mkdir(temp_path)
# ------------------------------------------------------------------------------------

# ----------- DEFINE REPORT DIRECTORIES ----------------------------------------------
report_path = string_path + 'analysis/report/'
report_path_figures = string_path + 'analysis/report/figures/' #pass in analyze sc script
report_path_tables = string_path + 'analysis/report/tables/' #pass in create_images script
compare_images_path = string_path + 'analysis/report/compare_images'

if not os.path.exists(compare_images_path):
	os.mkdir(compare_images_path)
# --------------------------------------------------------------------------------------

# ----------- ADD TO REPORT DIRECTORIES ------------------------------------------------
config_file = 'configs/prelim_configs.yaml'
new_config_file = report_path + 'prelim_configs.yaml'

interactive_report = 'src/rmd/Interactive_report.Rmd'
new_interactive_report = report_path + 'Interactive_report.Rmd'

new_sample_file = report_path + PROJECT + '_' + 'samples.sample_list'
new_user_gene_file = report_path + PROJECT + '_' + 'user_supplied_visualization_genes.txt'
new_user_regression_file = report_path + PROJECT + '_' + 'user_supplied_regression_genes.txt'

# cp config file, samples.sample_list, user gene_lists, & Interactive_report.Rmd
shutil.copyfile(config_file, new_config_file)
shutil.copyfile(interactive_report, new_interactive_report)
shutil.copyfile(sample_file, new_sample_file)

if USER_GENE_FILE != 'does_not_exist':
	shutil.copyfile(USER_GENE_FILE, new_user_gene_file)
if REGRESSION_FILE != 'does_not_exist':
	shutil.copyfile(REGRESSION_FILE, new_user_regression_file)
# -------------------------------------------------------------------------------------

# ----------- INITIAL SEURAT AND QC REPORT --------------------------------------------
#  will save merged seurat object here
storage_file = string_path + 'analysis/RDS/' + PROJECT + '_initial_seurat_object.qs'

if len(SAMPLE_LIST) == 1:
	storage_file = string_path + 'analysis/RDS/' + PROJECT + '_initial_seurat_object.qs'

# the source for the merging of samples will be the STARTING_DATA, unless soupX is run
SEURAT_CREATION_SOURCE = STARTING_DATA
if RUN_SOUPX == 'y':
	SEURAT_CREATION_SOURCE = 'soupX'

# will save qc images here
figure_outs = 'data/endpoints/' + PROJECT + '/analysis/figures/'

# initial_seurat output
f1 = figure_outs + PROJECT + '_qc_1_vln.png'
f2 = figure_outs + PROJECT + '_qc_2_vln.png'

# qc_report input (rmarkdown) and qc_report output (html)
qc_report_rmd = 'src/rmd/qc_report.Rmd'
qc_report_html = string_path + 'analysis/report/qc_report/' + PROJECT + '_qc_report.html'

initial_seurat_list = [f1, f2, storage_file, qc_report_html]
# -------------------------------------------------------------------------------------

# memory
# this would be a good space for 10X vs Flex (conditionally assign folder names)
#-------------------------------------------------------------------------------------
barcodes_files = 'data/endpoints/' + PROJECT + '/*/' + SEURAT_CREATION_SOURCE + '/barcodes.tsv.gz',
features_files = 'data/endpoints/' + PROJECT + '/*/' + SEURAT_CREATION_SOURCE + '/features.tsv.gz',

if SEURAT_CREATION_SOURCE == 'cellranger':
   barcodes_files = 'data/endpoints/' + PROJECT + '/*/' + SEURAT_CREATION_SOURCE + '/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
   features_files = 'data/endpoints/' + PROJECT + '/*/' + SEURAT_CREATION_SOURCE + '/outs/filtered_feature_bc_matrix/features.tsv.gz',
#-------------------------------------------------------------------------------------

# main analyzed seurat object
#-------------------------------------------------------------------------------------
analyzed_seurat_object = 'data/endpoints/' + PROJECT + '/analysis/RDS/' +  PROJECT + '_analyzed_seurat_object.qs'
#-------------------------------------------------------------------------------------

# dge, proportions, dimplots
# -------------------------------------------------------------------------------------
touch_file_create_images_DGE = create_images_DGE_log + PROJECT + '_create_images_DGE_dummy.txt',
#-------------------------------------------------------------------------------------

# List of final files (for testing and troubleshooting purposes)
# -------------------------------------------------------------------------------------
final_files = [qc_files, initial_seurat_list, analyzed_seurat_object, touch_file_create_images_DGE]  #files to be passed to biggie
# -------------------------------------------------------------------------------------

#--------------------TARGET RULES-----------------------------------
# if there is a matching script, it goes in src/scripts/
if RUN_CELLRANGER == 'y': 
	include: 
		"src/rules/cellranger.rules"
if (STARTING_DATA == 'fastq' or STARTING_DATA == 'cellranger') and RUN_MULTIQC == 'y':
	include: 
		"src/rules/multiQC.rules"
include:
	"src/rules/doubletFinder.rules"
include:
	"src/rules/soupX.rules"
include:
	"src/rules/create_initial_seurat.rules"
include:
	"src/rules/get_memory.rules"
include:
	"src/rules/analyze_seurat_object.rules"
include:
	"src/rules/create_images_DGE.rules"
#-------------------------------------------------------------------------------------

#--------------------MESSAGES-----------------------------------
onsuccess:
	print("The main controller pipeline completed with no errors.")
	shell("mail -s 'The main controller pipeline completed with no errors.' "+ config['contact']+" < {log}")

onerror:
	print("The main controller pipeline did not complete without errors."),
	shell("mail -s 'The main controller pipeline did not complete without errors, check the logs and try again.' "+ config['contact']+" < {log}")
#-------------------------------------------------------------------------------------

#--------------------RULES---------------------------------------
rule biggie:
	input:
		# initial_seurat_list,
		# memory_file
		# analyzed_seurat_object
		# final_files
		touch_file_create_images_DGE
#--------------------OUTPUT--------------------------------------
