set -e

# show citations
#-----------------------------------------------------------------------------
sh helper_scripts/citations.sh
#-----------------------------------------------------------------------------

# confirm sample list file exists
#-----------------------------------------------------------------------------
sample_list_file="samples.sample_list"

if [ ! -e "$sample_list_file" ]; then
	echo "------------------------------"
	echo "ACHTUNG ACHTUNG ACHTUNG"
	echo "------------------------------"
	echo "$sample_list_file does not exist."
	echo ""
	echo "Create $sample_list_file according to the instructions & run this script again."
	exit 1;
fi
#-----------------------------------------------------------------------------

# confirm config file exists
#-----------------------------------------------------------------------------
config_file="configs/prelim_configs.yaml"

if [ ! -e "$config_file" ]; then
	echo "------------------------------"
	echo "ACHTUNG ACHTUNG ACHTUNG"
	echo "------------------------------"
	echo "$config_file does not exist."
	echo ""
	echo "Create $config_file according to the instructions & run this script again."
	exit 1;
fi
#-----------------------------------------------------------------------------

# get project name, create dirs, create symbolic links
#-----------------------------------------------------------------------------
echo "Setting up project directory"

path="data/endpoints/"
project=`python3 helper_scripts/cache.py PROJECT:` #retrieves project name from config file

echo "Creating project directory (if it does not exist)"
path=$path$project

starting_data=`python3 helper_scripts/cache.py STARTING_DATA:` #retrieves starting data from config file
run_cellranger=`python3 helper_scripts/cache.py RUN_CELLRANGER:` #retrieves starting data from config file

mkdir -p $path
python3 helper_scripts/setup.py $project $starting_data $run_cellranger #makes project_name/sample_name/[matrix]or[cellranger] for each sample

python3 helper_scripts/make_symbolic_links.py $project $starting_data $run_cellranger #this will make a bash script
sh helper_scripts/make_symbolic_links.sh #this will make the symbolic links in each matrix directory

threads=`python3 helper_scripts/cache.py THREADS:` #retrieves project name from config file
rpath=`python3 helper_scripts/cache.py RPATH:` #retrieves rpath name from config file

# make .Renviron file to control libpaths
echo "R_LIBS_USER=$rpath" > .Renviron
#-----------------------------------------------------------------------------

# call Snakemake 
#-----------------------------------------------------------------------------
# snakemake --snakefile Snakefile --printshellcmds --dryrun 
# snakemake --snakefile Snakefile --printshellcmds --dryrun --rerun-triggers mtime
snakemake --cores $threads --snakefile Snakefile --printshellcmds 
#-----------------------------------------------------------------------------

# show citations again
#-----------------------------------------------------------------------------
sh helper_scripts/citations.sh
#-----------------------------------------------------------------------------

#check if final config yaml exists
#-----------------------------------------------------------------------------
final_config_file="configs/post_annotation_configs.yaml"
run_final=`python3 helper_scripts/cache_final.py RUN_FINAL_ANALYSIS:` #retrieves starting data from config file

if [ -e "$final_config_file" ] && [[ $run_final == "y" ]]; then
	echo "You have the final config file, let the magic begin"
	# snakemake --snakefile FinalSnakefile --printshellcmds --dryrun
	snakemake --cores $threads --snakefile FinalSnakefile --printshellcmds #--forceall 
fi
#-----------------------------------------------------------------------------
