set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# show citations
#-----------------------------------------------------------------------------
sh $SCRIPT_DIR/helper_scripts/citations.sh
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
project=`python3 $SCRIPT_DIR/helper_scripts/cache.py PROJECT:` #retrieves project name from config file

echo "Creating project directory (if it does not exist)"
path=$path$project

starting_data=`python3 $SCRIPT_DIR/helper_scripts/cache.py STARTING_DATA:` #retrieves starting data from config file
run_cellranger=`python3 $SCRIPT_DIR/helper_scripts/cache.py RUN_CELLRANGER:` #retrieves cellranger y/n from config file

mkdir -p $path
python3 $SCRIPT_DIR/helper_scripts/setup.py $project $starting_data $run_cellranger #makes project_name/sample_name/[matrix]or[cellranger] for each sample

python3 $SCRIPT_DIR/helper_scripts/make_symbolic_links.py $project $starting_data $run_cellranger #this will make a bash script
sh $SCRIPT_DIR/helper_scripts/make_symbolic_links.sh #this will make the symbolic links in each matrix directory

threads=`python3 $SCRIPT_DIR/helper_scripts/cache.py THREADS:` #retrieves project name from config file
rpath=`python3 $SCRIPT_DIR/helper_scripts/cache.py RPATH:` #retrieves rpath name from config file

# make .Renviron file to control libpaths
echo "R_LIBS_USER=$rpath" > .Renviron
#-----------------------------------------------------------------------------

# Get paths for Singularity bind mounts (preliminary analysis)
#-----------------------------------------------------------------------------
prelim_bind_mnts=""

# Gets paths from sample file to use as bind mounts to access sample data
echo -e "\n\n============= Sample file directories for Singularity bind mounts ============="
sample_bind_mnts=$(python3 $SCRIPT_DIR/helper_scripts/get_sample_paths.py)
prelim_bind_mnts+="$sample_bind_mnts"
echo -e "===================================================================================\n\n"

echo -e "========== CELLRANGER_REFERENCE directory for Singularity bind mounting ==========="
# Gets cellranger reference genome dir from config file, checks existtence, adds as bind mount
if [ "$run_cellranger" = "y" ]; then
	cellranger_reference=`python3 $SCRIPT_DIR/helper_scripts/cache.py CELLRANGER_REFERENCE:`
	if [[ -n "$cellranger_reference" && -d "$cellranger_reference" ]]; then
		echo -e "[PASS] CELLRANGER_REFERENCE path exists, directory will be bound -> $cellranger_reference"
		prelim_bind_mnts+=",$cellranger_reference"
	else
		echo -e "[WARN] CELLRANGER_REFERENCE path DOES NOT EXIST -> $cellranger_reference"
	fi
	else
		echo -e "[WARN] No directory path entered for CELLRANGER_REFERENCE"
fi
echo -e "===================================================================================\n\n"

echo -e "============= Preliminary analysis files for Singularity bind mounts =============="
# Check that any files in the prelim_configs.yaml exists, gets their directory for bind mounting
prelim_file_configs=("TRANSFERDATA_REF_FILE" "REGRESSION_FILE" "USER_GENE_FILE")
for config in "${prelim_file_configs[@]}"; do
	# Get the file path by calling cache.py with the config name
	file_path=$(python3 "$SCRIPT_DIR/helper_scripts/cache.py" "${config}:")

	if [[ -z "$file_path" ]]; then
  		echo -e "[WARN] No file path entered for $config"
  		continue
	fi

	# Get the absolute directory path of that file
	dir=$(realpath "$(dirname "$file_path")")

	# Check dir exsists
	if [[ -n "$dir" && -d "$dir" ]]; then
		[[ -n "$prelim_bind_mnts" ]] && prelim_bind_mnts+=","
		prelim_bind_mnts+="$dir"
		echo -e "[PASS] File exists for $config: $file_path, directory will be bound -> $dir"
	else
    	echo -e "Warning: Directory for $config does not exist: $config -> $dir"
  	fi
done
echo -e "===================================================================================\n\n"

# Remove duplicates
echo -e "================= Singularity bind mounts for preliminary analysis ================="
prelim_bind_mnts=$(echo "$prelim_bind_mnts" | tr ',' '\n' | awk '!seen[$0]++' | paste -sd ',' -)
echo -e "$prelim_bind_mnts"
echo -e "====================================================================================\n\n\n"
#-----------------------------------------------------------------------------

# call Snakemake (sans Singularity)
#-----------------------------------------------------------------------------
# snakemake --snakefile Snakefile --printshellcmds --dryrun 
# snakemake --snakefile Snakefile --printshellcmds --dryrun --rerun-triggers mtime
# snakemake --cores $threads --snakefile $SCRIPT_DIR/Snakefile
#-----------------------------------------------------------------------------

# call Snakemake with Singularity
#-----------------------------------------------------------------------------
snakemake --snakefile $SCRIPT_DIR/Snakefile \
	--cores $threads \
	--printshellcmds \
	--use-singularity \
	--singularity-args "-B $prelim_bind_mnts"
#-----------------------------------------------------------------------------

# show citations again
#-----------------------------------------------------------------------------
sh $SCRIPT_DIR/helper_scripts/citations.sh
#-----------------------------------------------------------------------------

#check if final config yaml exists
#-----------------------------------------------------------------------------
final_config_file="configs/post_annotation_configs.yaml"
run_final=`python3 $SCRIPT_DIR/helper_scripts/cache_final.py RUN_FINAL_ANALYSIS:` #retrieves starting data from config file

if [ -e "$final_config_file" ] && [[ $run_final == "y" ]]; then
	echo "You have the final config file, let the magic begin."

	# Get paths for Singularity bind mounts (post-annotation analysis)
	#-----------------------------------------------------------------------------
	postanno_bind_mnts=""

	echo -e "============= Post-annotation analysis files for Singularity bind mounts =============="
	# Gets cellranger reference genome dir from config file
	postanno_file_configs=("CLUSTER_ANNOTATION_FILE" "USER_ANALYZED_SEURAT_OBJECT" "FINAL_USER_GENE_FILE")

	for config in "${postanno_file_configs[@]}"; do
		# Get the file path by calling cache.py with the config name
		file_path=$(python3 "$SCRIPT_DIR/helper_scripts/cache_final.py" "${config}:")

		if [[ -z "$file_path" ]]; then
			echo "[WARN] No file path entered for $config"
			continue
		fi

		# Get the absolute directory path of that file
		dir=$(realpath "$(dirname "$file_path")")

		# Check dir exsists
		if [[ -n "$dir" && -d "$dir" ]]; then
			[[ -n "$postanno_bind_mnts" ]] && postanno_bind_mnts+=","
			postanno_bind_mnts+="$dir"
			echo "[PASS] File exists for $config: $file_path, directory will be bound -> $dir"
		else
			echo "Warning: Directory for $config does not exist: $config -> $dir"
		fi
	done
	echo -e "=======================================================================================\n\n"

	# Remove duplicates
	echo -e "============== Singularity bind mount list for post-annotation analysis ==============="
	postanno_bind_mnts=$(echo "$postanno_bind_mnts" | tr ',' '\n' | awk '!seen[$0]++' | paste -sd ',' -)
	echo -e "$postanno_bind_mnts"
	echo -e "=======================================================================================\n\n\n"
	#-----------------------------------------------------------------------------

	# snakemake --snakefile FinalSnakefile --printshellcmds --dryrun
	snakemake --snakefile FinalSnakefile \
		--cores $threads \
		--printshellcmds \
		--use-singularity \
		--singularity-args "-B $postanno_bind_mnts"
fi
#-----------------------------------------------------------------------------
