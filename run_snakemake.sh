set -e

for arg in "$@"; do
  case $arg in
    --sampledir=*)
      sampledir="${arg#*=}"
      shift
      ;;
    *)
      ;;
  esac
done

echo "Sample directory: $sampledir"

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
run_cellranger=`python3 $SCRIPT_DIR/helper_scripts/cache.py RUN_CELLRANGER:` #retrieves starting data from config file
cellranger_reference=`python3 $SCRIPT_DIR/helper_scripts/cache.py CELLRANGER_REFERENCE:` #retrieves cellranger reference genome from config file

mkdir -p $path
python3 $SCRIPT_DIR/helper_scripts/setup.py $project $starting_data $run_cellranger #makes project_name/sample_name/[matrix]or[cellranger] for each sample

python3 $SCRIPT_DIR/helper_scripts/make_symbolic_links.py $project $starting_data $run_cellranger #this will make a bash script
sh $SCRIPT_DIR/helper_scripts/make_symbolic_links.sh #this will make the symbolic links in each matrix directory

threads=`python3 $SCRIPT_DIR/helper_scripts/cache.py THREADS:` #retrieves project name from config file
rpath=`python3 $SCRIPT_DIR/helper_scripts/cache.py RPATH:` #retrieves rpath name from config file

# make .Renviron file to control libpaths
echo "R_LIBS_USER=$rpath" > .Renviron
#-----------------------------------------------------------------------------

# call Snakemake (sans Singularity)
#-----------------------------------------------------------------------------
# snakemake --snakefile Snakefile --printshellcmds --dryrun 
# snakemake --snakefile Snakefile --printshellcmds --dryrun --rerun-triggers mtime
# snakemake --cores $threads --snakefile $SCRIPT_DIR/Snakefile
#-----------------------------------------------------------------------------

# call Snakemake with Singualarity
#-----------------------------------------------------------------------------
snakemake --cores $threads --snakefile $SCRIPT_DIR/Snakefile --printshellcmds --use-singularity --singularity-args "-B $sampledir,$cellranger_reference"
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
	echo "You have the final config file, let the magic begin"
	# snakemake --snakefile FinalSnakefile --printshellcmds --dryrun
	snakemake --cores $threads --snakefile FinalSnakefile --printshellcmds #--forceall 
fi
#-----------------------------------------------------------------------------
