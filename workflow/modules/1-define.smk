import pandas as pd
from pathlib import Path
import os
import sys
from snakemake.utils import validate

# Load and validate config 
configfile: "inputs/config/config.yaml"
config_schema = "../schema/config.schema.yaml"
validate(config, config_schema)


common_config = config["common"]
pipeline_config = config["pipeline"]
deseq_config = config["DESeq2"]


# Check that group_facet is not in intgroup_to_plot
if deseq_config['group_facet'] in deseq_config['intgroup_to_plot']:
    sys.exit(f"Error! 'group_facet' should not be in 'intgroup_to_plot' elements.")

# Set up input directories
main_dir = common_config["projectdir"]
if main_dir is None:
    main_dir = os.getcwd()
main_dir = Path(main_dir)

genome_dir = Path(pipeline_config["genomedir"])
num_threads = pipeline_config["threads"]

input_dir = main_dir / "inputs"
raw_dir = input_dir / "raw" # should already exist
metadata_dir = input_dir / "metadata" # should already exist
metadata_file = metadata_dir / "metadata.txt" # should aready exist

# Validate metadata file
meta_pandas = pd.read_table(metadata_file)

validate(meta_pandas, "../schema/metadata.schema.yaml")

# Import and validate sample information

sample_id_col = pipeline_config["sample_id"]

SAMPLES = pd.read_table(metadata_file)[sample_id_col].tolist()
print("samples: " + str(SAMPLES))


# Check existence of reference files, break if not there
genome_filename = pipeline_config["genome_filename"]
annotation_filename = pipeline_config["annotation_filename"]

genome_filepath = genome_dir / genome_filename
annotation_filepath = genome_dir / annotation_filename

check_ref_fasta = os.path.exists(genome_filepath)

if check_ref_fasta == False:
    sys.exit(f"Error! You are missing the expected reference genome file: {genome_filepath}")

check_ref_genome = os.path.exists(annotation_filepath)

if check_ref_genome == False:
    sys.exit(f"Error! You are missing the expected reference annotation file: {annotation_filepath}")

if common_config["platform"] =="TempO-Seq":
    biospyder_filepath = common_config["biospyder_dbs"] + common_config["biospyder_manifest_file"]

    check_manifest = os.path.exists(biospyder_filepath)

    if check_manifest == False:
        sys.exit(f"Error! You are missing the expected biospyder manifest file: {biospyder_filepath}")


# Set up output directories

output_dir = main_dir / "output"
output_dir.mkdir(parents=True, exist_ok=True)

processed_dir = output_dir / "processed"
processed_dir.mkdir(parents=True, exist_ok=True)

trim_dir = processed_dir / "trimmed"
align_dir = processed_dir / "aligned"
quant_dir = processed_dir / "quantified"

trim_dir.mkdir(parents=True, exist_ok=True)
align_dir.mkdir(parents=True, exist_ok=True)
quant_dir.mkdir(parents=True, exist_ok=True)

qc_dir =  output_dir / "QC"
qc_dir.mkdir(parents=True, exist_ok=True)

log_dir = output_dir / "logs"
log_dir.mkdir(parents=True, exist_ok=True)

print("Using genome: " + str(genome_dir))

# Copy state of workflow for posterity.
shell("cp -rf workflow {log_dir}")
shell("cp -rf {input_dir}/config {log_dir}")
