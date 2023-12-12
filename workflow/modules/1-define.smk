import pandas as pd
from pathlib import Path
import os
from snakemake.utils import validate

# load and validate config 
configfile: "inputs/config/config.yaml"
config_schema = "../../workflow/schema/config.schema.yaml"
validate(config, config_schema)


common_config = config["common"]
pipeline_config = config["pipeline"]
deseq_config = config["DESeq2"]

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

validate(meta_pandas, "../../workflow/schema/metadata.schema.yaml")

# Import and validate sample information

sample_id_col = pipeline_config["sample_id"]

SAMPLES = pd.read_table(metadata_file)[sample_id_col].tolist()
print("samples: " + str(SAMPLES))

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
