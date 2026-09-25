import pandas as pd
from pathlib import Path
import os
from snakemake.utils import validate

# load and validate config 
configfile: "inputs/config/config.yaml"
config_schema = "../schema/config.schema.yaml"
validate(config, config_schema)


common_config = config["common"]
pipeline_config = config["pipeline"]
deseq_config = config["DESeq2"]

# Check that deseq_facet is not in intgroup_to_plot

if deseq_config['deseq_facet'] is not None and deseq_config['deseq_facet'] in deseq_config['intgroup_to_plot']:
    sys.exit(f"Error! The column used for 'deseq_facet' should not be in 'intgroup_to_plot' elements.")

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


# Import and validate metadata file
metadata_file = metadata_dir / common_config["metadata_file"] # should aready exist
meta_pandas = pd.read_table(metadata_file)

validate(meta_pandas, "../schema/metadata.schema.yaml")

# Import and validate sample information

sample_id_col = pipeline_config["sample_id"]

SAMPLES = pd.read_table(metadata_file)[sample_id_col].tolist()
print("samples: " + str(SAMPLES))

if common_config["platform"] == "DRUG-Seq":
    # Special information onlu
    LIBRARIES = meta_pandas["library_ID"].unique().tolist()
    LIBRARY_SAMPLES = {lib: meta_pandas[meta_pandas["library_ID"] == lib]["sample_ID"].tolist() 
                   for lib in LIBRARIES}
    
    # Create a dictionary mapping (library, sample_id) to barcode
    SAMPLE_BARCODE_MAP = {}
    for _, row in meta_pandas.iterrows():
        key = (row["library_ID"], row["sample_ID"])
        SAMPLE_BARCODE_MAP[key] = row["sample_barcode"]

    def get_barcode_for_sample(wildcards):
        """Get the barcode for a given sample_id in a library"""
        key = (wildcards.library, wildcards.sample)
        return SAMPLE_BARCODE_MAP.get(key, "")
    
    print(f"Found {len(LIBRARIES)} libraries: {LIBRARIES}")

    DEDUP_METHODS = pipeline_config["umi_dedup_method"].split(sep = " ")
    DEDUP_METHODS_FOR_COMPARISON = [m for m in DEDUP_METHODS if m != "NoDedup"] # To avoid stats on NoDedup vs NoDedup in rule quantify_dedup

    dedup_method_for_downstream = pipeline_config["umi_method_definitive"]
    
    print(f"UMI deduplication method(s) used by STARsolo: {DEDUP_METHODS}")

# Import and validate contrasts file
contrasts_dir = input_dir / "contrasts"
contrasts_file =  contrasts_dir / common_config["contrasts_file"]

contrasts_pandas = pd.read_table(contrasts_file, sep = "\t")
# Check if contrasts file is tab-delmited with two columns
if contrasts_pandas.shape[1] != 2:
    sys.exit(f"Error! The contrasts file should have two tab-delimited columns, but it has {contrasts_pandas.shape[1]} columns. \n Double check that your file is tab-separated, not space separated.")

# Check existence of reference files, break if not there
if common_config["platform"] == "DRUG-Seq":
    genome_filepath = os.path.join(genome_dir, pipeline_config["genome_filename"])
    annotation_filepath = os.path.join(genome_dir, pipeline_config["annotation_filename"])
    ercc_fasta = os.path.join(pipeline_config["erccdir"], pipeline_config["ercc_fasta"])
    ercc_gtf = os.path.join(pipeline_config["erccdir"], pipeline_config["ercc_gtf"])

    if pipeline_config["include_ercc"] == True:
        # Can't check for combined files yet, they'll be created by a rule
        genome_fasta = os.path.join(pipeline_config["genomedir"], f"{pipeline_config['genome_name']}_ERCC_combined.fa")
        genome_gtf = os.path.join(pipeline_config["genomedir"], f"{pipeline_config['genome_name']}_ERCC_combined.gtf")
        index_dir = os.path.join(pipeline_config["genomedir"], pipeline_config["genome_name"] + "_ERCC_" + "STARindex")
        index_genome = os.path.join(pipeline_config["genomedir"], pipeline_config["genome_name"] + "_ERCC_" + "STARindex", "Genome")
        
        if os.path.exists(ercc_fasta) == False:
            sys.exit(f"Error! You are missing the expected ERCC fasta file: {ercc_fasta}")
        
        if os.path.exists(ercc_gtf) == False:
            sys.exit(f"Error! You are missing the expected ERCC fasta file: {ercc_gtf}")
    else:
        genome_fasta = genome_filepath
        genome_gtf = annotation_filepath
        index_dir = os.path.join(pipeline_config["genomedir"], pipeline_config["genome_name"] + "_STARindex")
        index_genome = os.path.join(pipeline_config["genomedir"], pipeline_config["genome_name"] + "_STARindex", "Genome")
else:
    # TempO-Seq or RNA-seq do not use ERCC spike-ins
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
    if common_config["biospyder_dbs"] is None or common_config["biospyder_manifest_file"] is None:
        sys.exit(f"Error! You have biospyder information set to null in the config, but it is required for TempO-Seq analysis.")
    
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

sm_temp_dir = output_dir/ ".rodaf_internal"
sm_temp_dir.mkdir(parents=True, exist_ok=True)

print("Using genome: " + str(genome_dir))

# Copy state of workflow for posterity.
shell("cp -rf workflow {log_dir}")
shell("cp -rf {input_dir}/config {log_dir}")
