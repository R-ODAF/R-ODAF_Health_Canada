import pandas as pd
import os

indir = "/mnt/nextseqData/Runs"
indirs = ["220331_NB501400_0164_AHTWJ2BGXK_FASTQ/Dechlorane_Plus", "220406_NB501400_0165_AHTWC7BGXK_FASTQ", "220407_NB501400_0166_AHTYG5BGXK_FASTQ", "220408_NB501400_0167_AHTWHGBGXK_FASTQ", "220411_NB501400_0168_AHTW7LBGXK_FASTQ"] 

outdir = "data/raw"

metadata = pd.read_csv("data/metadata/metadata.csv")

for sample in metadata["Sample_ID"]:
    outfile = outdir + "/" + sample + ".fastq.gz"
    print(outfile)
    #os.system("rm {}".format(outfile))
    for i in indirs:
        file = indir + "/" + i + "/" + sample + "_S*_R1_001.fastq.gz"
        print(file)
        #os.system("cat {} >> {}".format(file, outfile))
