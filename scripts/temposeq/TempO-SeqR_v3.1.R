# TempO-SeqR Script Provided by BioSpyder
# Project-specific customization by Matthew J. Meier
# Command Line Arguments:
# 1: FASTA Reference File
# 2: Directory of FASTQ files to align
# 3: Number of CPUs to use

# Activate Conda environment for TempO-Seq and Required Software
# system2("conda activate temposeq")

run.alignment.analysis.fn <- function(fasta, queryDirectory, num_clusters){
  
  curr_dir <- getwd()
  logfile <- paste0(curr_dir, "/TempOSeq_log.txt")
  sink(logfile)
  
  # Import the QuasR library
  library(QuasR) 
  
  # Get full path for FASTA file (in case partial path provided to command line)
  fasta <- system2('readlink', paste0('-f ', fasta), stdout=T)
  print(paste0("FASTA: ", fasta))
  
  # Get path of reference directory (trim filename)
  ref_dir <- dirname(fasta)
  print(paste0("ref_dir: ", ref_dir))
  print(paste0("FASTQ_dir: ", queryDirectory))
  print(paste0("Removing old temp files..."))
  
  ###############################################################################
  # Setup folders:
  # Remove temporary files from previous alignments if they are present
  if(dir.exists("_STARtmp")){
    unlink("_STARtmp", recursive=TRUE)
  }
  print(paste0("Done."))
  print(paste0("Creating new temp directories..."))
  # Create temporary directories for temp-alignment files
  dir.create("_STARtmp")
  print(paste0("Changing directory to STARtmp..."))
  setwd("_STARtmp")
  print(paste0("Creating new temp directories"))
  dir.create("proba")
  dir.create("fastqs")
  print(paste0("Copying FASTQs"))
  
  fastq_paths <- list.files(path=queryDirectory,
                            recursive = F,
                            full.names = T,
                            all.files = F,
                            pattern = "fastq|fq",
                            ignore.case = T)
  write.table(fastq_paths,
              file = "fastq_list.txt",
              quote = F,
              row.names = F,
              col.names = F)
  #copy_fastqs <- paste0("ln -s ", queryDirectory, "* ./fastqs/")
  #copy_fastqs <- paste0("ln -s ", paste(fastq_paths, collapse=' '), " ./fastqs/")
  copy_fastqs <- paste0("cat fastq_list.txt | xargs -I % ln -s % ./fastqs/")
  
  print(head(fastq_paths))
  print(copy_fastqs)
  system(copy_fastqs)
  queryDirectory <- "fastqs/"
  
  # Get path of GTF annotation file
  annotFile <- gsub("\\.fa", ".gtf", fasta)
  print(paste0("Making cluster"))
  # Make the number of clusters you want to use
  cl <- makeCluster(num_clusters)
  # Number of CPUs per STAR instance
  cores <- floor(num_clusters/7)
  print(paste0("Getting list of files"))
  
  # Get list of FASTQ files
  # This globbing is not going to be foolproof. Very permissive.
  get_files <- paste0("ls ", queryDirectory, "> samplefile.txt")
  append_files <- paste0("sed -e 's|^|", queryDirectory, "|' -i samplefile.txt")
  system(get_files)
  print(get_files)
  print(append_files)
  
  # Add FASTQ path to file name
  system(append_files)
  print(paste0("Loading alignment function..."))
  
  ###############################################################################
  # Alignment Function
  aligning.fn <- function(x){
    # Test the FASTQ file: is it zipped or not?
    if(substr(x, (nchar(x)-2),nchar(x))==".gz"){
      # If the FASTQ is gzipped
      cmd.gz <- paste("STAR --genomeDir ",
                      ref_dir,
                      " --genomeLoad LoadAndKeep ",
                      " --limitBAMsortRAM 50000000000 ",
                      " --readFilesIn ", x,
                      " --readFilesCommand zcat ",
                      " --runThreadN ", cores,
                      " --outSAMtype BAM SortedByCoordinate ",
                      " --scoreDelOpen -10000 ",
                      " --scoreInsOpen -10000 ",
                      " --outFilterMultimapNmax 1 ",
                      " --outFilterMismatchNmax 2 ",
                      " --outSAMunmapped Within ",
                      " --outFileNamePrefix ", x,
                      sep="")
      print(cmd.gz)
      system(cmd.gz)
      
    }else{
      # If the FASTQ is not zipped
      cmd.notZipped <- paste("STAR --genomeDir ",
                             ref_dir,
                             " --genomeLoad LoadAndKeep ",
                             " --limitBAMsortRAM 50000000000 ",
                             " --readFilesIn ", x,
                             " --runThreadN ", cores,
                             " --outSAMtype BAM SortedByCoordinate ",
                             " --scoreDelOpen -10000 ",
                             " --scoreInsOpen -10000 ",
                             " --outFilterMultimapNmax 1 ",
                             " --outFilterMismatchNmax 2 ",
                             " --outSAMunmapped Within ",
                             " --outFileNamePrefix ", x,
                             sep="")
      print(cmd.notZipped)
      system(cmd.notZipped)
    }
    cmd.samtools <- paste("samtools index ",
                          paste(x,"Aligned.sortedByCoord.out.bam",sep=""),
                          sep="")
    system(cmd.samtools)
  }
  # End Alignment
  print(paste0("Alignment function loaded."))
  
  ###############################################################################
  # Import sample file
  print(paste0("Importing sample file names..."))
  sampleF <- read.table("samplefile.txt", header=F, as.is=T)
  sampleF <- sampleF$V1
  # Get sample names
  sampleFile2.bam <- as.data.frame(sampleF)
  
  # Run star Alignment
  print(paste0("Running STAR aligner..."))
  parallel::mclapply(sampleF, FUN=aligning.fn, mc.cores=7)
  
  sampleFile2.bam[,1] <- paste0(sampleFile2.bam[,1],
                                "Aligned.sortedByCoord.out.bam")
  sampleFile2.bam[,2] <- gsub("fastqs/||\\.fastq.*", "", sampleFile2.bam[,1])
  colnames(sampleFile2.bam) <- c("FileName", "SampleName")
  write.table(sampleFile2.bam, "sample.bam.txt", quote=F, row.names=F, sep="\t")
  sampleFile2.bam.file <- "sample.bam.txt"
  
  # Create proj file:
  # This will make the BAM files, uses BWA internally
  print(paste0("Running qAlign function..."))
  proj2 <- qAlign(sampleFile2.bam.file,
                  fasta,
                  paired="no",
                  cacheDir="proba",
                  clObj=cl)
  
  print(paste0("Loading rtracklayer and GenomicFeatures..."))
  library(rtracklayer)
  library(GenomicFeatures)
  
  txStart <- import.gff(annotFile, format="gtf")
  names(txStart) <- txStart@seqnames
  
  cnt <- qCount(proj2, txStart, clObj = cl)
  
  # Make dataframe as count table
  setwd("..")
  cnmat           <- as.data.frame(cnt, header=TRUE)
  colnames(cnmat) <- gsub("fastqs/", "", colnames(cnmat))
  cnmat$width     <- NULL
  write.csv(cnmat, "count_table.csv")
  cat("Counts Table Completed\n")
  
  # Get Mapped/Unmapped Counts and Output
  a_frame     <- data.frame(alignmentStats(proj2), header=TRUE, check.names=FALSE) 
  unmapped    <- a_frame$unmapped 
  mapped      <- colSums(cnmat) 
  two         <- rbind(mapped) 
  three       <- data.frame(two, check.names=FALSE)     
  three$width <- NULL 
  ll          <- rbind(unmapped=unmapped, mapped=three)
  
  write.csv(as.matrix(ll), "mapped_unmapped.csv")
  cat("Mapped/Unmapped Table Completed\n")
  
  sessioninfo::session_info()
  
  #stop logging
  sink()
  
} #End alignment/analysis function

###############################################################################
# Extract Arguments from command line
args = (commandArgs(TRUE))
fasta          <- args[1]
queryDirectory <- args[2]
num_clusters   <- as.numeric(args[3])

###############################################################################
## Run program 
tryCatch({suppressWarnings(run.alignment.analysis.fn(fasta, queryDirectory, num_clusters))},
         error=function(err){
           sink()
         })
