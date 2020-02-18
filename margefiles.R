WORKDIR <- 'data/'
args <- commandArgs(trailingOnly = TRUE)
WORKDIR <- args[1]


if(!dir.exists(WORKDIR)){
	stop("Please provide a dirctory with bam files.", call.=FALSE)
}

cat('Setting WORKDIR to:', WORKDIR, '\n' )
setwd(WORKDIR)
# Load library
library(data.table)

# Get a List of all files in directory named with a key word, say all `.csv` files
  filenames <- list.files("./", pattern="*.csv", full.names=TRUE)

# read and row bind all data sets
data <- rbindlist(lapply(filenames,fread),fill=TRUE)

# Save in a single file
write.csv(data, file="marged.csv", row.names=T)

