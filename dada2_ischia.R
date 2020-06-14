####### dada2 with ischia 16S samples following MK English pipeline
# dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html

library(dada2)
library(phyloseq)
path <- "/Users/Winni/Desktop/Mueller_Lab/2019_Ischia/16S_reads" # Change the directory containing the fastq files after unzipping.
list.files(path) 

# Read in the names of the fastq files, and perform string manipulation to get matched lists of the forward and reverse fastq files
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#visualizing the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])

#visualize the quality profile of the reverse reads
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
# Assign the filenames for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,180),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
out
#Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
#113134140 total bases in 461772 reads from 24 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#100649880 total bases in 559166 reads from 27 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#error rates should drop with increasing quality score, and the black and red lines should somewhat agree

#didn't dereplicate because I'm following tutorial v1.12, and the dereplication step is included in 1.8
# Dereplicate
#derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names
#names(derepRs) <- sample.names
#apply the core sample inference algorithm to the filtered and trimmed sequence data.

# Register entry points for exported C++ functions
#methods::setLoadAction(function(ns) {.Call('_dada2_RcppExport_registerCCallable', PACKAGE = 'dada2')})
#denoised
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]
# merge the forward and reverse reads together to obtain the full denoised sequence
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
mergers[[1]]

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# 92 x 2688

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#rows corresponding to samples, columns corresponding to sequence variants. 2688 ASVs, 92 samples
#this gives the length of the merged sequences

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#92 x 2662
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab) # percent of sequence reads that are not chimeras
# 99.78
seqtab2 <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% 245:256]
#remove sequences that are much longer or shorter than expected
dim(seqtab2)


# final step: track the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track
write.table(track, "dada2_readtracking.txt", quote=FALSE, sep="\t")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab2, "/Users/Winni/Desktop/Mueller_Lab/2019_Ischia/16S_reads/silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)

#also add species
taxa <- addSpecies(taxa, "/Users/Winni/Desktop/Mueller_Lab/2019_Ischia/16S_reads/silva_species_assignment_v132.fa.gz")

#save the file just in case
saveRDS(taxa, file="taxa.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
taxa.print
taxa
write.csv(taxa.print, file="ischia_taxa.csv")

#remove chloroplasts and mitochondria
reduce_taxa <- apply(ischia_taxa, 1, function(r) any(r %in% c("Chloroplast", "Mitochondria", "Eukaryota")))
#remove from count table
ischia_asv_reduced <- ischia_asv[,!reduce_taxa]
dim(ischia_asv)
dim(ischia_asv_reduced)
#new count table
write.csv(ischia_asv_reduced, "ischia_asv_nochloro_mito_euk.csv")
#remove from taxonomy table
ischia_taxa_reduced <- ischia_taxa[!reduce_taxa,]
dim(ischia_taxa)
dim(ischia_taxa_reduced)
#new taxonomy table
write.csv(ischia_taxa_reduced, "ischia_taxa_nochloro_mito_euk.csv")
