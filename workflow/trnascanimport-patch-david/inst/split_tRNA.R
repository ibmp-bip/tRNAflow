#!/usr/bin/env Rscript
message("Loading packages...")
suppressPackageStartupMessages(library(data.table, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(tRNAscanImport, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(tRNA, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(tools, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(plyr, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(stringr, quietly = T, warn.conflicts = F))
suppressPackageStartupMessages(library(argparse, quietly = T, warn.conflicts = F))

# tRNA parts to identify
STRUCTURE = c("anticodonStem", "Dprime5", "DStem", "Dloop", "Dprime3", 
              "acceptorStem", "anticodonLoop", "variableLoop", "TStem", 
              "Tloop", "discriminator")

# Laurence convention for naming things
# The names can be changed here if needed
# but do not change the values of the vector
HEADERS = c(
  "Acc-Stem1" = "acceptorStem.prime5",
  "P8-9"      = "Dprime5",
  "D-stem1"   = "DStem.prime5",
  "D-loop"    = "Dloop",
  "D-stem2"   = "DStem.prime3",
  "p26"       = "Dprime3",
  "Ac-stem1"  = "anticodonStem.prime5",
  "Ac-loop"   = "anticodonLoop",
  "Ac-stem2"  = "anticodonStem.prime3",
  "V-region"  = "variableLoop",
  "T-stem1"   = "TStem.prime5",
  "T-loop"    = "Tloop",
  "T-stem2"   = "TStem.prime3",
  "Acc-stem2" = "acceptorStem.prime3",
  "p73"       = "discriminator"
)

# ArgumentParser options --------------------------------------------------
parser <- ArgumentParser(description='Extract structures information from a tRNAscan structure')
parser$add_argument("infile", help="load a tRNAscan structure file")
parser$add_argument("--summary", action="store_true", help="output an additional summary file")
parser$add_argument("--progress", action="store_true", help="show progress bar while processing tRNA")
args <- parser$parse_args()
# -------------------------------------------------------------------------

infile <- args$infile
summary <- args$summary

if(file.access(infile) == -1)
  stop(sprintf("Specified file ( %s ) does not exist", infile))

#stopifnot(endsWith(infile, TRNASCAN_STRUCT_EXT))

# Function to extract structural information
extract_tRNAs_structure <- function(infile) {
  message("Loading input file: ", infile)
  gr <- try(import.tRNAscanAsGRanges(infile, trim.intron = T, min.size = 40))
  #print(head(gr))
  #dt <- as.data.table(as.data.frame(gr))
  message(paste("Analysing", length(gr), "tRNA(s)"))
  if(summary) {
    message("Generating summary...")
    summary.dt <- as.data.table(gettRNASummary(gr))
    summary.dt <- cbind(summary.dt, dt[, c(.id, seqnames, start, end, width, tRNA_type, tRNA_anticodon, tRNA_seq, tRNA_str)])
    fwrite(summary.dt, paste0(file_path_sans_ext(infile), ".struct_pieces.summary.csv"), row.names = F)
    message("Summary done.")
  }
  message("Identifying tRNAs structures now...")
  tRNA_structure.dt <- as.data.table(
    ldply(1:length(gr), function(i) {
      out <- tryCatch(
        expr = {
          tmp <- as.data.table(as.data.frame(gettRNAstructureSeqs(gr[i], padSequences = F)))
          # Recreate the ID of the tRNA --> Chr1.trn70
          tmp[, id := paste0(gr[i]@seqnames, ".trna", gr[i]$no)]
        },
        error = function(e){ 
          message("\nError in GRanges object at position:", i)
        },
        warning = function(w){
          message("\nWarning catched in GRanges object at position: ", i)
        }
      )
      return(out)
    }, .progress = ifelse(args$progress, "text", "none"))
  )
  message("Done!")
  # Get Laurence convention naming
  setnames(tRNA_structure.dt, HEADERS, names(HEADERS))
  # Put id as the first column
  setcolorder(tRNA_structure.dt, c("id", names(HEADERS)))
  # Save to csv
  fwrite(tRNA_structure.dt, paste0(file_path_sans_ext(infile), ".struct_pieces.csv"), row.names = F)
}

extract_tRNAs_structure(infile)

