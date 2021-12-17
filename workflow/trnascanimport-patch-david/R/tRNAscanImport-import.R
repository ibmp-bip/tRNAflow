#' @include tRNAscanImport.R
NULL

#' @name import.tRNAscanAsGRanges
#' @aliases import.tRNAscanAsGRanges tRNAscan2GFF tRNAscanID
#' 
#' @title Importing a tRNAscan output file as a GRanges object
#' 
#' @description
#' The function \code{import.tRNAscanAsGRanges} will import a tRNAscan-SE output
#' file and return the information as a GRanges object. The reported 
#' intron sequences are spliced from the result by default, but can also 
#' returned as imported.
#' 
#' The function \code{tRNAScan2GFF} formats the output of 
#' \code{import.tRNAscanAsGRanges} to be GFF3 compliant.
#' 
#' \code{tRNAscanID} generates a unique tRNA ID, which is like the format used 
#' in the SGD annotation 
#' 
#' \code{t*AminoAcidSingleLetter*(*Anticodon*)*ChromosomeIdentifier**optionalNumberIfOnTheSameChromosome*}
#'  
#' Example: tP(UGG)L or tE(UUC)E1.
#'
#' @references 
#' Chan, Patricia P., and Todd M. Lowe. 2016. “GtRNAdb 2.0: An Expanded Database
#' of Transfer Rna Genes Identified in Complete and Draft Genomes.” Nucleic
#' Acids Research 44 (D1): D184–9. doi:10.1093/nar/gkv1309.
#'
#' Lowe, T. M., and S. R. Eddy. 1997. “TRNAscan-Se: A Program for Improved
#' Detection of Transfer Rna Genes in Genomic Sequence.” Nucleic Acids Research
#' 25 (5): 955–64. 
#' 
#' @param input 
#' \itemize{
#' \item \code{import.tRNAscanAsGRanges}: a tRNAscan-SE input file
#' \item \code{tRNAscan2GFF}: a compatible \code{GRanges} object such as the 
#' output of \code{import.tRNAscanAsGRanges}
#' }
#' @param as.GFF3 optional logical for \code{import.tRNAscanAsGRanges}: returns 
#' a gff3 compatible GRanges object directly. (default: \code{as.GFF3 = FALSE})
#' @param trim.intron optional logical for \code{import.tRNAscanAsGRanges}: 
#' remove intron sequences. This changes the tRNA length reported. To retrieve
#' the original length fo the tRNA gene, use the \code{width()} function on the 
#' GRanges object. (default: \code{trim.intron = TRUE})
#'
#' @return a GRanges object
#' @export
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Structstrings DotBracketStringSet
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_trim str_locate_all str_detect
#' @importFrom rtracklayer export.gff3
#'
#' @examples
#' gr <- import.tRNAscanAsGRanges(system.file("extdata", 
#'                                file = "yeast.tRNAscan", 
#'                                package = "tRNAscanImport"))
#' gff <- tRNAscan2GFF(gr)
#' identical(gff,import.tRNAscanAsGRanges(system.file("extdata", 
#'                                file = "yeast.tRNAscan", 
#'                                package = "tRNAscanImport"),
#'                                as.GFF3 = TRUE))
import.tRNAscanAsGRanges <- 
  function(input, as.GFF3 = FALSE, trim.intron = TRUE, min.size = 40){
  # input check
  if(!assertive::is_a_bool(as.GFF3)){
    warning("'as.GFF3' is not a bool. Resetting 'as.GFF3' == FALSE.",
            call. = FALSE)
    as.GFF3 <- FALSE
  }
  if(!assertive::is_a_bool(trim.intron)){
    warning("'trim.intron' is not a bool. Resetting 'trim.intron' == TRUE.",
            call. = FALSE)
    trim.intron <- TRUE
  }  
  if(!assertive::is_a_number(min.size)){
    warning("'min.size' must be numeric. Resetting 'min.size' = 40.",
            call. = FALSE)
    min.size <- 40
  }
  # get tRNAscan as data.frame
  df <- .parse_tRNAscan_improved(input)
  
  # filter tRNA by size, because sometimes tRNAscan predicts 
  # very small tRNA and those lead to bugs in the package.
  df <- df[tRNA_length >= min.size]
  
  # optional: remove intron sequences
  if(trim.intron){
    df <- .cut_introns(df)
  }
  # Contruct GRanges object
  gr <- GRanges(df)
  S4Vectors::mcols(gr)$tRNA_seq <- 
    Biostrings::DNAStringSet(S4Vectors::mcols(gr)$tRNA_seq) 
  S4Vectors::mcols(gr)$tRNA_str <- 
    Structstrings::DotBracketStringSet(S4Vectors::mcols(gr)$tRNA_str)
  S4Vectors::mcols(gr)$tRNA_length <- 
    nchar(as.character(S4Vectors::mcols(gr)$tRNA_seq))
  # sort GRanges object
  gr <- gr[order(GenomeInfoDb::seqnames(gr), BiocGenerics::start(gr))]
  # convert to gff3 compatible GRanges object
  if(as.GFF3){
    gr <- tRNAscan2GFF(gr)
  }
  return(gr)
}

.parse_tRNAscan_improved <- function(file) {
  lines <- readLines(file) 
  breaks <- which(grepl("^[[:space:]]*$", lines))
  size <- breaks - c(0, breaks[1:length(breaks)-1])
  cuts <- as.factor(rep(1:length(breaks), size))
  # get tRNA blocks
  blocks <- split(lines, cuts)
  tail(blocks)
  res <- as.data.table(ldply(blocks, .parse_tRNAscan_block))
  res[, tRNA_length := as.numeric(tRNA_length)]
  # Add intron sequences
  res[!is.na(tRNAscan_intron.start), tRNAscan_intron.seq := substr(tRNA_seq, tRNAscan_intron.start, tRNAscan_intron.end)]
  return(res)
}

.parse_tRNAscan_block <- function(block) {
  
  if(length(block) < 5) {
    warning("A tRNA block size is not correct")
    next
  }
  
  pattern_tRNA_pos_length <- "([a-zA-Z0-9.:^*$@!+_?-|]+).trna([A-Z,-,_,0-9]+) \\(([0-9]+)-([0-9]+)\\).*Length: ([0-9]+) bp"
  pattern_type_anticodon_score <- "Type: ([A-z]{3}).*Anticodon: ([A-z]{3}) at ([0-9]+)-([0-9]+) .*Score: (.*)$"
  pattern_hmm <- "HMM Sc=([-.,0-9]+).*$"
  pattern_struct_score <- "Sec struct Sc=([-.,0-9]+).*$"
  pattern_sequence <- "Seq: ([A-z]+)$"
  pattern_structure <- "Str: ([<,>,.]+)$"
  pattern_possible_intron <- "Possible intron: ([0-9]+)-([0-9]+) \\(([0-9]+)-([0-9]+)\\).*$"

  info <- c(no = NA, chr = NA, start = NA, end = NA, strand = NA, tRNA_length = NA, tRNA_type = NA, 
            tRNA_anticodon = NA, tRNA_anticodon.start = NA, tRNA_anticodon.end = NA, tRNAscan_score = NA, 
            tRNA_seq = NA, tRNA_str = NA, tRNA_CCA.end = NA, tRNAscan_potential.pseudogene = NA, 
            tRNAscan_potential.truncation = NA, tRNAscan_intron.start = NA, tRNAscan_intron.end = NA, 
            tRNAscan_intron.locstart = NA,  tRNAscan_intron.locend = NA, tRNAscan_intron.seq = NA, 
            tRNAscan_hmm.score = NA, tRNAscan_sec.str.score = NA, tRNAscan_infernal = NA)
  
  info_size <- length(info)
  
  for(line in block) {
    if (str_detect(line, "Length:")) {
      tmp <- .regex_custom(line, pattern_tRNA_pos_length)
      info["no"] <- as.numeric(tmp[3])
      info["chr"] <- as.character(tmp[2])
      # Check the strand information
      info["start"] <- ifelse(as.numeric(tmp[4]) < as.numeric(tmp[5]), as.numeric(tmp[4]), as.numeric(tmp[5]))
      info["end"] <- ifelse(as.numeric(tmp[4]) < as.numeric(tmp[5]), as.numeric(tmp[5]), as.numeric(tmp[4]))
      info["strand"] <- ifelse(as.numeric(tmp[4]) < as.numeric(tmp[5]), "+", "-")
      info["tRNA_length"] <- as.numeric(tmp[6])
    } 
    if (str_detect(line, "Type:")) {
      tmp <- .regex_custom(line, pattern_type_anticodon_score)
      info["tRNA_type"] <- as.character(tmp[2])
      info["tRNA_anticodon"] <- as.character(tmp[3])
      info["tRNA_anticodon.start"] <- as.character(tmp[4])
      info["tRNA_anticodon.end"] <- as.character(tmp[5])
      info["tRNAscan_score"] <- as.character(tmp[6])
    }
    if (startsWith(line, "Str")) {
      info["tRNA_str"] <- as.character(.regex_custom(line, pattern_structure)[2])
    } 
    if (startsWith(line, "Seq")) {
      info["tRNA_seq"] = as.character(.regex_custom(line, pattern_sequence)[2])
    }
    if (str_detect(line, "HMM")) {
      info["tRNAscan_hmm.score"] = as.numeric(.regex_custom(line, pattern_hmm)[2])
    } 
    if (str_detect(line, "Sec struc")) {
      info["tRNAscan_sec.str.score"] = as.numeric(.regex_custom(line, pattern_struct_score)[2])
    } 
    if (startsWith(line, "Possible intron:")) {
      tmp = .regex_custom(line, pattern_possible_intron)
      info["tRNAscan_intron.start"] = tmp[2]
      info["tRNAscan_intron.end"] = tmp[3]
      info["tRNAscan_intron.locstart"] = tmp[4]
      info["tRNAscan_intron.locend"] = tmp[5]
    } 
    if (startsWith(line, "Possible pseudogene")) {
      info["tRNAscan_potential.pseudogene"] = TRUE
    }     
    if (startsWith(line, "Possible truncation")) {
      info["tRNAscan_potential.truncation"] = TRUE
    } 
  }
  # last three nucleotides must be CCA and it must be unpaired
  info["tRNA_CCA.end"] = ifelse(endsWith(info["tRNA_seq"], "CCA") & endsWith(info["tRNA_str"], "..."), TRUE, FALSE)
  
  # info should always be same length
  if(length(info) != info_size) warning("WARNING! Problem when parsing data at the following line:", info)
  return(info)
}

# # create data.frame from tRNAscan file
# .read_tRNAscan <- function(file){
#   # parse the information as a list of named lists
#   result <- .parse_tRNAscan(file)
#   # aggregate the data
#   result <- lapply(result, 
#                    function(trna){
#                      res <- list(no = as.numeric(trna$trna[3]),
#                                  chr = as.character(trna$trna[2]))
#                      # If on minus strand
#                      if( as.numeric(trna$trna[5]) < as.numeric(trna$trna[4])){
#                        res <- append(res,
#                                      list(start = as.numeric(trna$trna[5]),
#                                           end = as.numeric(trna$trna[4]),
#                                           strand = "-"))  
#                      } else {
#                        res <- append(res,
#                                      list(start = as.numeric(trna$trna[4]),
#                                           end = as.numeric(trna$trna[5]),
#                                           strand = "+"))  
#                      }
#                      res <- append(res,
#                                    list(tRNA_length = as.numeric(trna$trna[6]),
#                                         tRNA_type = as.character(trna$type[2]),
#                                         tRNA_anticodon = as.character(trna$type[3]),
#                                         tRNA_anticodon.start = as.integer(trna$type[4]),
#                                         tRNA_anticodon.end = as.integer(trna$type[5]),
#                                         tRNAscan_score = as.numeric(trna$type[6]),
#                                         tRNA_seq = as.character(trna$seq[2]),
#                                         tRNA_str = as.character(trna$str[2]),
#                                         tRNA_CCA.end = as.logical(.has_CCA_end(trna$seq[2], 
#                                                                                trna$str[2])),
#                                         # do not force type - optional data
#                                         tRNAscan_potential.pseudogene = 
#                                           ifelse(length(!is.na(trna$pseudogene[2])) != 0,
#                                                  !is.na(trna$pseudogene[2]),
#                                                  FALSE),
#                                         tRNAscan_intron.start = trna$intron[4],
#                                         tRNAscan_intron.end = trna$intron[5],
#                                         tRNAscan_intron.locstart = trna$intron[2],
#                                         tRNAscan_intron.locend = trna$intron[3],
#                                         tRNAscan_hmm.score = trna$hmm[2],
#                                         tRNAscan_sec.str.score = trna$secstruct[2],
#                                         tRNAscan_infernal = trna$infernal[2]))
#                      # if a field returns NULL because it is not set switch to NA, since this
#                      # will persist for data.frame creation
#                      res[vapply(res, is.null, logical(1))] <- NA
#                      return(res)
#                    })
#   # create data.frame
#   df <- lapply(names(result[[1]]), 
#                function(name){ 
#                  unlist(lapply(result, 
#                                function(x){
#                                  unlist(x[[name]])
#                                }))
#                })
#   names(df) <- names(result[[1]])
#   df <- data.frame(df,
#                    stringsAsFactors = FALSE)
#   # set data types
#   rownames(df) <- NULL
#   df$no <- as.integer(df$no)
#   df$tRNA_length <- as.integer(df$tRNA_length)
#   df$tRNAscan_potential.pseudogene <- 
#     as.logical(df$tRNAscan_potential.pseudogene)
#   df$tRNAscan_intron.start <- as.integer(df$tRNAscan_intron.start)
#   df$tRNAscan_intron.end <- as.integer(df$tRNAscan_intron.end)
#   df$tRNAscan_intron.locstart <- as.integer(df$tRNAscan_intron.locstart)
#   df$tRNAscan_intron.locend <- as.integer(df$tRNAscan_intron.locend)
#   df$tRNAscan_hmm.score <- as.numeric(df$tRNAscan_hmm.score)
#   df$tRNAscan_sec.str.score <- as.numeric(df$tRNAscan_sec.str.score)
#   df$tRNAscan_infernal <- as.numeric(df$tRNAscan_infernal)
#   df[is.na(df$tRNA_type),"tRNA_type"] <- "Und"
#   return(df)
# }


# regex wrapper
.regex_custom <- function(line,regex_string){
  unlist(regmatches(line, 
                    regexec(regex_string, 
                            line)))
}

# cuts out introns from sequence and structure
.cut_introns <- function(df){
  #df[!is.na(tRNAscan_intron.start), tRNAscan_intron.seq := substr(tRNA_seq, tRNAscan_intron.start, tRNAscan_intron.end)]
  df[!is.na(tRNAscan_intron.start), tRNA_seq := paste0(substr(tRNA_seq, 1, as.numeric(tRNAscan_intron.start) - 1), substr(tRNA_seq, as.numeric(tRNAscan_intron.end) + 1, tRNA_length))]
  df[!is.na(tRNAscan_intron.start), tRNA_str := paste0(substr(tRNA_str, 1, as.numeric(tRNAscan_intron.start) - 1), substr(tRNA_str, as.numeric(tRNAscan_intron.end) + 1, tRNA_length))]
  df[!is.na(tRNAscan_intron.start), tRNA_length := nchar(tRNA_seq)]
  return(df)
}

#' @rdname import.tRNAscanAsGRanges
#'
#' @export
tRNAscan2GFF <- function(input) {
  .check_trnascan_granges(input, TRNASCAN_FEATURES)
  tRNAscan <- input
  # patch GRanges object with necessary columns for gff3 comptability
  S4Vectors::mcols(tRNAscan)$tRNA_seq <- 
    as.character(S4Vectors::mcols(tRNAscan)$tRNA_seq)
  S4Vectors::mcols(tRNAscan)$tRNA_str <- 
    as.character(S4Vectors::mcols(tRNAscan)$tRNA_str)
  # generate unique tRNA ID
  # SGD like format is used 
  # t*AminoAcidSingleLetter*(*Anticodon*)*ChromosomeIdentifier*
  # *optionalNumberIfOnTheSameChromosome*
  # Example: tP(UGG)L or tE(UUC)E1
  S4Vectors::mcols(tRNAscan)$ID <- tRNAscanID(tRNAscan)
  S4Vectors::mcols(tRNAscan)$type <- "tRNA"
  S4Vectors::mcols(tRNAscan)$type <- 
    as.factor(S4Vectors::mcols(tRNAscan)$type)
  S4Vectors::mcols(tRNAscan)$source <- "tRNAscan-SE"
  S4Vectors::mcols(tRNAscan)$source <- 
    as.factor(S4Vectors::mcols(tRNAscan)$source)
  S4Vectors::mcols(tRNAscan)$score <- NA
  S4Vectors::mcols(tRNAscan)$score <- 
    as.numeric(S4Vectors::mcols(tRNAscan)$score)
  S4Vectors::mcols(tRNAscan)$phase <- NA
  S4Vectors::mcols(tRNAscan)$phase <- 
    as.integer(S4Vectors::mcols(tRNAscan)$phase)
  S4Vectors::mcols(tRNAscan)$score <- 
    as.integer(S4Vectors::mcols(tRNAscan)$phase)
  # arrange columns in correct order
  S4Vectors::mcols(tRNAscan) <- 
    cbind(S4Vectors::mcols(tRNAscan)[,c("source",
                                        "type",
                                        "score",
                                        "phase",
                                        "ID")],
          S4Vectors::mcols(tRNAscan)[,-which(colnames(
            S4Vectors::mcols(tRNAscan)) %in% 
              c("source",
                "type",
                "score",
                "phase",
                "ID"))])
  return(tRNAscan)
}

#' @rdname import.tRNAscanAsGRanges
#'
#' @export
tRNAscanID <- function(input){
  .check_trnascan_granges(input, TRNASCAN_FEATURES)
  tRNAscan <- input
  # create ids based on type, anticodon and chromosome
  chrom <- as.character(GenomeInfoDb::seqnames(tRNAscan))
  chromIndex <- unlist(lapply(seq_along(unique(chrom)), 
                              function(i){
                                rep(i,length(which(chrom == unique(chrom)[i])))
                              }))
  chromLetters <- .get_chrom_letters(length(unique(chromIndex)))
  # Modified version of AA code since the tRNAscan uses a slight deviation
  # from the one defined in Biostrings
  # swap values and names first
  aacode <- Biostrings::AMINO_ACID_CODE
  aacode_names <- names(aacode)
  aacode_values <- aacode
  aacode <- aacode_names
  names(aacode) <- aacode_values
  aacode <- c(aacode, 
              "SeC" = "U", 
              "Und" = "X",
              "fMe" = "M",
              "iMe" = "M",
              "Sup" = "X")
  aa <- unlist(lapply(tRNAscan$tRNA_type, 
                      function(type){
                        aacode[names(aacode) == type]
                      }))
  if( length(aa) != length(tRNAscan$tRNA_anticodon) ||
      length(aa) != length(chromLetters[chromIndex]) ){
    stop("Unknown tRNA type identifier: ",
         paste(tRNAscan$tRNA_type[!(tRNAscan$tRNA_type %in% names(aacode))],
               collapse = ", "),
         "\nKnown type identifier are: ",
         paste(names(aacode), collapse = "','"),
         "'. If this is a genuine identifier, please let us know.",
         call. = FALSE)
  }
  id <- paste0("t",
               aa,
               "(",
               tRNAscan$tRNA_anticodon,
               ")",
               chromLetters[chromIndex])
  # make ids unique if more than one tRNA of the same type is on the same 
  # chromosome
  uniqueID <- id[!duplicated(id)]
  uniqueID <- uniqueID[!(uniqueID %in% id[duplicated(id)])]
  pos <- match(unique(id[duplicated(id)]),id)
  dupID <- unlist(lapply(pos, function(i){
    x <- id[id == id[i]]
    ipos <- which(id == id[i])
    res <- paste0(x,seq(length(x)))
    names(res) <- ipos
    res
  }))
  id[as.numeric(names(dupID))] <- as.character(dupID)
  id
}

# get character values from "A" to "ZZZ" for example
.get_chrom_letters <- function(n){
  let <- list()
  add <- ""
  i <- 1
  while(n > 0){
    # get new letters and retrieve remaining n
    let[[i]] <- .get_chrom_letters2(n,add)
    n <- let[[i]]$remainder
    i <- i + 1
    # 
    batch <- (i-2) %/% length(LETTERS) + 1
    ln <- (i-1) - (batch-1)*length(LETTERS)
    if(n > 0){
      add <- let[[batch]][["let"]][ln]
    }
  }
  let <- unlist(lapply(seq_along(let), function(x){let[[x]][["let"]]}))
  let
}
.get_chrom_letters2 <- function(n, add = ""){
  ret <- list(let = c(),
              remaineder = 0)
  if(n > 0){
    l <- LETTERS[seq_len(n)]
    l <- l[!is.na(l)]
    let <- paste0(add,l)
    ret <- list(let = let,
                remainder = (n - length(let)))
  }
  ret
}
