# setwd("V:/MPA-BTB/MPA-PRG/exercise_09")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
#BiocManager::install("GenomicRanges")
#BiocManager::install("Biostrings")
#library("Seqinfo")
#library("Biostrings")

Overlap <- function(seq1, seq2){
  len1 <- length(seq1)
  len2 <- length(seq2)
  overlap_count <- 0
  for (i in 0:len1){
    #print(paste0(seq1[(len1-i):len1], collapse = ""))
    #print(paste0(seq2[1:(i+1)], collapse = ""))
    if (grepl(paste0(seq1[(len1-i):len1], collapse = ""), paste0(seq2[1:(i+1)], collapse = ""))){
      overlap_count <- i + 1
    }
  }
  return(overlap_count)
}

OverlapMat <- function(S){
  len <- length(S)
  SeqMat <- matrix(0,len, len)
  for (i in 1:len){
    seq1 <- strsplit(S[i], "")[[1]]
    for (j in 1:len){
      seq2 <- strsplit(S[j], "")[[1]]
      if (i != j){
        SeqMat[i,j] <- Overlap(seq1,seq2)
      }
    }
  }
  return(SeqMat)
}

GreedySuperstring <- function(S){
  len <- length(S)
  while (len > 1){
    overlapMat <- OverlapMat(S)
    if (max(as.numeric(unlist(overlapMat))) == 0){
      return(S)
    }
    else{
      maxim_indexes = which(m == max(m), arr.ind = TRUE)
      seq1 <- S[maxim_indexes[1]]
      seq2 <- S[maxim_indexes[2]]
      list.remove(S, seq1)
      list.remove(S, seq2)
    }
  }
  return(S)
}

seq1 <- strsplit("TTCA", "")[[1]]
seq2 <- strsplit("CATGC", "")[[1]]
Overlap(seq1, seq2)
S <- c("CATGC","CTAAGT","GCTA","TTCA","ATGCATC")
OverlapMat(S)
GreedySuperstring(S)
