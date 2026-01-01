# setwd("V:/MPA-BTB/MPA-PRG/exercise_09")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
#BiocManager::install("GenomicRanges")
#BiocManager::install("Biostrings")
#library(Seqinfo)
library(Biostrings)

Overlap <- function(seq1, seq2){
  len1 <- length(seq1)
  len2 <- length(seq2)
  max_overlap <- min(len1, len2)
  
  for (k in max_overlap:1){
    if (identical(seq1[(len1-k+1):len1], seq2[1:k])){
      return(k)
    }
  }
  return(0)
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
  while (length(S) > 1){
    overlapMat <- OverlapMat(S)
    max_ov <- max(overlapMat)
    
    if (max_ov == 0){
      return(paste(S, collapse = ""))
    }
    
    idx <- which(overlapMat == max_ov, arr.ind = TRUE)[1,]
    i <- idx[1]
    j <- idx[2]
    
    s1 <- S[i]
    s2 <- S[j]
    
    ov <- max_ov
    merged <- paste0(s1, substr(s2, ov + 1, nchar(s2)))
    
    S <- S[-c(i, j)]
    S <- c(S, merged)
  }
  
  return(S)
}


seq1 <- strsplit("TTCA", "")[[1]]
seq2 <- strsplit("CATGC", "")[[1]]
Overlap(seq1, seq2)
S <- c("CATGC","CTAAGT","GCTA","TTCA","ATGCATC")
OverlapMat(S)
GreedySuperstring(S)
