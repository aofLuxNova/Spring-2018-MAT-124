# MAT 124 Homework 3
# QUestion 6
# Write a R program that searches the amino acid sequence
# HGVYYDPSKDLIAEIQ
# in the HIV complete protein sequence found in the webpage
# https://www.ncbi.nlm.nih.gov/nuccore/4558520
# with protein_id = AAC82598.2
# Coded by: Xinzhou Wang
library(seqinr)
path = "E:\\Documents\\Programming\\R\\Spring-2018-MAT-124\\HIV1.fasta"
seq = "HGVYYDPSKDLIAEIQ"

checkAASeq = function(path, seq) {
  # read fasta file and specify type of sequence as amino acid sequence
  # as.string = TRUE for outputting sequence as string instead of list of vectors of chars
  # seqonly = TRUE for returning only the sequence with no annotations
  d <- seqinr::read.fasta(path, seqtype = c('AA'), as.string = TRUE, seqonly = TRUE)
  
  # Is seq in d
  l <- grepl(seq, d) # TRUE
  # How many seq are in d
  n <- grep(seq, d) # 1
  
  return(list(l,n))
}

checkAASeq(path,seq)
