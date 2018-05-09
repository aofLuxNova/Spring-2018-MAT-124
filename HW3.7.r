# MAT 124 Homework 3
# Question 7
# Write a R program for the optimal global sequence alignment
# between the sequences
# A = ATCGT
# B = TGGTG
# with the following scoring:
# match = +1
# mismatch = -1
# indels = -2
# Coded by: Xinzhou Wang
rm(list = ls())

A = "ATCGT"
B = "TGGTG"

### Contructing matrix-ready vectors of char
# A = 0, A, T, C, G, T
A. <- c(0, unlist(strsplit(A,'')))
B. <- c(0, unlist(strsplit(B,'')))

### Scoring scheme
match = 1
mismatch = -1
indel = -2

### Set up score matrix
s.m = matrix(NA, nrow = length(A.), ncol = length(B.))
# apply indel score*array index over first column. -1 since first cell is 0.
s.m[,1] <- sapply(1:length(A.)-1, function(x) {indel*x})
s.m[1,] <- sapply(1:length(B.)-1, function(x) {indel*x})

#  | b  T  G  G  T  G
# -|------------------
# a| 0 -2 -4 -6 -8 -10
# A|-2
# T|-4
# C|-6
# G|-8
# T|-10

### global alignment
for (a in 2:length(A.)) {
  for (b in 2:length(B.)) {
    # diagonal
    if (A.[a] == B.[b]) {
      s.match <- s.m[a-1,b-1] + match
    } else {
      s.match <- s.m[a-1,b-1] + mismatch
    }
    # indels on A., vertical
    s.vert <- s.m[a-1,b] + indel
    
    # indels on B., horizontal
    s.hori <- s.m[a,b-1] + indel
    
    # pick max out of three paths
    s.m[a,b] = max(c(s.match, s.vert, s.hori))
  }
}

### Traceback parameters
i = length(A.)
j = length(B.)
seq.A <- character()
seq.B <- character()

### Traceback algorithm
while (i > 1 && j > 1) {
  s.global = s.m[i-1,j-1]
  # Logic:
  # Take score from previous [diagonal/vertical/horizontal] cell,
  # Calculate path score,
  # If it matches with actual score, then it is the best path.
  
  # Diagonal
  if (A.[i] == B.[j]) {
    s.global = s.global + match
  } else {
    s.global = s.global + mismatch
  }
  if (s.global == s.m[i,j]) {
    seq.A <- c(A.[i], seq.A)
    seq.B <- c(B.[j], seq.B)
    i = i - 1
    j = j - 1
    next()
  }
  
  # Vertical
  s.global = s.m[i-1,j]
  if (s.global + indel == s.m[i,j]) {
    seq.A <- c(A.[i], seq.A)
    seq.B <- c('-', seq.B)
    i = i - 1
    next()
  }
  
  # Horizontal
  s.global = s.m[i,j-1]
  if(s.global + indel == s.m[i,j]) {
    seq.A <- c('-', seq.A)
    seq.B <- c(B.[j], seq.B)
    j = j - 1
    next()
  }
}

### Alignment sequence correction
check.A <- seq.A
check.B <- seq.B
# length of each sequence without indels after Traceback algorithm
nonindel.A = length(check.A[check.A != '-'])
nonindel.B = length(check.B[check.B != '-'])

# Traceback algorithm contructs the aligned sequences from the end.
# If the alignment has indels, it is possible that the beginning of
# one of the sequences has indels.
# However, this is ignored by the Traceback algorithm, so we have to
# add the missing sequence back.

# If nonindel is shorter than the whole sequence,
# calculate how many letters shorter
# and add back each letter through a while loop.
if (nonindel.A != length(A.)-1) {
  diff = length(A.)-1 - nonindel.A
  
  while (diff > 0) {
    seq.A <- c(A.[diff+1], seq.A)
    seq.B <- c('-', seq.B)
    diff = diff - 1
  }
  
} else if (nonindel.B != length(B.)-1) {
  diff = length(B.)-1 - nonindel.B
  
  while (diff > 0) {
    seq.A <- c('-', seq.A)
    seq.B <- c(B.[diff+1], seq.B)
    diff = diff - 1
  }  
}

### results
# score matrix
s.m
# Aligned A
seq.A
# Aligned B
seq.B
# Final Score
s.m[length(s.m)]
  