# MAT 124 Homework 3
# Question 8
# Write a R program for the optimal local sequence alignment
# between the sequences
# A = ACCTAAGG
# B = GGCTCAATCA
# with the following scoring:
# match = +2
# mismatch = -1
# indels = -2
# Coded by: Xinzhou Wang
rm(list = ls())

A = "ACCTAAGG"
B = "GGCTCAATCA"

### Contructing matrix-ready vectors of char
A. <- c(0, unlist(strsplit(A,'')))
B. <- c(0, unlist(strsplit(B,'')))

### Scoring scheme
match = 2
mismatch = -1
indel = -2

### Set up score matrix
s.m = matrix(NA, nrow = length(A.), ncol = length(B.))
# first row and col are 0
s.m[,1] <- c(0)
s.m[1,] <- c(0)

### local alignment
# forces every negative value to be 0 before picking max
for (a in 2:length(A.)) {
  for (b in 2:length(B.)) {
    # diagonal
    if (A.[a] == B.[b]) {
      s.match <- s.m[a-1,b-1] + match
    } else {
      s.match <- s.m[a-1,b-1] + mismatch
    }
    if (s.match < 0) {
      s.match = 0
    }
    # indels on A., vertical
    s.vert <- s.m[a-1,b] + indel
    if (s.vert < 0) {
      s.vert = 0
    }
    
    # indels on B., horizontal
    s.hori <- s.m[a,b-1] + indel
    if (s.hori < 0) {
      s.hori = 0
    }
    
    # pick max out of three paths
    s.m[a,b] = max(c(s.match, s.vert, s.hori))
  }
}

### Traceback parameters
# i = row of highest value
# j = col of highest value
i = which(s.m == max(s.m), arr.ind = TRUE)[1]
j = which(s.m == max(s.m), arr.ind = TRUE)[2]
seq.A <- character()
seq.B <- character()

### Traceback algorithm
while (i > 1 && j > 1) {
  s.local = s.m[i-1,j-1]
  # Logic:
  # Take score from previous [diagonal/vertical/horizontal] cell,
  # Calculate path score,
  # If it matches with actual score, then it is the best path.
  
  # Diagonal
  if (A.[i] == B.[j]) {
    s.local = s.local + match
  } else {
    s.local = s.local + mismatch
  }
  if (s.local == s.m[i,j]) {
    seq.A <- c(A.[i], seq.A)
    seq.B <- c(B.[j], seq.B)
    # If previous diagonal is 0, break loop
    if (s.m[i-1,j-1] == 0) {
      break
    }
    
    i = i - 1
    j = j - 1
    next()
  }
  
  # Vertical
  s.local = s.m[i-1,j]
  if (s.local + indel == s.m[i,j]) {
    seq.A <- c(A.[i], seq.A)
    seq.B <- c('-', seq.B)
    
    if (s.m[i-1,j] == 0) {
      break
    }
    
    i = i - 1
    next()
  }
  
  # Horizontal
  s.local = s.m[i,j-1]
  if(s.local + indel == s.m[i,j]) {
    seq.A <- c('-', seq.A)
    seq.B <- c(B.[j], seq.B)
    
    if (s.m[i,j-1] == 0) {
      break
    }
    
    j = j - 1
    next()
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
max(s.m)
