# Class 6: R Functions
Shazreh Hassan (PID: A13743949)

- [First silly function](#first-silly-function)
- [A second function](#a-second-function)
- [if else Functions](#if-else-functions)

All functions in R have at least 3 things:

- A **name**, we pick this and use it to call our function
- Input **arguments** (can be multiple)
- The **body** lines that do the work

## First silly function

Write a function to add some numbers:

``` r
#in this case, the default for y will be 1 for when user doesn't specify a value for y
add <- function(x, y=1) {
  x+y
}
```

Now we can call this function:

``` r
add(10,5)
```

    [1] 15

``` r
#can add to a vector
add(c(10,10),100)
```

    [1] 110 110

## A second function

Write a a function to generate random nucleotide sequences of a user
specified length:

The `sample()` function can be helpful here; it samples randomly from a
vector.

``` r
sample(c("A","C","G","T"), size=100, replace=TRUE)
```

      [1] "T" "C" "T" "G" "G" "T" "A" "T" "T" "T" "A" "T" "A" "G" "A" "C" "G" "A"
     [19] "G" "G" "C" "G" "G" "A" "C" "G" "G" "G" "T" "T" "G" "G" "G" "G" "T" "G"
     [37] "G" "C" "C" "T" "G" "T" "T" "G" "G" "C" "A" "C" "G" "A" "A" "A" "T" "A"
     [55] "G" "G" "G" "G" "C" "G" "T" "T" "A" "A" "C" "T" "G" "A" "G" "C" "C" "A"
     [73] "G" "G" "C" "A" "G" "T" "G" "T" "A" "A" "T" "T" "C" "T" "A" "A" "T" "A"
     [91] "T" "T" "A" "A" "A" "T" "G" "T" "G" "T"

I want a 1 element long character vector that looks like FASTA format:
“TATTTA” instead of having quotes and spaces in between

``` r
#default for collapse is to have a space
v <- sample(c("A","C","G","T"), size=50, replace=TRUE)
paste(v, collapse = "")
```

    [1] "ATTCTCTTACGGTACACATGGGAATTCCTTTGAGCTGTCCCGGTTTACGA"

Turn this into a function:

``` r
generate_dna <- function(size=50) {
  v <- sample(c("A","T","C","G"),size=size, replace = TRUE)
  paste(v, collapse="")
}
```

Test it:

``` r
generate_dna(60)
```

    [1] "ATGCGTTTCCCTGGAACTGAGGAATTTACATACCATATGCACGAGTAGCTACAGAGGGCA"

## if else Functions

``` r
if(TRUE) {
  cat("HELLO You!")
}
```

    HELLO You!

Add the ability to return a multi-element vector or a single element
fasta-like vector:

``` r
generate_fasta <- function(size=50, fasta=TRUE) {
  v <- sample(c("A","T","C","G"),size=size, replace = TRUE)
  s <- (paste(v, collapse=""))

if(fasta) {
  return(s)
} else {
  return(v)
}
}

generate_fasta(8)
```

    [1] "AAGTGAGG"

Now change from DNA sequence to protein sequence:

``` r
generate_protein <- function(size=50, fasta=TRUE) {
  v <- sample( c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),size=size, replace = TRUE)
  s <- (paste(v, collapse=""))

if(fasta) {
  return(s)
} else {
  return(v)
}
}
  
generate_protein(9) 
```

    [1] "LEYQCPYGP"

Use the new `generate_protein()` function to make random protein
sequences of length 6 to 12 (i.e. one length 6, one length 7, etc. up to
12).

Use a `for()` loop

``` r
# \n is for new line
lengths <- 6:12

for(i in lengths) {
  cat(">",sep="",i,"\n")
  aa <- generate_protein(i)
  cat(aa)
  cat("\n")
}
```

    >6
    EVKEII
    >7
    KDAWCRM
    >8
    TTDGQNGW
    >9
    AESKYSPDN
    >10
    SCPLVCRWDQ
    >11
    KDIPIPVAWIF
    >12
    GPYGNHIPAMNE

A better way to solve this is to use the `apply()` family of functions,
specifically the `sapply` function in this case.

``` r
sapply(lengths, generate_protein)
```

    [1] "VWCYRN"       "ICSGCNE"      "YCPWLNRK"     "DRVAPDRKR"    "SCSFTLWGIL"  
    [6] "NNSDMINPGGA"  "FWSSIWPWCTAP"
