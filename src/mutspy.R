#!/usr/bin/env Rscript

"
\033[1;4;32m
Mutation counts in bam files.
\033[0m

Usage:
  mutspy.R count [-dsrDP] (-b bam | -l bamlist) [-o out -m mapq] <mut>...
  mutspy.R -h



Options:
  -h                       Show this screen.
  -d                       Calculate mutation distribution along read.
  -s                       Count fw and bw strand separately.
  -r                       Count both bases of overlapping read pairs.
  -n                       Normalize counts to fraction of total in interval.
  -D                       Include duplicate reads (bitflag 0x400) [default: False]
  -P                       Print reads instead of counting [default: False]
  -b --bam=bam             Path to an indexed bam file.
  -l --bamlist=bamlist     Path to file listing indexed bam files one per line.
  -o --out=out             Path to file were output are appended [default: stdout].
  -m --mapq=mapq           Only include data with minimum this mapq [default: 10].
  -i --interval=interval   Size interval to consider [default: 1,1000]

Arguments:
  <mut>                   Mutations is in the form <seqname>:<pos>_<ref>/<alt>, 
                          e.g. chr10:114925316_G/A (SNP), chr10:114925316_GGT/G (deletion),
                          chr10:114925316_G/GAT (insertion), and  chr10:114925316_GGT/GAT 
                          (MNP). Several mutations can be given as 'mut1 mut2 mutn'.
                          Alternatively a file with lines of sitemuts can be given.
                          The mutation N expands into all four bases.
                          Count outputs are 'alt total' or, if stranded,
                          'alt.fw alt.bw total.fw total.bw'
                          Note: mutspy is reference agnostic, only knowns the read sequence.

Requires:
  Rsamtools, docopt.

Version:
  0.3.1

" -> doc

suppressWarnings(suppressMessages(require(docopt)))
suppressWarnings(suppressMessages(require(Rsamtools)))

### FUNCTIONS #####
cigar_parser <- function(c) {
  # Return matrix with one row per operation type (op). Each row consist of:
  # (op type as integer, the length of op, 0-based index in read where op ends,
  # relative 0-based index in reference where op ends).
  # An op-row as (1, 3, 12, 9) indicates 3I in read position 12 and relative reference 9.
  #
  # M, 0, moves q(uery) and r(reference)
  # I, 1, moves q
  # D, 2, moves r
  # S, 3, moves q
  # H, 4, moves neither q not r
  # =,X,N,P  Not Implemented, throw warn
  
  ops <- strsplit(c, "(?<=[MIDNSHPX=])(?=[0-9])", perl = TRUE)[[1]]
  nops <- sum(grepl("M|I|D|S|H", ops))
  res <- matrix(rep(NA, 4*nops), ncol = 4)
  
  q <- 0 #Coordinate on query read
  r <- 0 #Relative coordinate on reference starting from pos
  n <- 1 #op counter
  
  for(op in ops) {
    o <- sub("([[:digit:]]+)([MIDSH]+)", "\\2", op) #op
    l <- as.numeric(sub("([[:digit:]]+)([MIDSH]+)", "\\1", op)) #length of op
    
    if(o == "M") {
      type <- 0
      q <- q + l
      r <- r + l
      
    } else if(o == "I") {
      type <- 1
      q <- q + l
      
    } else if(o == "D") {
      type <- 2
      r <- r + l
      
    } else if(o == "S") {
      type <- 3
      q <- q + l
      
    } else if(o == "H") {
      type <- 4
            
    } else {
      warning(cat("Unsupported operation: ", ops[n], "\n"))
      stop()
    }
    #Update result matrix
    res[n, ] <- c(type, l, q, r)
    n <- n + 1
  }
  return(res)
}

mutation_splitter <- function(mutation) {
  sapply(1:4, function(x)gsub("(.+):([[:digit:]]+)_([GATCN]+)/([GATCN]+)", 
                              paste0("\\", x), mutation))
}

mutation_finder <- function(m, r) {
  #Return boolean of the indeces in reads r that support mutation m.
  #reads r, list output from scanBam(bam,...) (overlapping mutation for speed) with
  #at least pos, seq and cigar elements
  #mutation m, vector of chr, pos, ref, alt, and output from mutation_splitter
  
  #Complex
  if(nchar(m[3]) > 1 & nchar(m[4]) > 1) {
    count <-
      sapply(1:length(r$seq), function(read) {
        ops0 <- cigar_parser(r$cigar[read]) #All operations in read
        dist0 <- as.numeric(m[2]) - r$pos[read] #Genomic distance from read start pos to mutation
        i <- which(ops0[,4] >= dist0)[1] #Index of OP whose genomic coordinate overlaps mutation
        j <- ops0[i, 3] - (ops0[i, 4] - dist0) + 1 #Mutation index position in read
        
        #Is the read sequence from seq start+mutation length (allowing for diequal to mutation sequence
        substr(as.character(r$seq[[read]]), j, j + nchar(m[4]) - 1) == m[4] 
      })
    return(count)
  }

  # SNV
  if(nchar(m[3]) == nchar(m[4])) {
    count <-
      sapply(1:length(r$seq), function(read) {
        ops0 <- cigar_parser(r$cigar[read]) #All operations in read
        dist0 <- as.numeric(m[2]) - r$pos[read] #Genomic distance from read pos to mutation
        i <- which(ops0[,4] >= dist0)[1] #Index of op whose genomic coordinate overlaps mutation
        j <- ops0[i, 3] - (ops0[i, 4] - dist0) + 1 #Mutation index position in read
        
        #Is the read sequence equal to mutation sequence
        substr(as.character(r$seq[[read]]), j, j + nchar(m[3]) - 1) == m[4]
      })
    return(count)
  }
  
  #Simple INSERTION mutation
  if(nchar(m[3]) == 1 & nchar(m[4]) > 1) {
    #Initialize result count
    count <- rep(F, length(r$cigar))
    
    #Filter read for cigars with any insertions
    hasins <- which(grepl(".*I.*", r$cigar))
    
    if(length(hasins) > 0) {
      #If insertions found, parse those to look for the particular m ins
      isins <-
        sapply(hasins, function(x) {
          any(apply(cigar_parser(r$cigar[x]), 1, function(op) {
            op[1] == 1 & #op is an insertion.
              substr(r$seq[x], op[3] - op[2], op[3]) == m[4] & #op read seq match mutation
              op[4] + r$pos[x] - 1 == m[2] #Genomic position match mutation
          }))
        })
      count[ hasins[ isins ] ] <- T #Update count with any m ins
    }
    return(count)
  }
  
  #Simple DELETION mutation
  if(nchar(m[3]) > 1 & nchar(m[4]) == 1) {
    #Initialize result count
    count <- rep(F, length(r$cigar))
    
    #Filter read for cigars with any deletions (faster with large depths)
    hasdel <- which(grepl(".*D.*", r$cigar))
    
    if(length(hasdel) > 0) {
      #If deletions found, parse those to look for the particular m deletion
      isdel <- 
        sapply(hasdel, function(x) {
          any(apply(cigar_parser(r$cigar[x]), 1, function(op) {
            op[1] == 2 & #op is a deletion.
              op[2] == nchar(m[3]) - 1 & #op has same length as mutation
              op[4] - op[2] + r$pos[x] - 1 == m[2] #Genomic position match
          }))
        })
      count[ hasdel[ isdel ] ] <- T #Update count with any m del
    }
    return(count)
  }
}

position_finder <- function(m, r) {
  #Return read position of the start of mutation m (last ref for indels), NA if not found.
  #reads r, list output from scanBam(bam,...) (overlapping mutation for speed) with
  #at least pos, seq and cigar elements
  #mutation m, a vector chr, pos, ref, alt, and output from mutation_splitter
  
  # SNV or simple MNV mutation O ~ N 
  if(nchar(m[3]) == nchar(m[4])) {
    pos <-
      lapply(1:length(r$seq), function(x) {
        #x = 1
        dist0 <- as.numeric(m[2]) - r$pos[x] #Genomic distance from read start to mutation  
        unlist(apply( cigar_parser(r$cigar[x]), 1, function(op) {
          #op = cigar_parser(r$cigar[x])[4,]
          if(op[1] == 0 & 
             dist0 <=  op[4] &
             dist0 > (op[4] - op[2]) &
             substr(r$seq[x], dist0, dist0 + nchar(m[4]) - 1) == m[4]) {
            return(op[3])
          }
        })
        )
      })
  }
  
  #Simple INSERTION mutation
  if(nchar(m[3]) == 1 & nchar(m[4]) > 1) {
    #Initialize result count
    pos <- 
      lapply(1:length(r$seq), function(x) {
        unlist(apply(cigar_parser(r$cigar[x]), 1, function(op) {
          # op = cigar_parser(r$cigar[x])[4,]
          if(op[1] == 1 & #op should be an insertion.
             substr(r$seq[x], op[3] - op[2], op[3]) == m[4] &  #op read seq match mutation
             op[4] + r$pos[x] - 1 == m[2]) {  #Genomic position match mutation's position
            return(op[3])
          }
        }))
      })
  }
  
  #Simple DELETION mutation
  if(nchar(m[3]) > 1 & nchar(m[4]) == 1) {
    pos <- 
      lapply(1:length(r$seq), function(x) {
        unlist(apply(cigar_parser(r$cigar[x]), 1, function(op) {
          # op = cigar_parser(r$cigar[x])[4,]
          if( op[1] == 2 & #op is a deletion.
              op[2] == nchar(m[3]) - 1 & #op has same length as mutation
              op[4] - op[2] + r$pos[x] - 1 == m[2] ) { #Genomic position match
            return(op[3])
          }
        }))
      })
  }
  
  pos[sapply(pos, is.null)] <- NA #Replace NULL with NA for mutation no found
  return(as.numeric(unlist(pos))) #Return read pos
}

### PARSE COMMAND LINE OPTS ####
opt <- docopt(doc)

mapq <- as.numeric(opt$mapq)
if(is.null(opt$bamlist)) bams <- opt$bam
if(!is.null(opt$bamlist)) bams <- readLines(opt$bamlist)

if(opt$D) {
  duplicates <- NA
} else {
  duplicates <- F
}

if(length(opt$mut) > 0) {
  
    if(file.exists(opt$mut)) {
        mutations <- readLines(opt$mut)
    } else {
        mutations <- unlist(strsplit(opt$mut, " "))
    } 
  
    #Let mutation N such as _A/N expand into _A/A, _A/G, _A/T, A/C - i.e. all bases
    #NOTE that only one ":" is allowed (and before pos) in order not to expand N in e.g. chrN
    f1 <- function(str)unique(as.vector(sapply(c("A", "T", "C", "G"), function(x) sub("N(?!.*:)", x, str, perl = T))))
    while(any(grepl("N(?!.*:)", mutations, perl = T))) mutations <- f1(mutations)
    mutations <- unique(grep("^$", mutations, invert = T, value = T)) #Remove duplicates and ""s
   
    stopifnot(all(grepl("^.+:[[:digit:]]+_[AGCT]+/[AGCT]+$", mutations)))
}

output_file <- ifelse(opt$out == "stdout", "", opt$out)

### MAIN ####

flag <- scanBamFlag(isDuplicate = duplicates)

if(opt$count) {
  #Split into list of (seq, pos, ref, mut) per provided mutation in opt$mut
  for(bam in bams) {
    for(mutation in mutations) {
      m0 <- mutation_splitter(mutation)
      which <- GRanges(seqnames = m0[1], ranges = m0[2])
      param <- ScanBamParam(what = c("pos", "seq", "cigar", "strand", "qname"),
                            which = which,
                            mapqFilter = mapq,
                            flag = flag)
      reads <- Rsamtools::scanBam(bam, param = param)[[1]]
      
      #Remove reads that does not span complete mutation (e.g. a long insertion)
      i <- reads$pos + nchar(reads$seq) >= as.numeric(m0[2]) + nchar(m0[4])
      reads <- sapply(reads, "[", i, simplify = F)
      
      
      if(!opt$r) { 
        #if not counting overlapping read pairs
        #Remove arbitrarily 1 read among those with duplicated read names
        set.seed(1)
        
        dups <- reads$qname[duplicated(reads$qname)]
        if(length(dups) > 0){
          to_remove <- sort(sapply(dups, function(x) sample(which(reads$qname %in% x), 1)))
          reads <- lapply(reads, function(x)x[ c(1:length(x))[-to_remove] ])
        }
      }               
      
      if(length(reads$seq) == 0) {
        if(!opt$s) count <- "0 0"
        if(opt$s) count <- "0 0 0 0"
        
      } else {
        i <- mutation_finder(m = m0, r = reads)
        
        #Do a not stranded, or stranded (fw and bw wise)
        if(!opt$s) {
          count <- paste(sum(i), length(i), collapse = " ")
        }
        if(opt$s) {
          count <- 
            paste(paste(table(reads$strand[i])[c("+", "-")], collapse = " "), 
                  paste(table(reads$strand)[c("+", "-")], collapse = " "),
                  collapse = " ")
        }
      }
      
      #Write results of count
      if(opt$P) {
        #Print all reads
        suppressWarnings(
          write.table(cbind(as.data.frame(sapply(reads, "[", i)), 
                            "bam" = rep(bam, sum(i)), "mutation" = rep(mutation, sum(i))), 
                      file = output_file, sep = " ", row.names = F, quote = F, append = T)
        )
      } else {
        #Write counts
        write(paste(basename(bam), paste(m0, collapse = " "), count, collapse = " "),
              file = output_file,
              append = T)
      }
    }
  }
}

