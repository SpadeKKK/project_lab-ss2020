# Copyright (c) 2016,
# Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.

# functions to calculate the FDR and to manage PSM (tsv format) in general

# static values
# TODO: set as default parameters
CONTAMINANT.PATTERN = "^contaminant"
SCORE = "EValue"

# MS-GF+ merges hits to different proteins of the same spectrum with identical 
# peptide sequence and score into one line (with separater ";") and usually only 
# keeps the IDs part of the protein fasta header (splits after the first space).
# However, X!Tandem keeps the complete fasta header which may contain characters
# like ";" in the description. Thus, many characters can not be used.
# Should I use different separators per search engine?
MERGE.SEPARATER = PARAMETERS$merge.separater

psm.parse <- function(psms.org, fdr.cutoff = 0.01){
  # keep only best scoring psms per spectrum
  psms.best <- psm.best(psms.org)
  # merge non-unique psms (same score and peptide sequence)
  psms <- psm.merge(psms.best, export=TRUE)
  # calculate FDR
  fdr.results <- psm.fdr(psms, fdr.cutoff)

  # print some histograms
  # psm.hist(psms.org, fdr.results$score.cutoff, "psm_hist_org.pdf")
  # psm.hist(fdr.results$psms, fdr.results$score.cutoff, "psm_hist_fdr.pdf")

  # make sure psms are separated before decoy/conts removal
  psms.sep <- psm.sep(fdr.results$psms)
  # remove decoys and contaminants
  psms.rem1 <- psm.rem.decoys(psms.sep)
  psms.rem2 <- psm.rem.conts(psms.rem1)

  return(psms.rem2)
}


# read a file with psms in MSGF+ tsv style
psm.read <- function(psms.file){
  cat("read psm file\n")
  # get header
  first.line <- readLines(psms.file, n=1)
  first.line <- sub("#", "", first.line)
  col.names <- strsplit(first.line, "\t", fixed=TRUE)[[1]]
  
  # read psms
  # (disable quotes and comments.char, since they can be part of fasta header or spectra title)
  psms.raw <- read.table(psms.file, sep="\t", skip=1, stringsAsFactors=FALSE, quote = "", comment.char = "")

  colnames(psms.raw) <- col.names
  return(psms.raw)
}

# write psms to a file
psm.write <- function(psms, filename){
  # write header
  writeLines(paste("#", paste(colnames(psms), collapse="\t"), sep=""), filename)
  # write data
  write.table(psms, filename, sep="\t", append=TRUE,
              quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# merge psms with same score and sequence (resp. non-uniques)
# NOTE: molecular weight identical AAs are replaced (I->L, Q->K)
psm.merge <- function(psms, export=FALSE){

  # replace AAs with highly similar molecular weight
  replaceAAs <- function(data){
    data <- gsub(pattern="I", replacement="L", x=data)
    data <- gsub(pattern="Q", replacement="K", x=data)
    return(data)
  }
  
  # iteration over all spectra (title)
  collapseAllMatches <- function(allMatches){
    titles <- unique(allMatches$Title)
    my.list <- vector(mode="list", length=length(titles))
    for (i in 1:length(titles)){
      selection <- allMatches[which(allMatches$Title==titles[i]), ]
      my.list[[i]] <- collapseMatchesPerPeptide(selection)
    }
    return(do.call(rbind, my.list))
  }
  
  # iterate over all peptides per spectra
  collapseMatchesPerPeptide <- function(specMatches){
    peptides <- unique(specMatches$Peptide)
    my.list <- vector(mode="list", length=length(peptides))
    for (i in 1:length(peptides)){
      selection <- specMatches[which(specMatches$Peptide==peptides[i]), ]
      my.list[[i]] <- collapseMatchesPerScore(selection)
    }
    return(do.call(rbind, my.list))
  }
  
  # iterate over all scores per peptide per spectra
  collapseMatchesPerScore <- function(specMatches){
    scores <- unique(specMatches[,SCORE])
    my.list <- vector(mode="list", length=length(scores))
    for (i in 1:length(scores)){
      selection <- specMatches[which(specMatches[,SCORE]==scores[i]), ]
      my.list[[i]] <- mergeSharedPeptides(selection)
    }
    return(do.call(rbind, my.list))
  }
  
  # merge the final selection (same peptide and same score)
  mergeSharedPeptides <- function(specMatches){
    proteins <- paste(specMatches$Protein, collapse=MERGE.SEPARATER)
    final.psm <- specMatches[1, ]
    final.psm$Protein <- proteins
    return(final.psm)
  }
  
  cat("replace AAs (I->L, Q->K)\n")
  psms$Peptide <- replaceAAs(psms$Peptide)
  cat("merge shared peptides\n")
  psms.merged <- collapseAllMatches(psms)
  
  # order by score
  psms.merged <- psms.merged[order(psms.merged[,SCORE]), ]

  if (export){
    cat("export data\n")
    # export merged PSM set
    psm.write(psms.merged, "psms.merged.tsv")
  }

  cat(paste("obtained", nrow(psms.merged), "of", nrow(psms), "matches\n"))
  
  return(psms.merged)
}

# keep psm with best score per spectrum
psm.best <- function(psms, output.file="psms.best.tsv"){
  cat("select best scoring psm per spectrum\n")
  # prepare
  psms.sort <- psms[order(psms$Title, psms[,SCORE]), ]
  
  titles <- unique(psms.sort$Title)
  psms.list <- vector(mode="list", length=length(titles))
  
  # get best matches for each spectra
  for (i in 1:length(titles)){
    sub.res <- psms.sort[which(psms.sort$Title==titles[i]), ]
    min.score <- min(sub.res[,SCORE])
    is.min <- sub.res[,SCORE] == min.score
    psms.list[[i]] <- sub.res[is.min, ]
  }
  
  # combine results
  res.final <- do.call(rbind,  psms.list)
  
  # export
  psm.write(res.final, output.file)

  cat(paste("retained", nrow(res.final), "of", nrow(psms), "matches\n"))
  
  return(res.final)
}

# calculate FDR and remove spectra below threshold
# decoys are not removed anymore!
# deprecated: decoys are recognized by the patterns in DECOY.PATTERN and removed before output
psm.fdr <- function(psms, fdr.cutoff=0.05, output.file="psms.fdr.tsv"){

  cat("prepare data\n")
  # order by score
  psms.res <- psms[order(psms[,SCORE]), ]
  
  # for each PSM, in the Protein list move decoys behind targets
  psms.res$Protein <- psm.move.decoys(psms.res$Protein)
  
  # identify decoy only hits (i.e. first protein is a decoy)
  is.decoy <- psm.index.decoys(psms.res$Protein)
  
  cat("calculate FDR\n")
  # calc FDR per PSM
  fdr <- cumsum(is.decoy) / cumsum(!is.decoy)
  
  cat(paste("use cutoff ", fdr.cutoff, "\n", sep=""))
  # find cut (before first PSM with fdr > cutoff)
  match.count <- nrow(psms.res)
  fdr.cut <- match.count
  bad <- (fdr>fdr.cutoff)
  if (any(bad)) fdr.cut <- which(bad)[1]-1
  psms.filtered <- psms.res[1:fdr.cut, ]
  
  cat("export data\n")
  # export filtered PSM set
  psm.write(psms.filtered, output.file)
  
  cat(paste("retained", nrow(psms.filtered), "of", nrow(psms), "matches\n"))
  
  return(list(psms=psms.filtered, score.cutoff=psms.res[,SCORE][fdr.cut]))
}

# histogram over psm Scores, Targets vs Decoys
psm.hist <- function(psms, cut.off, filename = "fdr_histogram.pdf"){
  require(ggplot2)
  
  psms$Decoy <- psm.index.decoys(psms$Protein)
  
  ggplot(data = psms) + 
    geom_histogram(mapping = aes(x = eval(parse(text=SCORE)), fill = Decoy), bins = 50) + 
    geom_vline(xintercept = cut.off, colour = "red", linetype="dashed")
  
  ggsave(filename)
}

# separate psms with same score and sequence (resp. non-uniques)
psm.sep <- function(psms, output.file="psms.sep.tsv"){
  cat("separate shared peptides\n")
  psm.count <- nrow(psms)
  sep.list <- vector(mode="list", length=psm.count)
  
  for (i in 1:psm.count){
    psm <- psms[i, ]
    proteins <- strsplit(psm$Protein, split=MERGE.SEPARATER, fixed=T)[[1]]
    protein.count <- length(proteins)
    if (protein.count > 1){
      sub.list <- vector(mode="list", length=protein.count)
      for (j in 1:protein.count){
        sub.list[[j]] <- psm
        sub.list[[j]]$Protein <- proteins[j]
      }
      sep.list[[i]] <- do.call(rbind, sub.list)
    }
    else{
      sep.list[[i]] <- psm
    }
  }
  psms.sep <- do.call(rbind, sep.list)
  
  cat("export data\n")
  # export filtered PSM set
  psm.write(psms.sep, output.file)
  
  cat(paste("obtained", nrow(psms.sep), "of", nrow(psms), "matches\n"))
  
  return(psms.sep)
}

# remove contaminat hits from the psms
psm.rem.conts <- function(psms, 
                          output.file.no="psms.no-conts.tsv",
                          output.file.only="psms.only-conts.tsv"){
  
  cat("remove contaminant peptides\n")
  
  is.contaminant <- grepl(pattern=CONTAMINANT.PATTERN, psms$Protein)
  psms.only.conts <- psms[is.contaminant,]
  psms.no.conts <- psms[!is.contaminant,]
  
  cat("export data\n")
  # export filtered PSM set
  psm.write(psms.only.conts, output.file.only)
  psm.write(psms.no.conts, output.file.no)
  
  cat(paste("obtained", nrow(psms.no.conts), "of", nrow(psms), "matches\n"))
  
  return(psms.no.conts)
}