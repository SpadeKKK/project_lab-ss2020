# Copyright (c) 2016,
# Mathias Kuhring, KuhringM@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.

# additional functions to handle decoys (delete, save, load, merge...)

# patterns to match decoys 
# TODO: should be a user parameter!
DECOY.PATTERN = "^XXX|^RRR|:reversed$"

# index which PSM is as decoy
# provide PSMs separated or if non-unique hits are merged, 
# make sure targets occure before decoys (use psm.move.decoys)!
psm.index.decoys <- function(proteins){
  return(grepl(pattern=DECOY.PATTERN, proteins))
}

# sort proteins hits such that targets occure before decoys
psm.move.decoys <- function(proteins){
  sortProteins <- function(protein){
    splits <- strsplit(x=protein, split=";", fixed=T)[[1]]
    decoys <- grepl(pattern=DECOY.PATTERN, splits)
    splits.new <- c(splits[!decoys], splits[decoys])
    protein.new <- paste(splits.new, collapse=";")
    return(protein.new)
  }
  return( sapply(proteins, sortProteins) )
}

# remove decoy hits from the psms
psm.rem.decoys <- function(psms, 
                          output.file="psms.no-decoys.tsv",
                          decoy.file="psms.decoys.tsv"){
  
  cat("remove decoys\n")
  # remove psms with decoy hits
  is.decoy <- grepl(pattern=DECOY.PATTERN, psms$Protein)
  psms.targets <- psms[!is.decoy, ]
  psms.decoys <- psms[is.decoy, ]
  
  cat("export data\n")
  # export target PSM set
  psm.write(psms.targets, output.file)
  
  # export decoy PSM set
  psm.write(psms.decoys, decoy.file)
  
  cat(paste("retained", nrow(psms.targets), "of", nrow(psms), "matches\n"))
  
  return(psms.targets)
}
