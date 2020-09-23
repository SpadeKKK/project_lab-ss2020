# import taxids and original counts
taxids <- readLines(snakemake@input$taxids)
counts.org <- as.numeric(readLines(snakemake@input$counts))

tax <- data.frame(Taxid = taxids, Count = counts.org)
tax <- tax[order(tax$Count, decreasing = TRUE), ]

filters <- snakemake@params$type
for (filter in filters){
    # infer filter type
    if (filter == "none"){
        counts.sel <- 1:nrow(tax)
    } else if (grepl(pattern = "^top\\d+$", x = filter)){
        # select top N biggest taxa, ie. all taxa >= the Nth biggest taxa
        top <- as.integer(sub(pattern = "top", replacement = "", x = filter))
        min <- tax$Count[min(nrow(tax), top)]
        counts.sel <- tax$Count >= min
    } else if (grepl(pattern = "^ratio\\d+(\\.\\d+)?$", x = filter)){
        # select taxa with counts over a certain relative threshould
        ratio <- as.numeric(sub(pattern = "ratio", replacement = "", x = filter))
        counts.sel <- (tax$Count / sum(tax$Count)) >= ratio
    } else if (grepl(pattern = "^min\\d+$", x = filter)){
        # select taxa with a minimum of N hits
        min <- as.integer(sub(pattern = "min", replacement = "", x = filter))
        counts.sel <- tax$Count >= min
    } else if (grepl(pattern = "^max\\d+(\\.\\d+)?$", x = filter)){
        # select taxa with counts over a certain threshould relative to maximal count
        thresh <- as.numeric(sub(pattern = "max", replacement = "", x = filter))
        counts.sel <- tax$Count >= thresh * max(tax$Count)
    } else {
        stop(paste("ERROR: unknown filter type:", filter))
    }
    tax <- tax[counts.sel, ]
}

# export these taxa and associated original counts
writeLines(as.character(tax$Taxid), snakemake@output$taxids)
writeLines(as.character(tax$Count), snakemake@output$counts)
