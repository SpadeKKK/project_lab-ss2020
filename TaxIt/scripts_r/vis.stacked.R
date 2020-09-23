require(ggplot2)
source(paste0(snakemake@params$wfpath, "/scripts_r/vis.export.R"))


# import original taxids
taxids.org.file <- snakemake@input$org_taxids
taxids.org <- as.numeric(readLines(taxids.org.file))

# import and normalize original counts
counts.org.file <- snakemake@input$org_counts
counts.org <- as.numeric(readLines(counts.org.file))
counts.org <- counts.org / sum(counts.org)
counts.org.label <- "Original"
# create data frame
data.org <- data.frame(Taxid=as.factor(taxids.org),
                       Ratio=counts.org,
                       Abundance=c(rep(counts.org.label, length(counts.org))))

# import filtered/corrected taxids
taxids.file <- snakemake@input$red_taxids
taxids <- as.numeric(readLines(taxids.file))

# import and normalize filtered counts
counts.red.file <- snakemake@input$red_counts
counts.red <- as.numeric(readLines(counts.red.file))
counts.red <- counts.red / sum(counts.red)
counts.red.label <- "Filtered"
# create dataframe
data.red <- data.frame(Taxid=as.factor(taxids),
                       Ratio=counts.red,
                       Abundance=c(rep(counts.red.label, length(counts.red))))

# import corrected relative counts (weighted, pipasic or uniques)
counts.rel.file <- snakemake@input$cor_counts
counts.rel <- as.numeric(readLines(counts.rel.file))
counts.rel.label <- snakemake@params$cor
substr(counts.rel.label, 1, 1) <- toupper(substr(counts.rel.label, 1, 1))
# create dataframe
data.rel <- data.frame(Taxid=as.factor(taxids),
                       Ratio=counts.rel,
                       Abundance=rep(counts.rel.label, length(counts.rel)))

# merge data
counts <- rbind(data.org, data.rel)
filtered <- length(snakemake@params$fil) > 1 || snakemake@params$fil[1] != "none"
if (filtered){
    counts <- rbind(counts, data.red)
}

# define bar order
abundances <- c("Original", "Filtered", "Weighted", "Pipasic", "Uniques")
counts$Abundance <- factor(counts$Abundance, levels = rev(abundances))

# import names
taxnames.file <- snakemake@input$red_names
taxnames.tab <- read.table(taxnames.file, sep = "\t", header = FALSE, stringsAsFactor = FALSE, quote = "")
taxnames <- taxnames.tab$V2
names(taxnames) <- taxnames.tab$V1

# add taxa names and percentage labels
counts$Taxa = taxnames[as.character(counts$Taxid)]
counts$Label = paste0(round(counts$Ratio * 100, digits = 2), "%")
counts$Label[counts$Ratio < 0.05] <- ""

# reduce legend/colors to L top final taxa (with counts >0)
topL <- snakemake@params$topL
if (sum(counts.rel) > 0){
  topL.idx <- order(counts.rel, decreasing = TRUE)[1:min(topL,length(counts.rel))]
  topL.idx <- topL.idx[counts.rel[topL.idx] > 0]
  topL.taxa <- taxnames[as.character(taxids[topL.idx])]
} else if (filtered && sum(counts.red) > 0){
  topL.idx <- order(counts.red, decreasing = TRUE)[1:min(topL,length(counts.red))]
  topL.idx <- topL.idx[counts.red[topL.idx] > 0]
  topL.taxa <- taxnames[as.character(taxids[topL.idx])]
} else {
  topL.idx <- order(counts.org, decreasing = TRUE)[1:min(topL,length(counts.org))]
  topL.idx <- topL.idx[counts.org[topL.idx] > 0]
  topL.taxa <- taxnames[as.character(taxids.org[topL.idx])]
}
counts$Taxa[!counts$Taxa %in% topL.taxa] <- NA
legend.breaks <- sort(unique(counts$Taxa)) #, na.last = TRUE) # uncomment to label NA
legend.labels <- legend.breaks
# legend.labels[is.na(legend.labels)] <- paste("Not in final top", length(topL.idx)) # uncomment to label NA

# number of legend columns
taxa.names <- sort(unique(counts$Taxa)) # , na.last=TRUE) # uncomment to label NA
# taxa.names[is.na(taxa.names)] <- "NA" # uncomment to label NA
legend.ncol <- 4
for (ncol in 4:2){
  taxa.names.size <- c(nchar(taxa.names), rep(0, (-length(taxa.names)) %% ncol))
  if (sum(apply(matrix(taxa.names.size, ncol = ncol),2,max)) <= 100) break
  else legend.ncol <- ncol-1
}

# plot relative abundances (original/filtered vs corrected)
p <- ggplot(data = counts, mapping = aes(x = Abundance, y = Ratio, fill = Taxa, color = Taxa, label = Label)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE), alpha=I(0.5)) +
  coord_flip(ylim = c(0, 1)) +
  scale_y_continuous(expand=c(0.01,0)) +
  scale_fill_discrete(breaks=legend.breaks, labels=legend.labels) +
  scale_color_discrete(breaks=legend.breaks, labels=legend.labels) +
  theme(legend.position="bottom", legend.title=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank()) +
  geom_text(size = 3, position = position_stack(vjust = 0.5, reverse = TRUE), color = "black") +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
  theme(legend.margin=margin(t=-0.3, r=0, b=0, l=0, unit="cm")) +
  guides(fill=guide_legend(ncol=legend.ncol))

margin.height = 0.5 * length(unique(counts$Abundance))
set_panel_size(p = p, file = snakemake@output$plot_pdf,
               margin = unit(1, "mm"),width = unit(18, "cm"), height = unit(margin.height, "cm"))
set_panel_size(p = p, file = snakemake@output$plot_png,
               margin = unit(1, "mm"),width = unit(18, "cm"), height = unit(margin.height, "cm"))