require(ggplot2)

# import taxids
taxids.file <- snakemake@input$taxids
taxids <- as.numeric(readLines(taxids.file))

# import names
taxnames.file <- snakemake@input$names
taxnames.tab <- read.table(taxnames.file, sep = "\t", header = FALSE, stringsAsFactor = FALSE, quote = "")
taxnames <- taxnames.tab$V2
names(taxnames) <- taxnames.tab$V1

# import and normalize associated filtered counts
counts.org.file <- snakemake@input$org
counts.org <- as.numeric(readLines(counts.org.file))
counts.org <- counts.org / sum(counts.org)

# import associated corrected relative counts (weighted or pipasic)
counts.rel.file <- snakemake@input$rel
counts.rel <- as.numeric(readLines(counts.rel.file))

# plot relative abundances (original vs corrected), TOP 100
top100.idx <- order(counts.org, decreasing = TRUE)[1:min(100,length(counts.org))]
counts <- data.frame(Taxid=as.factor(c(taxids[top100.idx],taxids[top100.idx])),
                     Ratio=c(counts.org[top100.idx], counts.rel[top100.idx]),
                     Abundance=c(rep("Original", length(counts.org[top100.idx])),
                              rep("Corrected", length(counts.rel[top100.idx]))))

g <- ggplot(data = counts, mapping = aes(x = Taxid, y = Ratio, fill = Abundance)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(name = "Taxonomic Assignement", labels = taxnames) +
  scale_y_continuous(name = "Relative Abundance Of Assigned Spectra") +
  theme(axis.text.x = element_text(angle = -8, hjust = 0))

ggsave(filename = snakemake@output$plot, width = 25, height = 15, units = "cm")

# plot relative abundances (original vs corrected), TOP 10 original
top10.idx <- order(counts.org, decreasing = TRUE)[1:min(10,length(counts.org))]
counts <- data.frame(Taxid=as.factor(c(taxids[top10.idx],taxids[top10.idx])),
                     Ratio=c(counts.org[top10.idx], counts.rel[top10.idx]),
                     Abundance=c(rep("Original", length(counts.org[top10.idx])),
                              rep("Corrected", length(counts.rel[top10.idx]))))

g <- ggplot(data = counts, mapping = aes(x = Taxid, y = Ratio, fill = Abundance)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(name = "Taxonomic Assignement", labels = taxnames) +
  scale_y_continuous(name = "Relative Abundance Of Assigned Spectra") +
  theme(axis.text.x = element_text(angle = -8, hjust = 0))

ggsave(filename = snakemake@output$plot10org, width = 25, height = 15, units = "cm")

# plot relative abundances (original vs corrected), TOP 10 corrected
top10.idx <- order(counts.rel, decreasing = TRUE)[1:min(10,length(counts.rel))]
counts <- data.frame(Taxid=as.factor(c(taxids[top10.idx],taxids[top10.idx])),
                     Ratio=c(counts.org[top10.idx],counts.rel[top10.idx]),
                     Abundance=c(rep("Original", length(counts.org[top10.idx])),
                              rep("Corrected", length(counts.rel[top10.idx]))))

g <- ggplot(data = counts, mapping = aes(x = Taxid, y = Ratio, fill = Abundance)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(name = "Taxonomic Assignement", labels = taxnames) +
  scale_y_continuous(name = "Relative Abundance Of Assigned Spectra") +
  theme(axis.text.x = element_text(angle = -8, hjust = 0))

ggsave(filename = snakemake@output$plot10rel, width = 25, height = 15, units = "cm")

# select taxids
if (snakemake@params$mode == "meta"){
  taxids.left <- taxids[counts.rel > 0]
} else{
  taxids.left <- taxids[counts.rel > 0 & max(counts.rel)==counts.rel]
}

# export left taxids, i.e. potential taxonomic hits
writeLines(as.character(taxids.left), snakemake@output$taxids)
# export corresponding tax names
write.table(taxnames.tab[taxnames.tab$V1 %in% taxids.left, ], snakemake@output$names, sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
