#Reading tables
tandem <- read.table(file = 'tandem.output.tsv', sep = '\t', header = FALSE)
knn <- read.table(file = 'final_knn_result_k50.csv', sep = ',', header = FALSE)

#Fixing headers
names(tandem)[names(tandem) == "V1"] <- "Raw"
names(knn)[names(knn) == "V1"] <- "Raw"
names(tandem)[names(tandem) == "V2"] <- "Scan"
names(knn)[names(knn) == "V2"] <- "Scan"
names(tandem)[names(tandem) == "V3"] <- "Peptide"
names(knn)[names(knn) == "V3"] <- "Peptide"
names(tandem)[names(tandem) == "V4"] <- "Protein"
names(knn)[names(knn) == "V4"] <- "Protein"
names(tandem)[names(tandem) == "V5"] <- "Hyperscore"
names(knn)[names(knn) == "V5"] <- "Hyperscore"
names(tandem)[names(tandem) == "V6"] <- "EValue"
names(knn)[names(knn) == "V6"] <- "EValue"

#Hyperscore Values
hytdm <- tandem$V5
hyknn <- knn$V5


#Histo Plots
hist(tandem$V5,
     main="Hyperscore of XTandem Result",
     col="darkmagenta",
     freq=FALSE
)

hist(log(knn$V5+1),
     main="Hyperscore of KNN Result",
     col="chocolate",
     freq=FALSE
)

#Creating Table With excluding unwanted columns
#AT: Tandem Table
#Bk: KNN Table
AT <- data.frame(raw = tandem$Raw , scan = tandem$Scan , hyperscore = tandem$Hyperscore)
BK <- data.frame(raw = knn$Raw , scan = knn$Scan , hyperscore = knn$Hyperscore)
AT$raw = "Run1_U4_2000ng.0free.mgf"
BK$raw = "Run1_U4_2000ng.0free.mgf"

#Groupby each table to represent the highest score of each scan number
library(dplyr)
AT <- AT %>% group_by(scan)%>%slice(1)%>%ungroup()
BK <- BK %>% group_by(scan)%>%top_n(n = 1, wt = hyperscore)%>%ungroup()

#Create new table 'new' in which both AT and BK are joined togther based on the common scan numbers
new <- left_join(AT,BK,by=c("scan","raw"))
new <- new %>% filter(!is.na(hyperscore.y))


##HexBin Plots

TandemScore <- new$hyperscore.x
KnnScore1 <- new$hyperscore.y
KnnScore <- log(new$hyperscore.y+1)

cols <- colorRampPalette(c("darkorchid4","darkblue","green","yellow", "red") )

par(mar=c(5,7,4,2))
plot(hexbin(TandemScore,KnnScore1), 
    
     colorcut = seq(0,1,length.out=24),
     colramp = function(n) cols(24), 
     xlab = "TandemScore", 
     ylab = "KNNScore \n",
     main="Hyperscore Correlation of 'X!Tandem search' and 'KNN search' results") 




