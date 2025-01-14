# Make Figure 1: PCA and mapping
rm(list = ls())
library(leaflet)
#library(viridis)
library(colourvalues) # https://symbolixau.github.io/colourvalues/
library(extrafont)

# ------------------------------------------------------------------------------------- #
# Upload data
# ------------------------------------------------------------------------------------- #

# upload coverage information
setwd("/Users/rozenn/GaTech Dropbox/CoS/BioSci/BioSci-Ratcliff/Rozenn/4.PandoProject/11-field-lab-work/7-variants/manual_filtering/data/GL/forF1/")
dp <- read.table("cov_pando_friends_samples.txt", sep = "\t", header = T)
cov <- dp[,-1]
dim(cov) # 22888 x 413

indMean <- colMeans(cov) # Obtain mean per sample 
plot(sort(indMean), xlab="Sample", ylab="Mean coverage per SNP per individual", pch = 16)
abline(h=4, col = "red")

# mean coverage per SNP and % removed
snpMean <- rowMeans(cov) # Obtain mean per sample 
length( which(snpMean < 4 ) ) / dim(cov)[1] *100 # 54% of SNPs removed due to low coverage
(dim(cov)[2]-213) / dim(cov)[2] *100 # 48% of the individuals removed because of low coverage

# Keep individuals with more than 4 reads total
keepInd <- indMean > 4
length(which(keepInd=="TRUE"))

# sample info
ids <- read.table("combined_IDs.txt",  sep="\t", header=T)

# individuals to keep: coverage > 4 and initial dataset
indices <- which(keepInd=="TRUE" & ids$Group=="pando-friends")
length(indices) # 184 individuals

mean(indMean[indices])  # mean = 8.1

#  upload GLs
gl <- read.table("pntest_combined_variant_40miss_14diff_half_mean.txt", sep = " ") 
dim(gl) # 22888 x 413
gl_2 <- gl[,indices]
dim(gl_2) #  22888 x 184
ids_2 <- ids[indices,]
dim(ids_2) #  184 x 6

# SUMMARY: initially 413 individuals, we remove 117 individuals (newer dataset) and 112 individuals (too low coverage)

# ------------------------------------------------------------------------------------- #
# PCA
# ------------------------------------------------------------------------------------- #
res.pca <- prcomp(t(gl_2),center = TRUE) # PCA on point estimate for Pando and friends
sum_pca <- summary(res.pca)
ax1 <- round(sum_pca$importance[2]*100) # proportion of variance explained by first axis
ax2 <- round(sum_pca$importance[5]*100) # proportion of variance explained by second axis

pc1 <- res.pca$x[,1]
pc2 <- res.pca$x[,2]

# ------------------------------------------------------------------------------------- #
# frame dataset
# ------------------------------------------------------------------------------------- #
data <- data.frame(cbind(pc1, pc2))
data$colors <- "#02818a"
data$colors[data$pc1 > 0] <- "#6e016b"
# length(which(data$pc1 > 0)) #89
data$long <- ids_2$long
data$lat <- ids_2$lat

# Export dataset for further use
# setwd("/Users/rpineau3/Dropbox (GaTech)/BioSci-Ratcliff/Rozenn/4.PandoProject/5-Writing/5.Figures/F1/")
# data$ID <- ids_2$ID
# data$Group <- ids_2$Group
# write.csv(data, "PCA_figure1_data.csv", row.names = FALSE)

#setwd("/Users/rpineau3/Dropbox (GaTech)/BioSci-Ratcliff/Rozenn/4.PandoProject/5-Writing/5.Figures/F1") 
#pdf(file=paste("PCA_pando_friends_only.pdf", sep="" ), bg = "transparent", width=5, height=5, family = "Times New Roman")

# Base R plot
par(family="Times New Roman", cex.lab=1.2, cex.axis=1.2, cex.main=1, cex.sub=1)

plot(data$pc1, data$pc2, pch = 16, col = data$colors, 
       xlab = paste("PC1 (", ax1, "%)"), ylab = paste("PC2 (", ax2, "%)"), 
     bty = "n", xaxt='n', yaxt='n')
axis(1, at = c(-30,-20,-10,0,10, 20), lty = 0) # add axes values
axis(2, at = c(-20, -10,0,10, 20) , lty = 0, las = 2) # add axes values

grid(nx = NULL, ny = NULL, lty = 1, col = "gray", lwd = 0.2)  # add grid

points(data$pc1, data$pc2, pch = 16, col = data$colors)
abline(v = 0, lty = 2, col = "grey", lwd = 2)
abline(h = 0, lty = 2, col = "grey", lwd = 2)

#dev.off()

# ------------------------------------------------------------------------------------- #
# K-means clustering
# ------------------------------------------------------------------------------------- #

data_km <- kmeans(t(gl_2), 2)
length(data_km$cluster) #184
cols <- c("#02818a", "#6e016b")
ks <- as.matrix(data_km$cluster)
par(family="Times New Roman", cex.lab=1.2, cex.axis=1.2, cex.main=1, cex.sub=1)

plot(data$pc1, data$pc2, pch = 16, col = cols[ks[,1]], 
     xlab = paste("PC1 (", ax1, "%)"), ylab = paste("PC2 (", ax2, "%)"), 
     bty = "n", xaxt='n', yaxt='n')
axis(1, at = c(-30,-20,-10,0,10, 20), lty = 0) # add axes values
axis(2, at = c(-20, -10,0,10, 20) , lty = 0, las = 2) # add axes values

write(ks, "/Users/rozenn/GaTech Dropbox/CoS/BioSci/BioSci-Ratcliff/Rozenn/4.PandoProject/11-field-lab-work/7-variants/manual_filtering/data/GL/forF1/k_assignment.txt",
      sep = "\t", ncolumns = 1)

# ------------------------------------------------------------------------------------- #
# Mapping
# ------------------------------------------------------------------------------------- #

m <- leaflet()
m <- addTiles(m)
m <- addCircles(m, lng = data$long, lat = data$lat, radius = 10 , opacity = 1, color = data$colors) # fillColor = data$cols,  color = data$cols,
m <-  addScaleBar(m, position = c( "bottomleft"), options = scaleBarOptions(imperial = FALSE))
m
