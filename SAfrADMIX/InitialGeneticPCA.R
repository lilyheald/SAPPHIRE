rm(list = ls())
# Set working directory
setwd("~/sapphire/SAfrADMIX")

system("./plink --bfile SAfrADMIX --nonfounders --no-pheno         --allow-no-sex --recode --out ADAPTmap_TOP")

system("./plink --allow-no-sex --nonfounders --file ADAPTmap_TOP --distance-matrix --out dataForPCA")

# PCA
# Load data
dist_populations<-read.table("dataForPCA.mdist",header=F)

# Extract ethnicity names
fam <- data.frame(famids=read.table("dataForPCA.mdist.id")[,1])
# Extract individual names
famInd <- data.frame(IID=read.table("dataForPCA.mdist.id")[,2])

mds_populations <- cmdscale(dist_populations,eig=T,5)

# Extract the eigenvectors
eigenvec_populations <- cbind(fam,famInd,mds_populations$points)
eigenvec_populations <- eigenvec_populations%>%
  mutate(Ethnicity = substr(famids, start = 1, stop = 2))

# Proportion of variation captured by each eigenvector
eigen_percent <- round(((mds_populations$eig)/sum(mds_populations$eig))*100,2)

ggplot(data = eigenvec_populations) +
  geom_point(mapping = aes(x = `1`, y = `2`,color = Ethnicity), show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "SAPPHIRE genetic pca",
       x = paste0("Principal component 1 (",eigen_percent[1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[2]," %)")) +
  theme_minimal()

