# 4. Noise Filtering & Dimension Reduction  
# load the data from 3.5
data3.5 <- read.csv("data3_5.csv",sep = ",",header = TRUE,row.names = 1)
# check the number genes
dim(data3.5)
# filter 1:
filter1 = data3.5 > log(15,2)
data4.1 <- data3.5[rowSums(filter1) > .2*134,]
dim(data4.1)

# filter 2
# computer the variation of each gene and find the median
each_gene_var <- apply(data4.1,MARGIN = 1,var) 
all_var_median <- median(each_gene_var)
# select the genes that pass one tail chisq 
filter2 <- 133*(each_gene_var/all_var_median)>qchisq(0.99,133)
data4.2 <- data4.1[filter2,]
# check how many genes left
dim(data4.2)
# save the data for the biologist
write.csv(data4.2,sep = ",",file = "data4_5.csv")

# filter 3
# compute the mean and standard deviation for each gene
data4.2_mean <- apply(data4.2,MARGIN = 1,mean)
data4.2_sd <- apply(data4.2,MARGIN = 1,sd)
# compute the coeffitient of variation
data4.2_cv <- data4.2_sd/data4.2_mean
data4.3 <- data4.2[data4.2_cv>0.186,]
# check how many genes left after 3 filters
dim(data4.3)
# save the data for further analysis
write.csv(data4.3,sep=",",file = "data4_4.csv")
data4.4 <- read.csv("data4_4.csv",header = TRUE, sep = ",",row.names = 1)

# 5. Hierarchical clustering & subtype discovery
clusters <- hclust(dist(t(data4.4)))
# plot the dendrogram
plot(clusters,labels = FALSE,main = "",sub="")
# add a line can separate the sample into two clusters
abline(h=100,lty=2)
clusterCut <- cutree(clusters, k = 2, h = 100)
# summarize how many samples in each cluster
table(clusterCut)  
# create a new vector for separating the C3 and C4 by colors
colorow <- replicate(n=134,NA)
batch_data <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
# label C3 as red, and C4 as blue
colorow <- ifelse(batch_data$cit.coloncancermolecularsubtype=="C3","red","blue")
heatmap(as.matrix(data4.4),ColSideColors = colorow,labRow = FALSE,labCol = FALSE)
legend(x=99,y=142,legend=c("C3","C4"),fill = c("red","blue"))
# summarize the C3 and C4 subtype and each cluster
table(batch_data$cit.coloncancermolecularsubtype)
table(clusterCut)
# perform a fisher test 
fish <- fisher.test(matrix(c(0,75,57,2),byrow = T, nrow = 2))

# compute the t-test for each gene between two clusters  
t_data <- apply(as.matrix(data4.4),1,function(x) t.test(x=x[clusterCut==1],y=x[clusterCut==2]))
# extract the p_val and test statistics
p_value <- sapply(t_data,function(x) x$p.value)
t_stats <- sapply(t_data,function(x) x$statistic)
# compute the q value according to p value
adjust_p <- p.adjust(p_value,method = "fdr")
# build a new dataframe with probeset ID, t stas, p val, q val.
data5.4 <- data.frame("Probeset_ID" = c(row.names(data4.4)),
                      t_stats,p_value,adjust_p)
# save the data in a csv file
write.csv(data5.4,sep = ",",row.names = FALSE,file = "data5_4.csv")
# select the gene with q value less than 0.05
diff_gene <- data5.4$Probeset_ID[data5.4$adjust_p<0.05]
# check how many genes can pass screen
length(diff_gene)

# do the same test to the data from 4.2 for biologist
data4.5 <- read.csv("data4_5.csv",row.names = 1,sep = ",",header = T)
t_data_bio <- apply(as.matrix(data4.5),1,function(x) t.test(x=x[clusterCut==1],y=x[clusterCut==2]))
bio_pvalue <- sapply(t_data_bio,function(x) x$p.value)
bio_t_stats <- sapply(t_data_bio, function(x) x$statistic)
bio_adjustp <- p.adjust(bio_pvalue,method = "fdr")
# create a new dataframe with probeset_ID, t stats, p val, q val
data5.6 <- data.frame("Probeset_ID"=c(row.names(data4.5)),"t_stats"=bio_t_stats,
                      "p_value"=bio_pvalue,"q_value"=bio_adjustp)
# select the genes with q value less than 0.05
data5.6 <- data5.6[data5.6$q_value<0.05,]
# check how many genes can pass the filter
dim(data5.6)
# save the data 
write.csv(data5.6,row.names = F,file = "data5_6.csv",sep=",")