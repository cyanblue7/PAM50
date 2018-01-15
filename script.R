library(genefu);
library(reshape2);
library(ggdendro);

root.dir <- 'PAM50/Input';

setwd(root.dir);
mrna <- read.delim('pam50_tcga_breast_dataset.txt', as.is = TRUE);
data(pam50);

colnames(mrna)['NUF2'  == colnames(mrna)] <- 'CDCA1';
colnames(mrna)['NDC80' == colnames(mrna)] <- 'KNTC2';

rownames(mrna) <- mrna$sampleID;

# mrna.pam50 <- intrinsic.cluster(
# 	data=as.matrix(mrna[,pam50$centroids.map$probe]),
# 	annot=pam50$centroids.map,
# 	do.mapping=TRUE,
# 	std="robust",
# 	intrinsicg=pam50$centroids.map[ ,c("probe", "EntrezGene.ID")],
# 	number.cluster=5,
# 	mins=5,
# 	method.cor="spearman",
# 	method.centroids="mean",
# 	verbose=TRUE
# 	);

PAM50Preds <- molecular.subtyping(
	sbt.model = "pam50",
	data=as.matrix(mrna[,pam50$centroids.map$probe]),
	annot=pam50$centroids.map,
	do.mapping=TRUE
	);

# run clustering
pam50.dendro <- as.dendrogram(hclust(d = dist(x = t(mrna[,pam50$centroids.map$probe]))));

# Create dendrogram plot
dendro.plot <- ggdendrogram(
	data = pam50.dendro, rotate = TRUE
	) + theme(axis.text.y = element_text(size = 6));

# melt data to format for ggplot
mrna.long <- melt(t(mrna[,pam50$centroids.map$probe]));

# extract the column order in the dendrogram
pam50.order <- order.dendrogram(pam50.dendro)

# order the levels according to their position in the cluster
mrna.long$Var1 <- factor(
	x = mrna.long$Var1,
	levels = mrna.long$Var1[pam50.order], 
	ordered = TRUE
	);

# create heatmap plot
heatmap.plot <- ggplot(
	data = mrna.long, aes(x = Var2, y = Var1)
	) + geom_tile(aes(fill = value)) +
	scale_fill_gradient2() +
	theme(axis.text.y = element_blank(),
	axis.title.y = element_blank(),
	axis.ticks.y = element_blank(),
	axis.text.x = element_blank(),
	axis.title.x = element_blank(),
	axis.ticks.x = element_blank(),
	legend.position = "top")

# plot dendrogram and heatmap together
library("grid")
setwd(file.path(root.dir, '../Output'));
png('temp.png', res = 600, width = 6, height = 6, unit = 'in');
grid.newpage();
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 0.98))
print(dendro.plot, vp = viewport(x = 0.87, y = 0.416, width = 0.2, height = 0.98))
dev.off();
