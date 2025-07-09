#haplotype cluster
# =============================================================
# FULL SCRIPT: PCA + K-means + Hierarchical Clustering
# Includes: AMO exclusion, Elbow & Silhouette, and Dendrogram
# =============================================================

# Load required libraries
library(Biostrings)
library(ggplot2)
library(cluster)
library(factoextra)
library(dplyr)

# STEP 1: Read aligned FASTA and extract sequences
seqs <- readDNAStringSet("aligned_mt_overlap.fasta")
pop_names <- names(seqs)

# Convert to character matrix and extract variable sites
seq_matrix <- as.matrix(DNAStringSet(seqs))
rownames(seq_matrix) <- pop_names
var_sites <- apply(seq_matrix, 2, function(col) length(unique(col)) > 1)
seq_var <- seq_matrix[, var_sites]

# Convert to numeric (A=1, C=2, G=3, T=4)
nuc_to_num <- function(x) {
  nuc_map <- c("A"=1, "C"=2, "G"=3, "T"=4, "-"=NA, "N"=NA)
  return(unname(nuc_map[x]))
}
numeric_matrix <- t(apply(seq_var, 1, nuc_to_num))
rownames(numeric_matrix) <- pop_names

# Remove rows with missing values
clean_matrix <- numeric_matrix[complete.cases(numeric_matrix), ]

# STEP 2: PCA
pca <- prcomp(clean_matrix, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Pop <- rownames(clean_matrix)

# Load population region table
region_df <- read.csv("poolinfo.csv")
colnames(region_df)[2:3] <- c("Pop", "Region")
region_df <- region_df[, c("Pop", "Region")]
pca_df <- merge(pca_df, region_df, by = "Pop", all.x = TRUE)

# STEP 3: Elbow method (with AMO)
datamat <- pca_df[, c("PC1", "PC2", "PC3")]
wss <- numeric()
for (k in 1:10) {
  set.seed(123)
  km <- kmeans(datamat, centers = k)
  wss[k] <- km$tot.withinss
}
plot(1:10, wss, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method (with AMO)")

# STEP 4: Silhouette method (with AMO)
fviz_nbclust(datamat, kmeans, method = "silhouette",
             diss = dist(datamat), k.max = 10, nstart = 25) +
  labs(title = "Silhouette (with AMO)")

# STEP 5: Remove AMO, re-run elbow & silhouette
pca_df_filtered <- subset(pca_df, Pop != "AMO")
datamat_filtered <- pca_df_filtered[, c("PC1", "PC2", "PC3")]

# Elbow (no AMO)
wss2 <- numeric()
for (k in 1:10) {
  set.seed(123)
  km <- kmeans(datamat_filtered, centers = k)
  wss2[k] <- km$tot.withinss
}
plot(1:10, wss2, type = "b", pch = 19,
     xlab = "Number of Clusters (k)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method (excluding AMO)")

# Silhouette (no AMO)
fviz_nbclust(datamat_filtered, kmeans, method = "silhouette",
             diss = dist(datamat_filtered), k.max = 10, nstart = 25) +
  labs(title = "Silhouette (excluding AMO)")

# STEP 6: Final PCA plot with k = 5 (AMO included)
k_to_plot <- 5
set.seed(123)
km <- kmeans(datamat, centers = k_to_plot)
pca_df$cluster <- as.factor(km$cluster)

# PCA scatter plot (with small points and population labels)
ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, shape = Region)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_text(aes(label = Pop), hjust = -1, size = 1.5) +
  scale_shape_manual(values = c("Alaska" = 16, "BC" = 17)) +
  labs(
    title = paste("PCA + k-means clustering (k =", k_to_plot, ")"),
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()

# STEP 7: Proportional barplot of Region composition by cluster
cluster_region_prop <- pca_df %>%
  group_by(cluster, Region) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(cluster) %>%
  mutate(prop = count / sum(count))

ggplot(cluster_region_prop, aes(x = cluster, y = prop, fill = Region)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportion of Alaska vs BC in each Cluster (k = 5)",
    x = "Cluster",
    y = "Proportion",
    fill = "Region"
  ) +
  theme_minimal()

# STEP 8: Hierarchical clustering (on PCA space)
dist_mat <- dist(datamat)  # distance matrix on PC1–PC3
hc <- hclust(dist_mat, method = "complete")  # or try "average", "ward.D2"

# Plot dendrogram with region color labels
plot(hc, labels = pca_df$Pop, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "")

# Optional: color branches by Region (simple approach)
region_colors <- as.numeric(as.factor(pca_df$Region))
label_colors <- setNames(region_colors, pca_df$Pop)

# (Advanced coloring using dendextend can be added)

# Run k-means clustering (k = 5)
k_to_plot <- 5
set.seed(123)
km5 <- kmeans(datamat, centers = k_to_plot)
sil <- silhouette(km5$cluster, dist(datamat))

# Add population names
rownames(sil) <- pca_df$Pop  # ← 关键一步：用 Pop 名作为标签

# Plot with population names
fviz_silhouette(sil, label = TRUE) +
  labs(title = paste("Silhouette Plot for k =", k_to_plot))


# =============================================================
# Hierarchical Clustering on Aligned Sequences (No PCA, No k-means)
# =============================================================

# Load required libraries
library(Biostrings)
library(ggplot2)
library(dplyr)

# Step 1: Read aligned consensus sequences
seqs <- readDNAStringSet("aligned_mt_overlap.fasta")
pop_names <- names(seqs)

# Step 2: Convert to character matrix
seq_matrix <- as.matrix(DNAStringSet(seqs))
rownames(seq_matrix) <- pop_names

# Step 3: Compute Hamming distances between all sequences
hamming_distance <- function(s1, s2) {
  sum(s1 != s2)
}

n <- nrow(seq_matrix)
dist_mat <- matrix(0, n, n)
rownames(dist_mat) <- colnames(dist_mat) <- rownames(seq_matrix)

for (i in 1:n) {
  for (j in 1:n) {
    dist_mat[i, j] <- hamming_distance(seq_matrix[i, ], seq_matrix[j, ])
  }
}

# Convert to dist object for clustering
dist_obj <- as.dist(dist_mat)

# Step 4: Perform hierarchical clustering
hc <- hclust(dist_obj, method = "complete")  # Can try "average", "ward.D2"

# Step 5: Plot dendrogram
plot(hc, main = "Hierarchical Clustering of Aligned Sequences",
     xlab = "Samples", sub = "", cex = 0.8)

# Step 6: Sensitivity analysis: Cut tree at various heights
h_values <- seq(5, max(hc$height), by = 1)
cluster_counts <- sapply(h_values, function(h) {
  length(unique(cutree(hc, h = h)))
})

# Step 7: Plot sensitivity analysis
sensitivity_df <- data.frame(h = h_values, clusters = cluster_counts)

ggplot(sensitivity_df, aes(x = h, y = clusters)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(size = 2, color = "darkred") +
  labs(title = "Sensitivity Analysis: Cut Height vs Number of Clusters",
       x = "Cut Height (h)",
       y = "Number of Clusters") +
  theme_minimal()