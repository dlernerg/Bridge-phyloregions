library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(picante)
library(ape)
library(phytools)
library(RColorBrewer)
library(sf)
library(caret)
library(randomForest)
library(boot)

sf_use_s2(F)

source("functions.R")

###########################check how ido's paper has all the updated info

clus.globe <- wwf_simpl.b %>% 
  rename(cells = 1) %>%
  mutate(cluster = cells)
GBIF_overlap <- st_intersects(global.inland,GBIFpolygon_groups)


# Make phylotree ----------------------------------------------------------

tree_Sanchez_etal<-read.tree(file = "SanchezM (2020).tre")

tree <- tree_Sanchez_etal
my_data.spec <- GBIFpolygon_groups$specie %>% as.data.frame()
my_data.spec$genus <- stringr::word(my_data.spec[,1],1)

genus <- tree$tip.label
GBIF.indata <- which(my_data.spec$genus %in% genus)

genus.absent <- genus[which(!(genus %in% unique(my_data.spec$genus)))]
genus.present <- genus[which((genus %in% unique(my_data.spec$genus)))]
tree.new <- drop.tip(tree,genus.absent)


# Calculate complete pairwise phylogenetic turnover for all world grid cells --------------------------------------------

GBIF_overlap.m <- reshape2::melt(GBIF_overlap)
colnames(GBIF_overlap.m) <- c("GBIF.spec","position")
GBIF_overlap.m$value <- 1
GBIF_overlap.m$genus <- GBIFpolygon_groups$specie[GBIF_overlap.m$GBIF.spec] %>% stringr::word(.,1)
GBIF_overlap.m2 <- GBIF_overlap.m[,2:4] %>% group_by(position,genus) %>% distinct() %>% ungroup()
GBIF_overlap.m2$value <- 1
head(GBIF_overlap.m2)

spec.df_long <- GBIF_overlap.m2 %>%
  tidyr::pivot_wider(names_from = "position", values_from = "value",values_fn = sum) %>%
  as.data.frame()
rownames(spec.df_long) <- spec.df_long$genus
colnames(spec.df_long)
spec.df_long = spec.df_long[,2:ncol(spec.df_long)]

spec.df_long.clean <- spec.df_long
spec.df_long.clean[is.na(spec.df_long)] = 0
match.df <- spec.df_long.clean[match(phy$tip.label, rownames(spec.df_long.clean)),]

match.mat <- (match.df) %>% as.matrix()
match.mat.t <- t(match.mat)

raoD.val <- raoD(match.mat.t,phy)

# Maniputalte Matrix ------------------------------------------------------

globe.sample <- c(1:5851)

load("C:/Users/davidle.WISMAIN/Box/lab folder/OTT/phylogenetic trees/phylo dissimilarity/organized/raoD.all.RData")

dis.mat <- raoD.val$H
dis.mat2 <- delete.na(dis.mat) %>% as.matrix()
dis.mat2[which(is.infinite(dis.mat2))] = 1
#2. filter out globe cells without species 
a <- which(rowSums(dis.mat2, na.rm = T) != 0)

globe.sample <- c(1:5851)
cells.num <- globe.sample[as.numeric(names(a))]

#3. find which species and the number of species in each of the cells (num)
num <- find.num(cells.num)
num.f <- 5

dis.mat.f <- dis.mat2[-which(num <= num.f),-which(num <= num.f)]
cells.num.f <- cells.num[- which(num <= num.f)]

dis.mat_log <- log(dis.mat.f)
hc = hclust(dist(dis.mat_log), method = "complete")
k = 9
sub_grp <- cutree(hc, k = k)

#4. cut the tree given discrete k clusters

clus.globe <- data.frame(cells = cells.num.f, cluster = factor(sub_grp))


#5. visualize the different clusters

col.c=c("black", brewer.pal(12, "Set3"))
countries <- map_data("world")

k9 = ggplot() + 
  geom_polygon(data=countries, aes(long, lat, group=group), fill = "ivory2", col = "grey") + 
  geom_sf(global.inland[clus.globe$cells,], mapping = aes(fill = clus.globe$cluster), color = "black")+
  #scale_fill_manual(values = col.c[2:(k+1)])+
  theme_classic()


#6. Assign species to each cluster 

clus.choose <- clus.globe
spec.df <- find.spec(clus.choose)
spec.df$value <- 1

spec.df_long <- spec.df %>%
  tidyr::pivot_wider(names_from = "cluster", values_from = "value",values_fn = sum) %>%
  as.data.frame()
rownames(spec.df_long) <- spec.df_long$.
spec.df_long = spec.df_long[,2:ncol(spec.df_long)]

clus.globe2 <- clus.choose[,2] %>% as.data.frame()
names(clus.globe2) <- "cluster"

cluster.nm <- clus.globe2 %>%
  rstatix::convert_as_factor(cluster) %>%
  rstatix::reorder_levels(cluster, order = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14")) %>%
  group_by(cluster) %>%
  count() %>%
  ungroup()

spec.df.norm <- residues.fun_expected.corrected(spec.df_long)
spec.df.norm.t <- t(spec.df.norm)


# create a matrix for analysis --------------------------------------------------------------

trait <- spec.df.norm %>% as.matrix()
phy <- tree.new

trait[which(is.na(trait))] = 0
match.df <- trait[match(phy$tip.label, rownames(trait)),]

isna <- which(is.na(rowSums(match.df)))
if(is.integer0(isna)){
  match.df <- match.df
}else{
  phy <- drop.tip(phy, tip = phy$tip.label[isna])
  match.df <- match.df[-isna,]
}


match.df <- t(match.df)


# alpha MPD ---------------------------------------------------------------

#create a phylogenetic distance matrix
phydist <- cophenetic.phylo(phy)
mpd <-  ses.mpd(match.df,phydist, null.model = "taxa.labels",runs = 500, abundance.weighted = T)
mpd$cluster <- mapping

mpd.plot_avg <- avg.df.forplot(mpd,"mpd.obs.z")

mpd$significance <- ifelse(mpd$mpd.obs.p < 0.05, "clustered", ifelse(mpd$mpd.obs.p>0.95,"dispersed", "NS"))

mpd.plot <- ggplot(mpd, aes(clus,mpd.obs.z))+
  geom_point(aes(col = significance)) +
  geom_hline(yintercept = -1.64)+
  geom_hline(yintercept = 1.64)+
  scale_color_manual(values = c("red","black","grey"))+
  #scale_x_continuous(labels = c(1:k), breaks = c(1:k)) +
  theme_classic()+
  theme(legend.position = "none")

# BETA: Distance between clusters -----------------------------------------------

phydist.mean <- data.frame()
phydist <- cophenetic.phylo(phy)

match.mat <- t(match.df) %>% as.matrix()
match.mat[which(is.na(match.mat))] = 0

match.mat.t <- t(match.mat)
raoD.val <- raoD(match.mat.t,phy)
gplots::heatmap.2(raoD.val$H, trace = "none", key = F)

raoD.val.pcoa <- ape::pcoa(raoD.val$H)

raoD.pcoa.df <- raoD.val.pcoa$vectors %>% as.data.frame()


ggpubr::ggscatter(raoD.pcoa.df, x = "Axis.1", y = "Axis.2", 
                  label = colnames(raoD.val$H),
                  #label = names_biome$names,
                  size = 1,
                  repel = T,
                  title = "PCoA Phylogenetic Distance")


# Random Forest -----------------------------------------------------------

load("meansBioClim.RData")
load("altitude.RData")
load("P_olsen_df.RData")

meansBioClim2$bio20 <- altitude$m
meansBioClim2$P <- P_olsen_df
meansBioClim3 <- meansBioClim2[complete.cases(meansBioClim2),]

df.random.forest <- clus.choose
lat_lon <- st_coordinates(global.inland$cent[df.random.forest$cells,]) %>% as.data.frame()
df.random.forest$lon <- lat_lon$X
df.random.forest$lat <- lat_lon$Y
df.random.forest <- cbind(df.random.forest,meansBioClim2[df.random.forest$cells,])
df.random.forest <- df.random.forest[complete.cases(df.random.forest),]

training_indices <- createDataPartition(df.random.forest$cluster, p = 0.8, list = FALSE, times = 1)

training_set <- df.random.forest[training_indices, ]
testing_set <- df.random.forest[-training_indices, ]


global_model <- randomForest(cluster ~  bio1 + bio3 + bio4 + bio12 + bio20 + P + lat + lon, data = training_data)
global_pred <- predict(global_model, testing_data)
global_accuracy <- sum(global_pred == testing_data$cluster) / length(global_pred)

feature_importance <- importance(global_model)
varImpPlot(global_model)
partialPlot(x = global_model,pred.data = df.random.forest,x.var = "lat")

accuracies <- data.frame(model = c("env","geo","global"), accuracies = c(env_accuracy,geo_accuracy,global_accuracy))
ggplot(accuracies, aes(x = model, y = accuracies)) + 
  geom_point() + 
  theme_minimal() +
  labs(title = "Model Accuracy vs Number of Trees",
       x = "Number of Trees",
       y = "Accuracy")

conf_matrix <- confusionMatrix(global_pred, testing_data$cluster)
print(conf_matrix)

# Calculating class-specific metrics
precision <- conf_matrix$byClass[, "Precision"]
recall <- conf_matrix$byClass[, "Recall"]
F1 <- conf_matrix$byClass[, "F1"]
balanced <- conf_matrix$byClass[, "Balanced Accuracy"]
# Print class-specific metrics
print(precision)
print(recall)

group_names <- c(1:9) %>% as.character()

metrics_data <- data.frame(
  Group = rep(group_names, 4),
  Metric = rep(c("Precision", "Recall", "F1-Score", "Balanced"),each = length(group_names)),
  Value = c(precision, recall, F1, balanced)
)

metrics_data <- metrics_data[-which(metrics_data$Metric == "Balanced"),]
ggplot(metrics_data, aes(x = Group, y = Value, group = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Metric, scales = "free") +
  theme_minimal() +
  labs(y = "Metric Value", x = "Metric", title = "Performance Metrics by Evolutionary Group") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



n_bootstrap <- 100
bootstrap_results <- vector("list", n_bootstrap)

for(i in 1:n_bootstrap) {
  boot_indices <- sample(1:nrow(df.random.forest), replace = TRUE)
  bootstrap_data <- df.random.forest[boot_indices, ]
  
  training_indices <- createDataPartition(y = bootstrap_data$phylo, p = 0.8, list = FALSE)
  
  training_indices <- sample(1:nrow(bootstrap_data), nrow(bootstrap_data) * 0.8)
  training_data <- bootstrap_data[training_indices, ]
  testing_data <- bootstrap_data[-training_indices, ]
  
  boot_model <- randomForest(phylo ~ bio1 + bio3 + bio4 + bio12 + bio20 + P + lat + lon, data = training_data)
  
  global_pred <- predict(boot_model, testing_data)
  conf_matrix <- confusionMatrix(global_pred, testing_data$phylo)
  
  bootstrap_results[[i]] <- list(
    Model = boot_model,
    ConfusionMatrix = conf_matrix
  )
}

results_df <- data.frame(Iteration = integer(),
                         Precision = numeric(),
                         Recall = numeric(),
                         F1 = numeric(),
                         class = character()
)

for (i in seq_along(bootstrap_results)) {
 
  cm <- bootstrap_results[[i]]$ConfusionMatrix
  precision <- cm$byClass[,"Precision"]
  recall <- cm$byClass[,"Recall"]
  f1 <- cm$byClass[,"F1"]
  iteration_data <- data.frame(
    Iteration = rep(i, length(precision)),  
    class = names(precision),               
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
  # Bind this iteration's data to the results dataframe
  results_df <- rbind(results_df, iteration_data)
}

results_long <- results_df_phylo %>%
  pivot_longer(cols = -c(Iteration,class), names_to = "Metric", values_to = "Value")

# Generate the plot
plot <- ggplot(results_long, aes(x = Iteration, y = Value, color = Metric)) +
  geom_line() + 
  geom_smooth(method = "loess", se = FALSE) +
  facet_grid(Metric ~ class, scales = "free_y") +  
  labs(title = "Bootstrap Model Performance Metrics",
       x = "Bootstrap Iteration",
       y = "Metric Value") +
  theme_minimal()

ggsave(filename = "bootstrap_raw_phylo.svg", width = 14, height = 9)

results_long_f <- results_long[complete.cases(results_long),]
results_summary <- results_long_f %>%
  group_by(Metric, class) %>%
  summarize(
    Mean = mean(Value),
    LowerCI = quantile(Value, probs = 0.025),
    UpperCI = quantile(Value, probs = 0.975),
    .groups = 'drop'
  )

plot_summary <- ggplot(results_summary, aes(x = class, y = Mean, color = Metric)) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2) +
  geom_point() +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "Bootstrap Metric Summaries by Class",
       x = "Class",
       y = "Metric Mean Value") +
  theme_minimal()

ggsave(filename = "bootstrap_interval_phylo.svg", width = 8, height = 5)
