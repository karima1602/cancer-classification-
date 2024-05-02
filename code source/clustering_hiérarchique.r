
#Importer les données:


data2 <-  read.csv('C:/Users/asus/OneDrive/Desktop/Projet/data/lung_cancer.csv',
                   sep=';', row.names='id_sample')
dim(data2)



head(data2)


# Extraire les données d’expression et les étiquettes

#Données d'expression de 50 gènes
X <- subset(data2, select = sapply(data2, is.numeric))
print(paste("X", dim(X)[1], dim(X)[2], sep=" "))




#Etiquettes des échantillons
y <- data2$class
cat("y", length(y), length(unique(y)),
    "['", paste(unique(y), collapse = "' '"), "']", sep=" ")

# Appliquer une normalisation centrée-réduite:

#Charger la bibliothèque
library(stats)

#Normalisation centrée-réduite
X_scaled <- scale(X)
X_scaled <- as.data.frame(X_scaled)
rownames(X_scaled) <- rownames(X)
colnames(X_scaled) <- colnames(X)
print(paste("X_scaled", dim(X_scaled)))

#Définir les couleurs pour chaque classe:


#Définir les couleurs pour chaque classe
class_color <- c('NTL' = 'aquamarine',
                 'ADK' = 'darkslateblue',
                 'SQC' = 'darkorange')

#Etiquettes des échantillons
y <- data2$class
#print(paste("y", length(y), unique(y)))

# Convertir les étiquettes en couleurs
y_color <- class_color[y]


#Réaliser un clustering hiérarchique:


#Installer la bibliothèque si ce n'est pas déjà fait
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

#Charger la bibliothèque
library(pheatmap)

#Paramètres du clustering hiérarchique
metric <- 'euclidean'
method <- 'ward.D2'

sample_rows <- sample(rownames(X_scaled), size = 25)
sample_cols <- sample(colnames(X_scaled), size = 15)
X_scaled_ <- X_scaled[sample_rows, sample_cols]
#Créer une matrice de données transposée (car pheatmap utilise les colonnes pour le clustering)
X_scaled_transposed <- t(X_scaled)

#Assurer que y_color a la bonne longueur (le nombre de colonnes dans X_scaled)
y_color <- rep(y_color, each = ncol(X_scaled))

#Créer le heatmap avec clustering hiérarchique
clustergrid <- pheatmap(
  X_scaled_transposed,
  clustering_distance_rows = metric,
  clustering_method = method,
  col = y_color,
  width = 18,
  height = 17
)




#Récupérer la figure
fig <- clustergrid

#Sauvegarder la figure en tant qu'image PNG
ggsave(file = "C:/Users/asus/OneDrive/Desktop/Projet/hierarchical_clustering2.png",
plot = fig,dpi = 300, width = 8, height = 7, units = "in", device = "png")





