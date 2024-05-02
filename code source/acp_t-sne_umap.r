#Importer les données:
data <- read.csv('C:/Users/asus/OneDrive/Desktop/Projet/data/breast_cancer.csv',
                 sep=';', row.names='id_sample')
dim(data)
head(data)

#Afficher les occurences de chaque niveau de la colonne **pam50**:


group_sizes <- table(data$pam50)
print(group_sizes)


library(ggplot2)
ggplot(data, aes(x = pam50)) +
  geom_bar() +
  labs(title = "Nombre d'occurrences par groupe",
       x = "pam50",
       y = "Nombre d'occurrences")


# 2.Séparer les données d'expression et les étiquettes


library(dplyr)
#Données d'expression de 50 gènes
X <- select_if(data, is.numeric)
cat('X', dim(X), '\n')



#Etiquettes correspondantes (sous-types moléculaires)
y <- data$pam50
cat('y', length(y), '\n')

#Afficher les valeurs d'expression:

library(reshape2)
library(tidyr)
#Triez les colonnes par moyenne
sort_by_mean <- colMeans(X)
sort_by_mean <- sort_by_mean[order(sort_by_mean)]
#Créez un data.frame trié
X_sorted <- X[, names(sort_by_mean)]
#Créez le graphique en boîte avec ggplot2
ggplot(melt(X_sorted), aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "Graphique en boîte trié par moyenne",
       x = "Gènes",
       y = "Expression") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# 4.Appliquer une normalisation centrée-réduite aux données:


#Normalisation centrée-réduite (z-score)
X_scaled <- scale(X)
#Convertir la matrice résultante en un data.frame avec les mêmes noms de colonnes et indices
X_scaled <- as.data.frame(X_scaled)
colnames(X_scaled) <- colnames(X)
rownames(X_scaled) <- rownames(X)
#head(X_scaled)
ggplot(data = gather(X_scaled), aes(x = key, y = value)) +
  geom_boxplot() +
  labs(title = "Graphique en boîte des données normalisées",
       x = "Gènes",
       y = "Expression") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#Faire une analyse en composantes principales (ACP)
##Calcul de l’ACP:

#Installez et chargez les bibliothèques FactoMineR et factoextra si ce n'est pas déjà fait
#install.packages("FactoMineR")
#install.packages("factoextra")
library(FactoMineR)
library(factoextra)

#Réalisez l'ACP sur les données X_scaled
pca_result <- PCA(X_scaled, graph = FALSE)

#Obtenez les composantes principales (PC)
X_pca <- as.data.frame(pca_result$ind$coord)

#Renommez les colonnes des PC
colnames(X_pca) <- paste0("PC", 1:ncol(X_pca))

#Affichez les premières lignes du data.frame avec les composantes principales
head(X_pca)

##Calcul de la variance expliquée:

pca <- prcomp(X_scaled)
explained_variance <- tibble(
  PC = seq_along(pca$sdev),
  Explained_Variance_Ratio = pca$sdev^2 / sum(pca$sdev^2) * 100
)
# Afficher la variance expliquée
print(explained_variance)

head(explained_variance)

#Créer le graphique à barres avec ggplot2
ggplot(data = explained_variance, aes(x = factor(PC), y = Explained_Variance_Ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Variance expliquée par chaque composante principale",
       x = "Composante Principale",
       y = "Pourcentage de Variance Expliquée") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



#Visualisation des deux premières composantes principales de l’ACP
library(ggplot2)
# Create a scatter plot of PC1 vs PC2
ggplot(X_pca, aes(x = PC1, y = PC2)) +
  geom_point(color = 'gray') +
  labs(title = 'Scatter Plot des PC1 et PC2', x = 'PC1', y = 'PC2')
```
Les points se rassemblent naturellement en clusters, notamment deux clusters sont nettement visibles. Présentons chaque sous-type moléculaire en différente couleur pour comprendre à quoi correspondent ces clusters.


library(ggplot2)

# Ajouter les composantes principales au DataFrame data
data_with_pca <- cbind(data, PC1 = X_pca$PC1, PC2 = X_pca$PC2)

# Créer une palette de couleurs pour les classes
dict_colors <- c('luminal-A' = 'forestgreen', 'luminal-B' = 'royalblue', 'HER2-enriched' = 'orange', 'basal-like' = 'crimson')
data_with_pca$color <- dict_colors[data_with_pca$pam50]

# Créer le graphique en nuage de points
ggplot(data_with_pca, aes(x = PC1, y = PC2, color = color)) +
  geom_point() +
  labs(title = "Graphique en nuage de points avec couleurs de classe", x = "PC1", y = "PC2") +
  theme_minimal()


#Visualisation des trois premières composantes principales de l’ACP:



#Create 3D scatter plot
install.packages("rgl")
library(rgl)
plot3d(X_pca[, 1], X_pca[, 2], X_pca[, 3], col = y_colors, size = 30, xlab = paste0("PC1 - ", explained_variance[1], "%"),
ylab = paste0("PC2 - ", explained_variance[2,], "%"), zlab = paste0("PC3 - ", explained_variance[3,], "%"))
rgl.viewpoint(theta = 45, phi = 15)



#Visualiser les données avec la méthode t-SNE:
## Projection 2D:




#Installer la bibliothèque si ce n'est pas déjà fait
if (!requireNamespace("Rtsne", quietly = TRUE)) {
  install.packages("Rtsne")
}

#Charger la bibliothèque
library(Rtsne)

#Définir les paramètres de t-SNE
tsne_params <- Rtsne(X_scaled, check_duplicates = FALSE, pca = TRUE)

#Obtenir les coordonnées t-SNE
X_tsne <- tsne_params$Y

#Renommer les colonnes
colnames(X_tsne) <- paste0("DIM", 1:ncol(X_tsne))

#Afficher les premières lignes
head(X_tsne)

#Convertir la matrice en data.frame
X_tsne_df <- as.data.frame(X_tsne)

#Ajouter la classe y à X_tsne_df
X_tsne_df$y <- y

#Créer le graphique de dispersion t-SNE en 2D avec ggplot2
tsne_plot <- ggplot(X_tsne_df, aes(x = DIM1, y = DIM2, color = y)) +
  geom_point() +
  labs(title = "Graphique de dispersion t-SNE",
       x = "Dimension 1",
       y = "Dimension 2",
       color = "Classe")

#Afficher le graphique
print(tsne_plot)


#Charger le package Rtsne
library(Rtsne)

tsne_result <- Rtsne(X_scaled, dims = 2)

#Accéder à la divergence de Kullback-Leibler
kl_divergence <- tsne_result$itercost

#Obtenir le résultat de la dernière itération
resultat_derniere_iteration <- kl_divergence[length(kl_divergence)]

#Afficher le résultat de la dernière itération
print(resultat_derniere_iteration)

```

#Visualiser les données avec la méthode UMAP:
## Projection 2D:



library(umap)

# Appliquer UMAP
set.seed(0)  # Pour la reproductibilité
embedding <- umap(X_scaled, n_components = 2, random_state = 0, n_threads = -1)

# Créer un dataframe
X_umap <- data.frame(DIM1 = embedding$layout[, 1], DIM2 = embedding$layout[, 2])

dict_colors <- c('luminal-A' = 'forestgreen', 'luminal-B' = 'royalblue', 'HER2-enriched' = 'orange', 'basal-like' = 'red')
y_colors <- as.character(sapply(y, function(yi) dict_colors[yi]))

# Tracer le graphique de dispersion
plot(X_umap$DIM1, X_umap$DIM2, col = y_colors, xlab = "DIM1", ylab = "DIM2", main = "UMAP Plot")

ggplot(X_umap, aes(x = DIM1, y = DIM2, color = as.factor(y))) +
  geom_point() +
  labs(title = "UMAP Plot",
       x = "DIM1",
      y = "DIM2")



