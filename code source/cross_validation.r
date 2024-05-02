#Importer les données
data <- read.csv('C:/Users/asus/OneDrive/Desktop/Projet/data/breast_cancer.csv',
                 sep=';', row.names='id_sample')
head(data)


## Afficher les dimensions des données
#Définir la graine pour la reproductibilité
set.seed(123) 
cat('data', dim(data),'\n')

#Afficher les occurences de chaque niveau de la colonne pam50
group_sizes <- table(data$pam50)
print(group_sizes)

#Séparer les données d’expression et les étiquettes
library(dplyr)
#Données d'expression de 50 gènes
X <- select_if(data, is.numeric)
cat('X', dim(X), '\n')

#Extraire une colonne
#Etiquettes correspondantes (sous-types moléculaires)
y <- data$pam50
cat('y', length(y), '\n')

#Créer une validation croisée:
#Initialisation
 library(ggplot2)
 library(lattice)
 library(caret)
 num_folds <- 3
 folds <- createFolds(y, k = num_folds, list = TRUE, returnTrain = TRUE)

#Boucler sur les plis
  for (fold in 1:num_folds) {
    #Extraire les indices d'entraînement et de test
    train_indices <- unlist(folds[[fold]])
    test_indices <- setdiff(seq_along(y), train_indices)
    #Données d'entraînement et de test
    X_train <- X[train_indices, , drop = FALSE]
    X_test <- X[test_indices, , drop = FALSE]
    #Afficher les résultats pour chaque itération
    cat('Fold', fold, 'Train:', dim(X_train), 'Test:', dim(X_test), '\n')
  }


library(e1071)
 #Définir la graine pour la reproductibilité
 set.seed(0)
 #Créer une fonction pour calculer l'accuracy
 accuracy <- function(y_true, y_pred) {
   mean(y_true == y_pred)
 }

y <- factor(y)
#Initialiser le modèle SVM
classifier <- svm(pam50 ~ ., data = data.frame(X, pam50 = y),
                  kernel = "linear", class.weights = NULL)
#Initialiser la validation croisée
num_folds <- 3
folds <- createFolds(y, k = num_folds, list = TRUE, returnTrain = TRUE)
#Initialiser un vecteur pour stocker les précisions
accuracy_vector <- numeric(length = length(folds))

#Boucler sur les plis  
   for (iteration in seq_along(folds)) {
     #Extraire les indices d'entraînement et de test
      train_indices <- unlist(folds[[iteration]])
      test_indices <- setdiff(seq_along(y), train_indices)
   
     #Données d'entraînement et de test
     X_train <- X[train_indices, , drop = FALSE]
     y_train <- y[train_indices]
     X_test <- X[test_indices, , drop = FALSE]
     y_test <- y[test_indices]
   
     #Normaliser les données
     scaler <- scale
     X_train_scaled <- scaler(X_train)
     X_test_scaled <- scaler(X_test, center = 
                  attr(X_train_scaled,  "scaled:center"), scale =  
                  attr(X_train_scaled, "scaled:scale"))
   
     #Entraîner le modèle SVM
     model <- svm(x = X_train_scaled, y = y_train,
                  kernel = "linear", class.weights = NULL)
   
     #Prédiction sur les données de test
     y_pred_test <- predict(model, newdata = X_test_scaled)
   
     #Calculer l'accuracy
     accuracy_vector[iteration] <- accuracy(y_test, y_pred_test)
   
     #Afficher les résultats pour chaque itération
     cat('Iteration', iteration,
         'Accuracy =', format(accuracy_vector[iteration],
          digits = 8), '\n')
   }



  #Afficher la moyenne des précisions
  cat('Mean accuracy', format(mean(accuracy_vector), digits = 3), '\n')


  

  #Créer un modèle SVM avec prétraitement
   svm_model <- svm(pam50 ~ ., data = data.frame(X, pam50 = y),
                    kernel = "linear", class.weights = NULL)
  #Créer un pipeline avec le modèle SVM et le prétraitement
  pipeline <- caret::train(
   pam50 ~ .,
   data = data.frame(X, pam50 = y),
   method = "svmLinear",
   trControl = trainControl(method = "cv", number = 3),  # 3-fold CV
   preProcess = c("center", "scale"),  #Centrer et réduire les données
   class.weights = "balanced"  #Poids de classe équilibrés
  )
print(pipeline)


