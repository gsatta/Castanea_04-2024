# Voglio fare una LDA per ogni inoculum

library(readr)
library(dplyr)
library(ggplot2)

# passiamo a notazione numerica
options(scipen=999)

#### 1 - creazione di tre dataframe, uno per ogni inoculum

# Carica il df result_df salvato precedentemente
result_df <- read_csv("./DATAFRAMES/result_df_weight.csv")

# Elimina la colonna indesiderata $...1
result_df$...1 <- NULL

# Converti la colonna "plant" in caratteri (stringhe)
result_df$plant <- as.character(result_df$plant)

result_df <- result_df %>%
  mutate(combination = paste(inoculum, treatment, sep = "-")) %>%
  select(plant, inoculum, treatment, combination, everything())

# Crea un nuovo dataframe uguale al precedente
df <- result_df

# df$lenPerVol <- NULL

df$crossings <- NULL

df$tips <- NULL

df$forks <- NULL

# Trasformazione della variabile treatment e inoculum in fattore
df$treatment <- as.factor(df$treatment)
df$inoculum <- as.factor(df$inoculum)
df$combination <- as.factor(df$combination)

# Seleziona solo le colonne numeriche da standardizzare
df_numerici <- df %>%
  select(length, avgDiam, rootVolume, FRL, CRL, FRS, CRS, FVOL, lenPerVol, weight)

# Standardizza le colonne numeriche
df_scaled <- as.data.frame(scale(df_numerici))

# Aggiungi le colonne non numeriche al dataframe standardizzato
df_scaled <- cbind(df[, -(5:14)], df_scaled)

head(df_scaled)

# Converti i valori di NI = 1, AB21 = 2, AB1 = 3
# Assegna i valori desiderati alla colonna Inoculum
df_scaled$inoculum <- ifelse(df_scaled$inoculum == "NI", 1,
                             ifelse(df_scaled$inoculum == "AB21", 2,
                                    ifelse(df_scaled$inoculum == "AB1", 3, df_scaled$inoculum)))


# Elimina i dari mancanti altrimenti non funziona la multivariata

df=na.omit(df)

df_scaled=na.omit(df_scaled)

###############################

ni_df <- subset(df_scaled, inoculum == 1)

ni_df$inoculum <- NULL

ab21_df <- subset(df_scaled, inoculum == 2) 

ab21_df$inoculum <- NULL

ab1_df <- subset(df_scaled, inoculum == 3) 

ab1_df$inoculum <- NULL


######
library(mda)
library(ggalt)

# Definisci la lista dei dataframe
df_list <- list(ni_df = ni_df, ab21_df = ab21_df, ab1_df = ab1_df)

# Ciclo for per eseguire le operazioni su ciascun dataframe
for (df_name in names(df_list)) {
  # Ottieni il dataframe corrente
  df <- df_list[[df_name]]
  
  # Rimuovi la colonna 'inoculum'
  df$inoculum <- NULL
  
  # Effettua il test di Kaiser
  # print(EFAtools::KMO(df[, -(1:3)], use = "na.or.complete", cor_method = "kendall"))
  
  # Esegue la PCA
  pca <- princomp(df[, -(1:3)], cor = TRUE)
  
  # Plot della PCA
  plot(pca)
  
  # Stampa la varianza spiegata
  print(pca$sdev^2)
  
  # Prendi gli scores direttamente dalla funzione
  pca_scores <- data.frame(treatments = df$treatment, pca$scores)
  
  # Plot della PCA
  ggplot(pca_scores, aes(x = Comp.1, y = Comp.2, col = treatments)) +
    geom_point() +
    ggtitle(paste("PCA", unique(df$treatment), "Plot"))
  
  # Esegui LDA
  lda_result <- lda(treatments ~ ., data = pca_scores[, -c(10:11)])
  
  # Assegna i valori predetti
  pred_lda <- predict(lda_result, dimen = 4)
  valori <- pred_lda$x
  
  # Estrai il nome del dataframe corrente
  df_name_clean <- gsub("_df$", "", df_name)
  
  # Combina i valori LDA con i valori PCA
  lda_plot <- cbind(pca_scores, LD1 = valori[,1], LD2 = valori[,2])
  
  # Crea il plot utilizzando ggplot2
  lda_ggplot <- ggplot(lda_plot, aes(LD1, LD2)) +
    geom_point(aes(color = treatments)) +
    geom_encircle(aes(color = treatments), s_shape = 1, expand = 0.05, linetype = "solid") +  # Aggiungi questa linea
    ggtitle(paste(df_name_clean, "- Linear Discriminant Analysis (LDA)")) +
    theme(legend.position = "right")  # Nascondi la legenda
  
  # Stampare esplicitamente il grafico LDA
  print(lda_ggplot)
  
}


