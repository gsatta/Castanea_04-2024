# Voglio fare una LDA per ogni inoculum
{
library(readr)
library(dplyr)
library(ggplot2)
library(MASS)
library(dplyr)
library(gridExtra)

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
  mutate(combination = paste(inoculum, treatment, sep = "-"))

result_df <- result_df %>%
  relocate(combination, .before = length)

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
  dplyr::select(length, avgDiam, rootVolume, FRL, CRL, FRS, CRS, FVOL, lenPerVol, weight)

# Standardizza le colonne numeriche
df_scaled <- as.data.frame(scale(df_numerici))

# Aggiungi le colonne non numeriche al dataframe standardizzato
df_scaled <- cbind(df[, -(5:14)], df_scaled)

head(df_scaled)

# Elimina i dari mancanti altrimenti non funziona la multivariata

df=na.omit(df)

df_scaled=na.omit(df_scaled)

}

# library(GGally)
# 
# ggscatmat(df_scaled, columns = 5:14, color = "inoculum")
# 
# 
# ggscatmat(df_scaled, columns = 5:14, color = "treatment")


###############################

inocula <- unique(df_scaled$inoculum)
df_list_inc <- list()

for (i in inocula) {
  temp_df <- subset(df_scaled, inoculum == i)
  temp_df$inoculum <- NULL
  df_list_inc[[i]] <- temp_df
}

lda_plots_list_in <- list()

######
library(mda)
library(ggalt)

# Ciclo for per eseguire le operazioni su ciascun dataframe
for (i in names(df_list_inc)) {
  # Ottieni il dataframe corrente
  df <- df_list_inc[[df_name]]
  
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
  
  # Estrai il nome del dataframe corrente
  df_name_clean <- gsub("_df$", "", df_name)
  
  # Plot della PCA
  PCA_ggplot <- ggplot(pca_scores, aes(x = Comp.1, y = Comp.2, col = treatments)) +
    geom_point() +
    ggtitle(paste("PCA", df_name_clean, "Plot"))
  
  # Esegui LDA
  lda_result <- lda(treatments ~ ., data = pca_scores[, c(1:10)])
  
  # Assegna i valori predetti
  pred_lda <- predict(lda_result, dimen = 4)
  valori <- pred_lda$x

  # Combina i valori LDA con i valori PCA
  lda_plot <- cbind(pca_scores, LD1 = valori[,1], LD2 = valori[,2])
  
  # Crea il plot utilizzando ggplot2
  lda_ggplot <- ggplot(lda_plot, aes(LD1, LD2)) +
    geom_point(aes(color = treatments)) +
    geom_encircle(aes(color = treatments), s_shape = 1, expand = 0.05, linetype = "solid") +  # Aggiungi questa linea
    ggtitle(paste(df_name_clean, "- Linear Discriminant Analysis (LDA)")) +
   
     theme(legend.position = "right")  # Nascondi la legenda
  
  # Salva il grafico LDA nella lista
  lda_plots_list_in[[i]] <- lda_ggplot
  
  # Salva il grafico LDA in locale
  ggsave(filename = paste0("./GRAPHS/LDA/INOCULUM/", df_name_clean, "_LDA_plot.png"), plot = lda_ggplot, width = 8, height = 6)
  
  # Salva il grafico PDA in locale
  ggsave(filename = paste0("./GRAPHS/PCA/INOCULUM/", df_name_clean, "_PCA_plot.png"), plot = PCA_ggplot, width = 8, height = 6)
  
  # Stampare esplicitamente il grafico PCA e LDA
  print(PCA_ggplot)
  print(lda_ggplot)
}

# Organizza i grafici LDA in una griglia
lda_grid_in <- do.call(grid.arrange, lda_plots_list_in)

# Salva il grafico della griglia LDA
ggsave(filename = "./GRAPHS/LDA/INOCULUM/LDA_inoculum.png", plot = lda_grid_in, width = 15, height = 12)

##################################
# LDA FOR TREATMENTS
#############################

# Creare una lista vuota per memorizzare i dataframe
treatment <- unique(df_scaled$treatment)
df_list_treat <- list()

# Ciclo attraverso ogni trattamento unico
for (i in treatment) {
  temp_df <- subset(df_scaled, treatment == i)
  temp_df$treatment <- NULL
  df_list_treat[[i]] <- temp_df
}

# Inizializza una lista per memorizzare i grafici LDA
lda_plots_list_treat <- list()

######
library(mda)
library(ggalt)

# Ciclo for per eseguire le operazioni su ciascun dataframe
for (i in seq_along(df_list_treat)) {
  # Ottieni il dataframe corrente
  df_name <- names(df_list_treat)[i]
  df <- df_list_treat[[df_name]]
  
  # Rimuovi la colonna 'treatment'
  df$treatment <- NULL
  
  # Effettua il test di Kaiser
  # print(EFAtools::KMO(df[, -(1:3)], use = "na.or.complete", cor_method = "kendall"))
  
  # Esegue la PCA
  pca <- princomp(df[, -(1:3)], cor = TRUE)
  
  # Plot della PCA
  plot(pca)
  
  # Stampa la varianza spiegata
  print(pca$sdev^2)
  
  # Prendi gli scores direttamente dalla funzione
  pca_scores <- data.frame(inoculum = df$inoculum, pca$scores)
  
  # Estrai il nome del dataframe corrente
  df_name_clean <- gsub("_df$", "", df_name)
  
  # Plot della PCA
  PCA_ggplot <- ggplot(pca_scores, aes(x = Comp.1, y = Comp.2, col = inoculum)) +
    geom_point() +
    ggtitle(paste("PCA", df_name_clean, "Plot"))
  
  # Esegui LDA
  lda_result <- lda(inoculum ~ ., data = pca_scores[, c(1:10)])
  
  # Assegna i valori predetti
  pred_lda <- predict(lda_result, dimen = 4)
  valori <- pred_lda$x
  
  # Combina i valori LDA con i valori PCA
  lda_plot <- cbind(pca_scores, LD1 = valori[,1], LD2 = valori[,2])
  
  # Crea il plot utilizzando ggplot2
  lda_ggplot <- ggplot(lda_plot, aes(LD1, LD2)) +
    geom_point(aes(color = inoculum)) +
    geom_encircle(aes(color = inoculum), s_shape = 1, expand = 0.05, linetype = "solid") +  # Aggiungi questa linea
    ggtitle(paste(df_name_clean, "- Linear Discriminant Analysis (LDA)")) +
    
    theme(legend.position = "right")  # Nascondi la legenda
  
  # Salva il grafico LDA nella lista
  lda_plots_list_treat[[i]] <- lda_ggplot
  
  # Salva il grafico LDA in locale
  ggsave(filename = paste0("./GRAPHS/LDA/TREATMENT/", df_name_clean, "_LDA_plot_", i, ".png"), plot = lda_ggplot, width = 8, height = 6)
  
  ggsave(filename = paste0("./GRAPHS/PCA/TREATMENT/", df_name_clean, "_PCA_plot_", i, ".png"), plot = PCA_ggplot, width = 8, height = 6)
  
  # Stampare esplicitamente il grafico LDA
  print(lda_ggplot)
}

# Organizza i grafici LDA in una griglia
lda_grid_treat <- do.call(grid.arrange, lda_plots_list_treat)

# Salva il grafico della griglia LDA
ggsave(filename = "./GRAPHS/LDA/TREATMENT/LDA_treatments.png", plot = lda_grid_treat, width = 15, height = 12)


