## ----setup, include=FALSE--------------------------------------------------------------------
library(knitr)
# knitr options
knitr::opts_chunk$set(echo = FALSE, results = 'show', message = FALSE, warning = FALSE)


## --------------------------------------------------------------------------------------------
#Leer datos
intestinal_data <- read.table("ST000002_AN000002_clean.csv", header=TRUE, sep="\t", row.names = 1)
info <-readLines("ST000002_AN000002_dataset_info.md")
metabolite_metadata <- read.table("ST000002_AN000002_metadata.txt", header=TRUE, sep="\t", row.names = 1)


## --------------------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")


## --------------------------------------------------------------------------------------------
#Cargar paquete `SummarizedExperiment`
library(SummarizedExperiment)

#Convertir los datos de metabolitos a una matriz numérica, convirtiendo los datos a tipo numérico
metabolite_data <- as.matrix(as.data.frame(lapply(intestinal_data[-1,], as.numeric)))

#Obtener los metadatos de las muestras (grupo al que pertenecen) a partir del dataset
sample_metadata <- data.frame(
  Groups = t(intestinal_data[1,]), #Transponer la primera fila para obtener el vector de grupos
  row.names = colnames(intestinal_data)
)

#Crear el objeto `SummarizedExperiment`
se <- SummarizedExperiment(
  assays = list(counts = metabolite_data), #Matriz de datos
  rowData = metabolite_metadata, #Metadatos de los metabolitos
  colData = sample_metadata #Metadatos de las muestras
)

#Agregar la información adicional sobre el dataset en 'metadata' del objeto `SummarizedExperiment`
metadata(se)$dataset_info <- paste(info, collapse = "\n")

#Mostrar el objeto `SummarizedExperiment` resultante
se



## --------------------------------------------------------------------------------------------
#Guardar en un archivo .Rda
save(se, file = "Contenedor_SummarizedExperiment.Rda")


## --------------------------------------------------------------------------------------------
#Clase del objeto `SummarizedExperiment`
class(se)


## --------------------------------------------------------------------------------------------
#Dimensiones del objeto `SummarizedExperiment`
dim(se)


## --------------------------------------------------------------------------------------------
#Matriz de datos del objeto `SummarizedExperiment`
head(assays(se)$counts)


## --------------------------------------------------------------------------------------------
#Metadatos de las filas del objeto `SummarizedExperiment`
rowData(se)


## --------------------------------------------------------------------------------------------
#Metadatos de las columnas del objeto `SummarizedExperiment`
colData(se)


## --------------------------------------------------------------------------------------------
#Metadatos del dataset del objeto `SummarizedExperiment`
metadata(se)


## --------------------------------------------------------------------------------------------
#Obtener el texto completo de `dataset_info`
info_text <- metadata(se)$dataset_info

#Dividir el texto en líneas
info_lines <- strsplit(info_text, "\n")[[1]]

# Filtrar líneas que contienen pares clave-valor
info_kv_lines <- grep(":", info_lines, value = TRUE)

#Separar cada línea en clave y valor usando `strsplit`
info_split <- strsplit(info_kv_lines, "\t|: ")

#Crear un `data.frame` con claves y valores
info_table <- data.frame(
  Clave = sapply(info_split, `[`, 1),
  Valor = sapply(info_split, function(x) paste(x[-1], collapse = " ")),
  stringsAsFactors = FALSE
)

#Ver la tabla generada
knitr::kable(info_table, format = "pipe")



## --------------------------------------------------------------------------------------------
#Obtención de las estadísticas básicas de cada muestra (redondeado)
round(apply(assays(se)$counts,2, summary))


## --------------------------------------------------------------------------------------------
#Histogramas de expresión de cada muestra
par(mfrow=c(3,3)) #Mostrar 3x3 histogramas por ventana
for (i in 1:ncol(assays(se)$counts))
  hist(assays(se)$counts[,i], main = colnames(assays(se)$counts)[i], col = "lightgreen")


## --------------------------------------------------------------------------------------------
#Gráfico de densidad de todas las muestras
par(mfrow = c(1, 1))  #Restablece la ventana de gráficos a una sola gráfica
plot(density(assays(se)$counts[, 1]), main = "Gráfico de densidad de expresión de los metabolitos en las muestras",
     xlab = "Expresión de metabolitos", ylab = "Densidad")
for (i in 2:ncol(assays(se)$counts)) {
  lines(density(assays(se)$counts[, i]), col = i)
}
#Leyenda para el gráfico de densidad
legend("topright", legend = colnames(assays(se)$counts), col = 1:ncol(assays(se)$counts), lwd = 2, cex = 0.7)



## --------------------------------------------------------------------------------------------
#Diagramas de cajas (Boxplot) para todas las muestras
par(mfrow = c(1, 1))  # Restablece la ventana de gráficos a una sola gráfica
colors <- c(rep("orange", 6), rep("skyblue", 6)) #Colores para diferenciar entre grupos de muestras
boxplot(assays(se)$counts, col=colors, main="Valores de expresión de metablotios en las muestras (2 grupos)", ylab="Expresión", las=2, cex.axis=0.8)


## --------------------------------------------------------------------------------------------
#Diagramas de cajas (Boxplot) para todas las muestras con escala logarítmica
par(mfrow = c(1, 1))  # Restablece la ventana de gráficos a una sola gráfica
colors <- c(rep("coral", 6), rep("blue", 6)) #Colores para diferenciar entre grupos de muestras
boxplot(log(assays(se)$counts), col=colors, main="Valores de log de expresión de metablotios en las muestras (2 grupos)", ylab="Expresión", las=2, cex.axis=0.8)


## --------------------------------------------------------------------------------------------
#Cálculo de las componentes principales (PC)
pc <- prcomp(t(assays(se)$counts), scale. = TRUE)
summary(pc)


## --------------------------------------------------------------------------------------------
#Generar el gráfico de las dos primeras componentes principales
plot(pc$x[, 1], pc$x[, 2],
     #Para generar los nombres de los ejes, obtenemos las proporciones de varianza de las dos PCs
     xlab = paste("PC1 (", round(summary(pc)$importance[2, 1] * 100, 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(summary(pc)$importance[2, 2] * 100, 2), "%)", sep = ""),
     col = colors, lwd = 2,
     main = "Gráfico de las dos primeras componentes principales")

#Obtenemos los nombres de las columnas del dataset para añadirlos al gráfico de PCs
sample_names <- paste0(colnames(assays(se)$counts), 1:6)
text(pc$x[,1],pc$x[,2],sample_names, pos=3, cex=0.6)


## --------------------------------------------------------------------------------------------
#Generar el gráfico de las cinco primeras componentes principales
pca_scores <- as.data.frame(pc$x[, 1:5]) #Puntuaciones de las primeras cinco componentes

#Se genera el gráfico por pares
pairs(pca_scores, 
      main = "Gráfico de pares de las primeras cinco componentes crincipales",
      labels = paste("PC", 1:5,"(",round(summary(pc)$importance[2, 1:5] * 100, 2), "%)", sep=""),
      col=colors, lwd = 2)


## --------------------------------------------------------------------------------------------
#Dendograma de la agrupación jerárquica (cluster) de las muestras
#Cálculo de las distancias y composición del cluster
clust_dist <- hclust(dist(t(assays(se)$counts)), method = "complete")
plot(clust_dist, hang=-1)


## ----code, file="Jimenez_Guijarro_Beatriz_PEC1.R", eval=FALSE--------------------------------
## 

