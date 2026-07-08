# TFM - New Analysis with raw counts 
# Latest date: 25/06/2026
# Author: Marina Delgado Valero
# Purpose:
# - Differential expression analysis reported in my Master's thesis.
# - This script corresponds to the methodological validation of the exploratory analysis performed in ProyectoTFM.R.
# --------------------------------------------
# Data files:
# data_mrna_seq_read_counts.txt # Conteos crudos
# data_mrna_seq_tpm.txt        # Gene expression (normalizado por longitud de gen y profundidad de secuenciación)
# data_clinical_patient.txt    # Clinical data (está la supervivencia (OS_STATUS y OS_MONTHS))
# data_mutations.txt           # Vital para identificar quién tiene mutación en EGFR o KRAS

setwd(dir = "/home/mdelgadov1/TFM_Marina/luad_tcga_gdc")

library(org.Hs.eg.db)
library(dplyr)
library(survival)
library(tidyr)

# 1. Cargar datos crudos
raw_counts <- read.delim("data_mrna_seq_read_counts.txt", check.names = FALSE)
clinical <- read.delim("data_clinical_patient.txt", skip = 4)

# 2. colapsar genes duplicados (si los hay)
raw_counts_colapsado <- aggregate(. ~ Entrez_Gene_Id, data = raw_counts, FUN = sum)
rownames(raw_counts_colapsado) <- raw_counts_colapsado$Entrez_Gene_Id
raw_counts_final <- raw_counts_colapsado[, -1]

# 3. Limpiar nombres de pacientes 
colnames(raw_counts_final) <- substr(colnames(raw_counts_final), 1, 12)

# 4. Sincronizar con clinical_final
raw_counts_subset <- raw_counts_final[, clinical_final$PATIENT_ID]
## Comprobación
print(paste("¿Coinciden los pacientes de los conteos con los de la clínica?:", 
            all(colnames(raw_counts_subset) == clinical_final$PATIENT_ID)))

# 5. Asegurar que son enteros (limma-voom lo requiere)
raw_counts_subset <- round(raw_counts_subset)

library(limma)
library(edgeR)

#### 1. DIFFERENTIAL EXPRESSION ANALYSIS (FINAL DEA) ####
## Este analisis compara tumores EGFR-mutado y KRAS-mutado usando datos de conteo crudos RNA-seq y el flujo limma-voom.
# 1. Definir los pacientes del DEA (EGFR vs KRAS)
pacientes_dea_new <- which(clinical_final$Grupo_Estudio %in% c("EGFR_mut", "KRAS_mut"))
# Usar raw_counts_subset (que ya limpiamos para que tenga los mismos pacientes)
counts_dea <- raw_counts_subset[, pacientes_dea_new]
grupo_dea_new <- factor(clinical_final$Grupo_Estudio[pacientes_dea_new], levels = c("KRAS_mut", "EGFR_mut"))

# 2. Preparar la anotación desde cero (usando el objeto org.Hs.eg.db)
library(org.Hs.eg.db)

## Seleccionamos información
genes_info_new <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = rownames(raw_counts_subset),
                                    columns = c("SYMBOL", "GENENAME"), 
                                    keytype = "ENTREZID")

# 3. Limpieza de duplicados 
genes_info_new <- genes_info_new[!duplicated(genes_info_new$ENTREZID), ]

## Forzamos a que 'genes_info_new' siga el MISMO orden que 'raw_counts_subset'
genes_info_new <- genes_info_new[match(rownames(raw_counts_subset), genes_info_new$ENTREZID), ]
## Comprobación final
all(rownames(raw_counts_subset) == genes_info_new$ENTREZID)

# 4. Crear objeto DGEList
dge <- DGEList(counts = counts_dea)

# 5. Filtrado de baja expresión (Para quitar genes que no aportan nada)
design_new <- model.matrix(~grupo_dea_new)
keep <- filterByExpr(dge, design_new)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# 6. Normalización (Calcular factores que corrigen diferencias de secuenciación)
dge <- calcNormFactors(dge)

# 7. VOOM (Transformación necesaria para usar limma con conteos)
v <- voom(dge, design_new, plot=TRUE) 

# 8. Ajuste del modelo lineal
fit_new <- lmFit(v, design_new)
fit_new <- eBayes(fit_new)

# 9. Obtención de resultados (Top Table)
res_dea_new <- topTable(fit_new, coef = 2, number = Inf)

# 10. Integrar la anotación que ya había preparada
# Usar 'match' con los rownames de res_dea para alinear los símbolos
res_dea_new$SYMBOL <- genes_info_new$SYMBOL[match(rownames(res_dea_new), genes_info_new$ENTREZID)]
res_dea_new$GENENAME <- genes_info_new$GENENAME[match(rownames(res_dea_new), genes_info_new$ENTREZID)]

# 11. Ordenar por significancia estadística 
res_dea_new <- res_dea_new[order(res_dea_new$adj.P.Val), ]

## Mira los 10 genes que más diferencian entre EGFR y KRAS
head(res_dea_new, 10)

# 12.Volcano plot
library(ggplot2)

## Definir los colores y etiquetas
res_dea_new$significado <- "No Sig"
res_dea_new$significado[res_dea_new$logFC > 1 & res_dea_new$adj.P.Val < 0.05] <- "Up en EGFR"
res_dea_new$significado[res_dea_new$logFC < -1 & res_dea_new$adj.P.Val < 0.05] <- "Up en KRAS"

## Seleccionar los 15 mejores para evitar saturar el gráfico
top_15_new <- head(res_dea_new, 15)

##Crear el gráfico
ggplot(res_dea_new, aes(x = logFC, y = -log10(adj.P.Val), color = significado)) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("Up en EGFR" = "red", 
                                "Up en KRAS" = "blue", 
                                "No Sig" = "grey")) +
  theme_minimal() +
  labs(title = "Expresión Diferencial: EGFR mut vs KRAS mut",
       subtitle = "TCGA-LUAD | Genes con FDR < 0.05",
       x = "log2 Fold Change", y = "-log10 FDR") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Usamos geom_text (que viene en ggplot2 por defecto)
  geom_text(data = top_15_new, 
            aes(label = SYMBOL), 
            vjust = -0.7,  # Desplaza el texto un poco hacia arriba del punto
            size = 3.5, 
            color = "black",
            fontface = "bold")

# 13. Tabla top 15 genes
## Extraemos exactamente los mismos 15 genes escritos en el gráfico
tabla_final_new <- head(res_dea_new[, c("SYMBOL", "logFC", "adj.P.Val")], 15)

## Añadimos la etiqueta de color/grupo para que se entienda el Volcano
tabla_final_new$Categoria <- "Gris (Cambio leve)"
tabla_final_new$Categoria[tabla_final_new$logFC > 1 & tabla_final_new$adj.P.Val < 0.05] <- "Rojo (Up en EGFR)"
tabla_final_new$Categoria[tabla_final_new$logFC < -1 & tabla_final_new$adj.P.Val < 0.05] <- "Azul (Up en KRAS)"

## Recortamos p-valores para que no salgan tantos decimales
tabla_final_new$adj.P.Val <- formatC(tabla_final_new$adj.P.Val, format = "e", digits = 2)
tabla_final_new$logFC <- round(tabla_final_new$logFC, 2)

## Mostrar la tabla
print(tabla_final_new)

# 14. Pheatmap final
library(pheatmap)

## Top 30 genes más significativos
top_30_ids_new <- rownames(res_dea_new)[1:30]
top_30_symbols_new <- res_dea_new$SYMBOL[1:30]

## Comprobar que sale 0
sum(duplicated(res_dea_new$SYMBOL[1:30]))

## Matriz voom
matriz_hp_new <- v$E[top_30_ids_new, ]

## nombres de genes
rownames(matriz_hp_new) <- res_dea_new$SYMBOL[1:30]

## anotación de columnas
annotation_col_new <- data.frame(
  Mutacion = grupo_dea_new
)

rownames(annotation_col_new) <- colnames(matriz_hp_new)

## heatmap
pheatmap(
  matriz_hp_new,
  scale = "row",
  annotation_col_new = annotation_col_new,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 30 genes diferencialmente expresados",
  color = colorRampPalette(c("navy","white","firebrick3"))(100)
)

#### 2. GEN ONTOLOGY ANALYSIS ####
library(limma)

# 1. Seleccionar los genes significativos (FDR < 0.05)
## Separar los que suben en EGFR de los que suben en KRAS
genes_up_egfr_new <- rownames(res_dea_new[res_dea_new$logFC > 0 & res_dea_new$adj.P.Val < 0.05, ])
genes_up_kras_new <- rownames(res_dea_new[res_dea_new$logFC < 0 & res_dea_new$adj.P.Val < 0.05, ])

# 2. Ejecutar el análisis de Gene Ontology (GO)
go_result_new <- goana(list(EGFR = genes_up_egfr_new, KRAS = genes_up_kras_new), 
                   species = "Hs")

# 3. Ver los procesos biológicos (BP) más importantes para cada grupo
top_go_egfr_new <- topGO(go_result_new, ontology = "BP", sort = "EGFR", number = 10)
top_go_kras_new <- topGO(go_result_new, ontology = "BP", sort = "KRAS", number = 10)

# 4. Mostrar resultados
print("--- TOP 10 PROCESOS BIOLÓGICOS (EGFR) ---")
print(top_go_egfr_new[, c("Term", "N", "EGFR", "P.EGFR")])

print("--- TOP 10 PROCESOS BIOLÓGICOS (KRAS) ---")
print(top_go_kras_new[, c("Term", "N", "KRAS", "P.KRAS")])

#### PRINCIPAL COMPONENT ANALYSIS (PCA) ####
# PCA is used as an unsupervised exploratory analysis to evaluate whether 
# global gene expression profiles separate EGFR-mutated and KRAS-mutated tumours

#A) PCA total
# 1. PCA rápido con la función incorporada de limma
pca_res_new <- plotMDS(v, top = 500, plot = FALSE) 
# Usamos los 500 genes que más varían

# 2. Crear un dataframe para ggplot
df_pca_new <- data.frame(PC1 = pca_res_new$x, PC2 = pca_res_new$y, Grupo = grupo_dea_new)

# 3. Graficar
library(ggplot2)
ggplot(df_pca_new, aes(x = PC1, y = PC2, color = Grupo)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: Separación entre EGFR y KRAS",
       subtitle = "Basado en los 500 genes con mayor varianza")

# --- NUEVO PCA: Basado solo en genes diferencialmente expresados (FDR < 0.05) ---

# 1. Identificar genes con FDR < 0.05 (significativos)
genes_significativos_new <- rownames(res_dea_new[res_dea_new$adj.P.Val < 0.05, ])

# 2. Extraer los datos de expresión solo para esos genes del objeto 'v' (voom)
# 'v$E' contiene los datos transformados log2-CPM
matriz_significativos_new <- v$E[genes_significativos_new, ]

# 3. Realizar el PCA sobre esta matriz
# Usamos t(matriz) porque plotMDS espera genes en filas y muestras en columnas
pca_sig_new <- plotMDS(matriz_significativos_new, plot = FALSE)

# 4. Crear el dataframe para ggplot
df_pca_sig_new <- data.frame(PC1 = pca_sig_new$x, PC2 = pca_sig_new$y, Grupo = grupo_dea_new)

# 5. Graficar
ggplot(df_pca_sig_new, aes(x = PC1, y = PC2, color = Grupo)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA: Separación basada en genes significativos",
       subtitle = paste("Genes con FDR < 0.05 (n =", length(genes_significativos_new), ")"),
       x = "PC1",
       y = "PC2") +
  scale_color_manual(values = c("KRAS_mut" = "#F8766D", "EGFR_mut" = "#00BFC4"))


#--------JUSTIFICACIÓN -------------
res_dea_new[res_dea_new$SYMBOL %in%
          c("MALSU1","PSPH","SLC52A1","LINC01132"),
        c("SYMBOL","logFC","adj.P.Val")]
