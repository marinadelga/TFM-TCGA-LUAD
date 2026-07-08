# ---------------------------------------------
# TFM - Analysis
# Latest date: 31/03/2026
# Author: Marina Delgado
# --------------------------------------------
# Data files:
# data_mrna_seq_tpm.txt        # Gene expression (normalizado por longitud de gen y profundidad de secuenciación)
# data_clinical_patient.txt    # Clinical data (está la supervivencia (OS_STATUS y OS_MONTHS))
# data_mutations.txt           # Vital para identificar quién tiene mutación en EGFR o KRAS
# -------------------------------------------

setwd(dir = "/home/mdelgadov1/TFM_Marina/luad_tcga_gdc")

library(org.Hs.eg.db)
library(dplyr)
library(survival)
library(tidyr)

# 1. Cargar datos
expression <- read.delim("data_mrna_seq_tpm.txt", check.names = FALSE)
clinical <- read.delim("data_clinical_patient.txt", skip = 4)

# 2. Preparar la matriz
## Calculamos la media de expresión de los genes repetidos
exp_colapsada <- aggregate(. ~ Entrez_Gene_Id, data = expression, FUN = mean)

## Ponemos el EntrezID como nombre de fila y quitamos la columna sobrante
rownames(exp_colapsada) <- exp_colapsada$Entrez_Gene_Id
exp_final_previa <- exp_colapsada[, -1]

## Limpiar nombres de pacientes (Recortar el -01A)
colnames(exp_final_previa) <- substr(colnames(exp_final_previa), 1, 12)

# 3. Sincronizar con la tabla Clínica
comunes <- intersect(colnames(exp_final_previa), clinical$PATIENT_ID)
expression_final <- exp_final_previa[, comunes]
clinical_final <- clinical[match(comunes, clinical$PATIENT_ID), ]

## Comprobación de seguridad: debe salir TRUE
print(paste("¿Sincronización total?:", all(colnames(expression_final) == clinical_final$PATIENT_ID)))

# 4. Limpiar Supervivencia 
clinical_final$OS_EVENT <- ifelse(clinical_final$OS_STATUS == "1:DECEASED", 1, 0)
clinical_final$OS_MONTHS <- as.numeric(clinical_final$OS_MONTHS)

## Quitamos pacientes con tiempo 0 o NA (darían error en la curva)
keep_surv <- !is.na(clinical_final$OS_MONTHS) & clinical_final$OS_MONTHS > 0
expression_final <- expression_final[, keep_surv]
clinical_final <- clinical_final[keep_surv, ]

## Comprobaciones: ¿Son iguales? ¿Se creó OS_EVENT? ¿OS_MONTHS es numérico?
print(all(colnames(expression_final) == clinical_final$PATIENT_ID))
head(clinical_final[, c("OS_STATUS", "OS_EVENT")], 10)
class(clinical_final$OS_MONTHS)

### Datos clínicos ###
# 1. Calcular Edad media
mean_age <- mean(as.numeric(clinical_final$AGE), na.rm = TRUE)
# 2. Sexo
tab_sex <- table(clinical_final$SEX)
# 3. Resumen de Estadios (simplificando nombres)
## Ajustar para que solo diga I, II, III, IV
clinical_final$STAGE_SIMPLE <- gsub("Stage ", "", clinical_final$PATH_STAGE)

clinical_final$STAGE_SIMPLE <- case_when(
  clinical_final$STAGE_SIMPLE %in% c("I", "IA", "IB") ~ "Estadio I",
  clinical_final$STAGE_SIMPLE %in% c("II", "IIA", "IIB") ~ "Estadio II",
  clinical_final$STAGE_SIMPLE %in% c("IIIA", "IIIB") ~ "Estadio III",
  clinical_final$STAGE_SIMPLE == "IV" ~ "Estadio IV",
  TRUE ~ "Sin clasificar"
)
table(clinical_final$STAGE_GROUP)

# 4. Tratamiento previo
## Saber si el paciente recibió tratamiento antes de la toma de muestra
table(clinical_final$PRIOR_TREATMENT)
# 5. Tabaquismo
## Crear una columna lógica de "Fumador" combinando variables
### Definimos fumador si años > 0 O si paquetes-año > 0
clinical_final <- clinical_final %>%
  mutate(
    es_fumador = case_when(
      (as.numeric(SMOKER_YEARS) > 0 | as.numeric(SMOKING_PACK_YEARS) > 0) ~ "Fumador",
      (as.numeric(SMOKER_YEARS) == 0 & as.numeric(SMOKING_PACK_YEARS) == 0) ~ "No fumador",
      TRUE ~ "No reportado" # Aquí caen los que tienen NA en ambas
    )
  )

# 6. Seguimiento
mean_os <- mean(as.numeric(clinical_final$OS_MONTHS), na.rm = TRUE)

## Mostrar resultados
cat("Edad media:", round(mean_age, 1), "\n")
cat("Sexo:\n"); print(round(prop.table(tab_sex)*100, 1))
print(table(clinical_final$SEX))
cat("Estadio:\n"); print(round(prop.table(tab_stage)*100, 1))
print(table(clinical_final$STAGE_SIMPLE))
cat("Seguimiento medio:", round(mean_os, 1), "\n")
table(clinical_final$es_fumador)
print(table(clinical_final$PRIOR_TREATMENT))

#### 1. EXPLORACIÓN DE MUTACIONES ####
# 1. Cargar el archivo de mutaciones
## Usamos comment.char = "#" para que ignore las líneas de metadatos
mut_data <- read.delim("data_mutations.txt", 
                       comment.char = "#", 
                       check.names = FALSE)

## Verificamos que se ha cargado bien
dim(mut_data)
head(mut_data[, 1:5])

# 2. Identificar los genes más frecuentes
library(dplyr)
library(ggplot2)

## Contar pacientes únicos por gen
frecuencia_muts <- mut_data %>%
  group_by(Hugo_Symbol) %>%
  summarise(n_pacientes = n_distinct(Tumor_Sample_Barcode)) %>%
  arrange(desc(n_pacientes))

## Ver el Top 10
print(head(frecuencia_muts, 10))

## Hacer una gráfica sencilla del Top 15
ggplot(head(frecuencia_muts, 15), aes(x = reorder(Hugo_Symbol, n_pacientes), y = n_pacientes)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Genes más mutados en conjuto TCGA-LUAD",
       x = "Gen", y = "Número de pacientes")

## Misma gráfica top 15 pero en procentaje
### Calculamos el total de pacientes en la cohorte final
total_pacientes <- nrow(clinical_final)
### Creamos el dataframe para el gráfico con la columna de porcentaje
datos_grafico <- head(frecuencia_muts, 15) %>%
  mutate(Porcentaje = (n_pacientes / total_pacientes) * 100)

### Gráfico en procentajes:
ggplot(datos_grafico, aes(x = reorder(Hugo_Symbol, Porcentaje), y = Porcentaje)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Genes más mutados en cohorte TCGA-LUAD",
       x = "Gen", 
       y = "Porcentaje de pacientes (%)") +
  # Añadir etiquetas con el porcentaje encima de las barras
  geom_text(aes(label = sprintf("%.1f%%", Porcentaje)), hjust = -0.1, size = 3)

# 3. Ver frecuencia de los genes driver (según literatura)
drivers_pulmon <- c("TP53", "KRAS", "EGFR", "STK11", "KEAP1", "NF1", "BRAF", "SETD2", "MET", "ALK")

tabla_justificacion <- frecuencia_muts %>%
  filter(Hugo_Symbol %in% drivers_pulmon) %>%
  mutate(Porcentaje = (n_pacientes / ncol(expression_final)) * 100)

print(tabla_justificacion)

## Identificar IDs de pacientes (12 caracteres)
egfr_ids <- mut_data %>% filter(Hugo_Symbol == "EGFR") %>% pull(Tumor_Sample_Barcode) %>% substr(1, 12) %>% unique()
kras_ids <- mut_data %>% filter(Hugo_Symbol == "KRAS") %>% pull(Tumor_Sample_Barcode) %>% substr(1, 12) %>% unique()
stk11_ids <- mut_data %>% filter(Hugo_Symbol == "STK11") %>% pull(Tumor_Sample_Barcode) %>% substr(1, 12) %>% unique()

## Clasificar en la tabla clínica
clinical_final <- clinical_final %>%
  mutate(Grupo_Estudio = case_when(
    PATIENT_ID %in% egfr_ids ~ "EGFR_mut",
    PATIENT_ID %in% kras_ids ~ "KRAS_mut",
    PATIENT_ID %in% stk11_ids ~ "STK11_mut", # Añadimos este por curiosidad
    TRUE ~ "Otros_WT"
  ))
## Ver resultado final grupos
table(clinical_final$Grupo_Estudio)

#### 2. CURVA DE KAPLAN-MEIER ####
library(survival)

# 1. Crear el objeto de supervivencia
## Usar el tiempo (OS_MONTHS) y el evento numérico (OS_EVENT)
fit_km <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ Grupo_Estudio, data = clinical_final)

# 2. Convertir el objeto fit_km a un data.frame que ggplot entienda
df_plot <- data.frame(
  time = fit_km$time,
  surv = fit_km$surv,
  upper = fit_km$upper,
  lower = fit_km$lower,
  group = rep(names(fit_km$strata), fit_km$strata)
)

## Quitar el "Grupo_Estudio=
df_plot$group <- gsub("Grupo_Estudio=", "", df_plot$group)

# 3. Dibujar la gráfica
ggplot(df_plot, aes(x = time, y = surv, color = group)) +
  geom_step(size = 1) + # geom_step es lo que crea el efecto "escalera"
  scale_color_manual(values = c("EGFR_mut" = "yellow", 
                                "KRAS_mut" = "blue", 
                                "Otros_WT" = "grey32", 
                                "STK11_mut" = "red")) +
  theme_minimal() +
  labs(title = "Curva de Supervivencia por Grupo Mutacional",
       subtitle = "Cohorte TCGA-LUAD",
       x = "Meses",
       y = "Probabilidad de Supervivencia",
       color = "Mutación") +
  ylim(0, 1) # La supervivencia va de 0 a 1

# 4. Buscar el "p-value" en la consola
surv_test <- survdiff(Surv(OS_MONTHS, OS_EVENT) ~ Grupo_Estudio, data = clinical_final)
print(surv_test)

### EXTRA: Diagrama de Venn
# 1. Obtener los IDs que sobrevivieron al filtro
ids_finales <- clinical_final$PATIENT_ID

# 2. Redefinir los grupos usando SOLO esos IDs
egfr_ids  <- mut_data %>% filter(Hugo_Symbol == "EGFR") %>% pull(Tumor_Sample_Barcode) %>% substr(1, 12) %>% unique() %>% intersect(ids_finales)
kras_ids  <- mut_data %>% filter(Hugo_Symbol == "KRAS") %>% pull(Tumor_Sample_Barcode) %>% substr(1, 12) %>% unique() %>% intersect(ids_finales)
stk11_ids <- mut_data %>% filter(Hugo_Symbol == "STK11") %>% pull(Tumor_Sample_Barcode) %>% substr(1, 12) %>% unique() %>% intersect(ids_finales)

library(VennDiagram)
# 3.Generar el diagrama
venn.plot <- venn.diagram(
  x = list(
    EGFR = egfr_ids,
    KRAS = kras_ids,
    STK11 = stk11_ids
  ),
  filename = NULL,
  fill = c("yellow", "blue", "red"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.1
)
grid.newpage()
grid.draw(venn.plot)

#### 3. ANÁLISIS DE EXPRESIÓN DIFERENCIAL EXPLORATORIO PARA SELECCION DE GENES CANDIDATOS ####
# Este análisis se realizó sobre datos log2(TPM + 1) y se utilizó únicamente como etapa exploratoria para seleccionar los genes
# candidatos con los que posteriormente se desarrolló la firma transcriptómica.
# El análisis de expresión diferencial definitivo incluido en el TFM se realizó de forma independiente utilizando conteos crudos
# (raw counts) y el método limma-voom (script: DEARaw.R), que constituye el análisis metodológicamente válido reportado en los
# resultados.

# 1. Filtrar los pacientes que nos interesan (EGFR y KRAS)
pacientes_interes <- which(clinical_final$Grupo_Estudio %in% c("EGFR_mut", "KRAS_mut"))
counts_subset <- expression_final[, pacientes_interes]

# 2. Crear el factor de grupo (importante el orden de los niveles)
## Poner KRAS primero para que sea el "control" y el resultado nos diga qué hace EGFR
grupo <- factor(clinical_final$Grupo_Estudio[pacientes_interes], 
                levels = c("KRAS_mut", "EGFR_mut"))

# 3. Preparar la anotación 
library(org.Hs.eg.db)
genes_info <- AnnotationDbi::select(org.Hs.eg.db, 
                     keys = rownames(expression_final),
                     columns = c("ENTREZID", "SYMBOL"), 
                     keytype = "ENTREZID")

## Limpiar duplicados
genes_info <- genes_info[!duplicated(genes_info$ENTREZID), ]
## Forzar a que 'genes_info' siga el MISMO orden que 'expression_final'
genes_info <- genes_info[match(rownames(expression_final), genes_info$ENTREZID), ]
## Comprobación final
all(rownames(expression_final) == genes_info$ENTREZID)

library(limma)

# 4. Log-transformación
log_tpm <- log2(expression_final + 1)

## Solo pacientes de EGFR y KRAS para la comparación
pacientes_dea <- which(clinical_final$Grupo_Estudio %in% c("EGFR_mut", "KRAS_mut"))
log_tpm_subset <- log_tpm[, pacientes_dea]
grupo_dea <- factor(clinical_final$Grupo_Estudio[pacientes_dea], levels = c("KRAS_mut", "EGFR_mut"))

# 5. Modelo Lineal
design <- model.matrix(~ grupo_dea)
fit <- lmFit(log_tpm_subset, design)
fit <- eBayes(fit, trend = TRUE)

## Tabla de resultados
res_dea <- topTable(fit, coef = 2, number = Inf)

## Añadir los Symbols 
res_dea$SYMBOL <- genes_info$SYMBOL[match(rownames(res_dea), genes_info$ENTREZID)]

## Mirar resultados
head(res_dea)

# 6.Volcano plot
library(ggplot2)

##Definir los colores y etiquetas
res_dea$significado <- "No Sig"
res_dea$significado[res_dea$logFC > 1 & res_dea$adj.P.Val < 0.05] <- "Up en EGFR"
res_dea$significado[res_dea$logFC < -1 & res_dea$adj.P.Val < 0.05] <- "Up en KRAS"

## Seleccionar los 15 mejores para evitar saturar el gráfico
top_15 <- head(res_dea, 15)

## Crear el gráfico
ggplot(res_dea, aes(x = logFC, y = -log10(adj.P.Val), color = significado)) +
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
  geom_text(data = top_15, 
            aes(label = SYMBOL), 
            vjust = -0.7,  # Desplaza el texto un poco hacia arriba del punto
            size = 3.5, 
            color = "black",
            fontface = "bold")

# 7. Tabla top 15 genes
## Extraemos exactamente los mismos 15 genes escritos en el gráfico
tabla_final <- head(res_dea[, c("SYMBOL", "logFC", "adj.P.Val")], 15)

## Añadir la etiqueta de color/grupo para que se entienda el Volcano
tabla_final$Categoria <- "Gris (Cambio leve)"
tabla_final$Categoria[tabla_final$logFC > 1 & tabla_final$adj.P.Val < 0.05] <- "Rojo (Up en EGFR)"
tabla_final$Categoria[tabla_final$logFC < -1 & tabla_final$adj.P.Val < 0.05] <- "Azul (Up en KRAS)"

## Recortar p-valores para que no salgan tantos decimales
tabla_final$adj.P.Val <- formatC(tabla_final$adj.P.Val, format = "e", digits = 2)
tabla_final$logFC <- round(tabla_final$logFC, 2)

## Mostrar la tabla
print(tabla_final)

# 8. Pheatmap 
library(pheatmap)

##. Seleccionar los 30 genes más significativos del DEA
top_30_ids <- rownames(res_dea)[1:30]
top_30_symbols <- res_dea$SYMBOL[1:30]

##. Extrar los datos de expresión (log_tpm) de esos genes y esos pacientes
matriz_hp <- log_tpm_subset[top_30_ids, ]
rownames(matriz_hp) <- top_30_symbols

##. Crear una barra de color para identificar los grupos
df_col <- data.frame(Mutacion = grupo_dea)
rownames(df_col) <- colnames(matriz_hp)

##. Dibujar el Heatmap
pheatmap(matriz_hp, 
         scale = "row",           # Normalizar por filas para ver contrastes
         annotation_col = df_col, # Añadir la barrita de EGFR vs KRAS
         show_colnames = FALSE,   # Quitar nombres de pacientes (son muchos)
         main = "Heatmap: Top 30 Genes Diferenciales",
         color = colorRampPalette(c("navy", "white", "red"))(100))

#### 4. GEN ONTOLOGY ANALYSIS ####
library(limma)

# 1. Seleccionar los genes significativos (FDR < 0.05)
## Separamos los que suben en EGFR de los que suben en KRAS
genes_up_egfr <- rownames(res_dea[res_dea$logFC > 0 & res_dea$adj.P.Val < 0.05, ])
genes_up_kras <- rownames(res_dea[res_dea$logFC < 0 & res_dea$adj.P.Val < 0.05, ])

# 2. Ejecutar el análisis de Gene Ontology (GO)
go_result <- goana(list(EGFR = genes_up_egfr, KRAS = genes_up_kras), 
                   species = "Hs")

# 3. Ver los procesos biológicos (BP) más importantes para cada grupo
top_go_egfr <- topGO(go_result, ontology = "BP", sort = "EGFR", number = 10)
top_go_kras <- topGO(go_result, ontology = "BP", sort = "KRAS", number = 10)

# 4. Mostrar resultados
print("--- TOP 10 PROCESOS BIOLÓGICOS (EGFR) ---")
print(top_go_egfr[, c("Term", "N", "EGFR", "P.EGFR")])

print("--- TOP 10 PROCESOS BIOLÓGICOS (KRAS) ---")
print(top_go_kras[, c("Term", "N", "KRAS", "P.KRAS")])

#### 5. RISK SCORE ####
# 1. Selección de genes a partir del top_30 DEA
library(survival)

ids_top30 <- rownames(res_dea)[1:30]

##Crear una función que testea la supervivencia gen por gen
# Usar log_tpm para que la escala sea estadísticamente correcta
lista_cribado <- lapply(ids_top30, function(id_gen) {
  simbolo <- res_dea$SYMBOL[rownames(res_dea) == id_gen]
  df_temp <- data.frame(
    Time  = clinical_final$OS_MONTHS,
    Event = clinical_final$OS_EVENT,
    Exp   = as.numeric(log_tpm[id_gen, clinical_final$PATIENT_ID])
  )
  df_temp <- na.omit(df_temp)
  # Ejecutar el Modelo de Cox Univariante
  fit_cox <- coxph(Surv(Time, Event) ~ Exp, data = df_temp)
  sum_cox <- summary(fit_cox)
  # Devolver los datos clave
  return(data.frame(
    ID = id_gen,
    Symbol = simbolo,
    P_val = sum_cox$logtest["pvalue"],      # Significancia estadística
    Coef = sum_cox$coefficients[1],         # Dirección del efecto (positivo o negativo)
    Hazard_Ratio = sum_cox$conf.int[1],     # >1 Riesgo, <1 Protector
    Lower_CI = sum_cox$conf.int[3],         # Intervalo de confianza inferior
    Upper_CI = sum_cox$conf.int[4]          # Intervalo de confianza superior
  ))
})

##. Unir todo en una tabla y filtrar por P-valor < 0.05
tabla_supervivencia_total <- do.call(rbind, lista_cribado)
genes_finales_firma <- tabla_supervivencia_total[tabla_supervivencia_total$P_val < 0.05, ]

## Ordenar por p-valor para ver los más potentes arriba
genes_finales_firma <- genes_finales_firma[order(genes_finales_firma$P_val), ]

## Mostrar el resultado
print("=== RESULTADOS DEL CRIBADO DE SUPERVIVENCIA (P < 0.05) ===")
print(genes_finales_firma)

# 2. Extraer los datos de los 4 genes ganadores
genes_ids <- c("115416", "5723", "55065", "100506810")
coefs <- c(0.4446, 0.2028, -0.1932, -0.2740) 
names(coefs) <- genes_ids

# 3. Calcular el Score sumando (Expresión * Coeficiente)
## Usamos log_tpm para que la escala sea la misma que en el test de Cox
matriz_firma <- t(log_tpm[genes_ids, clinical_final$PATIENT_ID])
clinical_final$Risk_Score <- as.matrix(matriz_firma) %*% coefs

# 4. Cox continuo  
clinical_final$Risk_Score_Z <- as.numeric(scale(clinical_final$Risk_Score))
fit_continuo <- coxph(
  Surv(OS_MONTHS, OS_EVENT) ~ Risk_Score_Z,
  data = clinical_final
)
summary(fit_continuo)

# 5. Curva de Kaplan-Meier
## Clasificar por la mediana
umbral <- median(clinical_final$Risk_Score, na.rm = TRUE)
clinical_final$Firma_Riesgo <- ifelse(clinical_final$Risk_Score > umbral, "Alto Riesgo", "Bajo Riesgo")

## Ver cuántos tenemos en cada grupo
table(clinical_final$Firma_Riesgo)

## Paquetes
library(survival)
library(ggplot2)
library(broom)
library(dplyr) # Para limpiar los nombres fácilmente

##. Ajustar el modelo
fit_km <- survfit(Surv(OS_MONTHS, OS_EVENT) ~ Firma_Riesgo, data = clinical_final)

##. Convertir a tabla y LIMPIAR los nombres de las columnas
df_tidy <- tidy(fit_km) %>%
  mutate(strata = gsub("Firma_Riesgo=", "", strata))

##. Calcular p-valor
sdiff <- survdiff(Surv(OS_MONTHS, OS_EVENT) ~ Firma_Riesgo, data = clinical_final)
p_val <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
p_label <- paste0("Log-rank p = ", format.pval(p_val, digits = 3))

## Gráfica
ggplot(df_tidy, aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_step(size = 1) +  
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  geom_point(data = filter(df_tidy, n.event == 0), shape = 3, size = 1.5) +
  scale_color_manual(values = c("Alto Riesgo" = "red", "Bajo Riesgo" = "blue")) +
  scale_fill_manual(values = c("Alto Riesgo" = "red", "Bajo Riesgo" = "blue")) +
  theme_minimal() +
  labs(title = "Validación de la Firma de 4 Genes",
       subtitle = p_label,
       x = "Meses de Seguimiento",
       y = "Probabilidad de Supervivencia") +
  ylim(0, 1)

# 6. Diferencia de riesgo entre grupos de estudio
## Calcular el modelo de Cox para los grupos
fit_cox_grupos <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ Firma_Riesgo, data = clinical_final)
summary(fit_cox_grupos)

library(ggplot2)

ggplot(clinical_final, aes(x = Grupo_Estudio, y = Risk_Score, fill = Grupo_Estudio)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) + # Esto dibuja los puntos de cada paciente
  theme_minimal() +
  labs(title = "Diferencia de Riesgo entre EGFR y KRAS",
       x = "Mutación (Grupo)",
       y = "Puntuación de Riesgo (Firma de 4 genes)")

# 7. ¿Qué predice mejor la muerte? ¿El Grupo_Estudio (mutaciones) o la Firma_Riesgo?
fit_comparativo <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ Grupo_Estudio + Risk_Score_Z, 
                         data = clinical_final)
summary(fit_comparativo)

# 8. Factores de confusión: Firma, Edad y Estadio
clinical_final$PATH_STAGE2 <- as.character(clinical_final$PATH_STAGE)

clinical_final$PATH_STAGE2[clinical_final$PATH_STAGE2 %in% c("Stage I","Stage IA","Stage IB")] <- "Stage I"
clinical_final$PATH_STAGE2[clinical_final$PATH_STAGE2 %in% c("Stage II","Stage IIA","Stage IIB")] <- "Stage II"
clinical_final$PATH_STAGE2[clinical_final$PATH_STAGE2 %in% c("Stage IIIA","Stage IIIB")] <- "Stage III"
table(clinical_final$PATH_STAGE2)
##Asegurarse de que PATH_STAGE2 sea un factor
clinical_final$PATH_STAGE2 <- as.factor(clinical_final$PATH_STAGE2)

fit_final <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ Risk_Score_Z + AGE + PATH_STAGE2, 
                   data = clinical_final)

summary(fit_final)
table(clinical_final$PATH_STAGE2)

library(forestmodel)
## Generar el Forest Plot 
forest_plot <- forest_model(fit_final)
print(forest_plot)

# 9. ¿Varía el riesgo con el Estadio? (Si el p-valor es > 0.05, no hay relación)
## Elimina los niveles que ya no tienen datos
data_limpio <- clinical_final[
  !is.na(clinical_final$PATH_STAGE2) & 
    clinical_final$PATH_STAGE2 != "" & 
    clinical_final$PATH_STAGE2 != " ", 
]

## Convertir a factor 
data_limpio$PATH_STAGE2 <- as.factor(data_limpio$PATH_STAGE2)
## Krustal test
resultado_kruskal <- kruskal.test(Risk_Score ~ PATH_STAGE2, data = data_limpio)
p_val_text <- sprintf("%.6f", resultado_kruskal$p.value)                         
## Graficar
ggplot(data_limpio, aes(x = PATH_STAGE2, y = Risk_Score, fill = PATH_STAGE2)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribución del Risk Score por Estadio",
    subtitle = paste("Kruskal-Wallis, p-value =", p_val_text),
    x = "Estadio Patológico", 
    y = "Puntuación de Riesgo"
  ) +
  theme(legend.position = "none")

# 10. ¿Varía el riesgo con la Edad? (Si el p-valor es > 0.05, no hay relación)
# Filtrar para el análisis de correlación
data_cor <- clinical_final[!is.na(clinical_final$AGE) & !is.na(clinical_final$Risk_Score), ]
# Ahora hacer correlación y gráfico solo con los que queden limpios
cor.test(data_cor$Risk_Score, data_cor$AGE)

ggplot(data_cor, aes(x = AGE, y = Risk_Score)) +
  geom_point(alpha = 0.4, color = "darkblue") +  # Puntos de cada paciente
  geom_smooth(method = "lm", color = "red") +    # Línea de tendencia (Regresión)
  theme_minimal() +
  labs(title = "Relación entre Edad y Risk Score",
       subtitle = "Correlación de Pearson: -0.14 (p = 0.0011)",
       x = "Edad (años)", y = "Puntuación de Riesgo")

#### 6. MARCADORES MICROAMBIENTE TUMORAL (TME) ####
# 1. Buscar los IDs de los marcadores del microambiente
genes_tme_info <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = c("CD274", "ACTA2", "CD163"), 
                                        columns = c("ENTREZID", "SYMBOL"), 
                                        keytype = "SYMBOL")

# 2. Extraer la expresión de log_tpm (usando los ENTREZID)
## clinical_final$PATIENT_ID para asegurar que el orden sea el mismo
tme_data <- as.data.frame(t(log_tpm[genes_tme_info$ENTREZID, clinical_final$PATIENT_ID]))
colnames(tme_data) <- genes_tme_info$SYMBOL

# 3. Unir estos nuevos datos a tu tabla clinical_final
clinical_final <- cbind(clinical_final, tme_data)

# A. Relación con Fibroblastos (ACTA2)
ggplot(clinical_final, aes(x = Firma_Riesgo, y = ACTA2, fill = Firma_Riesgo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("Alto Riesgo" = "red", "Bajo Riesgo" = "blue")) +
  labs(title = "Presencia de Fibroblastos (ACTA2) según Firma de Riesgo",
       subtitle = "TCGA-LUAD",
       y = "Expresión log2(TPM)", x = "Grupo de Riesgo")

## Test estadístico para ACTA2
wilcox.test(ACTA2 ~ Firma_Riesgo, data = clinical_final)

##  Test estadístico para CD274 (PD-L1)
wilcox.test(CD274 ~ Firma_Riesgo, data = clinical_final)

# B. Relación con Evasión Inmune (CD274 / PD-L1)
ggplot(clinical_final, aes(x = Firma_Riesgo, y = CD274, fill = Firma_Riesgo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_minimal() +
  scale_fill_manual(values = c("Alto Riesgo" = "red", "Bajo Riesgo" = "blue")) +
  labs(title = "Expresión de PD-L1 (CD274) según Firma de Riesgo",
       subtitle = "Marcador de Evasión Inmunológica",
       y = "Expresión log2(TPM)", x = "Grupo de Riesgo")

## Correlación con PD-L1
cor.test(clinical_final$Risk_Score, clinical_final$CD274, method = "spearman")

## Gráfico Risk Score vs PD-L1 (CD274)
ggplot(clinical_final, aes(x = Risk_Score, y = CD274)) +
  geom_point(alpha = 0.5, color = "blue", size = 2) + # Puntos con transparencia
  geom_smooth(method = "lm", color = "red", fill = "grey30") + # Línea de tendencia
  theme_minimal() +
  labs(
    title = "Correlación entre Risk Score y Evasión Inmune",
    subtitle = paste0("Spearman Rho = 0.119 | p-value = 0.0077"),
    x = "Puntuación de Riesgo (Firma 4 genes)",
    y = "Expresión de PD-L1 (CD274) [log2 TPM]"
  ) +
  # Añadir el texto estadístico dentro del gráfico
  annotate("text", x = max(clinical_final$Risk_Score)*0.7, 
           y = min(clinical_final$CD274), 
           label = "p < 0.01", color = "red", fontface = "bold")

## Correlación con ACTA2
cor.test(clinical_final$Risk_Score, clinical_final$ACTA2, method = "spearman")

#### 7. COMPROBACIONES EXTRA ####
# 1. Correlación interna de la firma
## Extraer los valores de expresión de los 4 genes para todos los pacientes
genes_ids <- c("115416", "5723", "55065", "100506810")
exp_firma <- as.data.frame(t(log_tpm[genes_ids, clinical_final$PATIENT_ID]))

## Cambiar los nombres de las columnas a los símbolos para que se entienda mejor
## 115416=MALSU1, 5723=PSPH, 55065=SLC52A1, 100506810=LINC01132
colnames(exp_firma) <- c("MALSU1", "PSPH", "SLC52A1", "LINC01132")

##. Calcular la correlación entre el Risk_Score y los genes
cor_firma <- cor(cbind(Risk_Score = clinical_final$Risk_Score, exp_firma), method="spearman")
print(cor_firma)

# 2. Curva ROC
library(survivalROC)
## Calcular el AUC a los 3 años (36 meses)
roc_3years <- survivalROC(Stime = clinical_final$OS_MONTHS,  
                          status = clinical_final$OS_EVENT,      
                          marker = clinical_final$Risk_Score,     
                          predict.time = 36, 
                          method = "KM")

## Ver el valor del AUC
print(paste("AUC a los 3 años:", round(roc_3years$AUC, 3)))

## Dibujar la curva
plot(roc_3years$FP, roc_3years$TP, type = "l", col = "red", lwd = 3,
     xlab = "1 - Especificidad", ylab = "Sensibilidad",
     main = "Curva ROC de la Firma (3 años)")
abline(0, 1, lty = 2, col = "grey")
text(0.6, 0.2, paste("AUC =", round(roc_3years$AUC, 3)), cex = 1.2, font = 2)

# 3. Modelo que sume el Estadio y el Risk_Score
modelo_combinado <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ Risk_Score + PATH_STAGE, data = clinical_final)

##. Calcular el C-index del modelo combinado
cindex_solo_firma <- summary(coxph(Surv(OS_MONTHS, OS_EVENT) ~ Risk_Score, data = clinical_final))$concordance[1]
cindex_combinado <- summary(modelo_combinado)$concordance[1]

print(paste("C-index de la Firma sola:", round(cindex_solo_firma, 3)))
print(paste("C-index Combinado (Firma + Estadio):", round(cindex_combinado, 3)))

# 4. Análisis de Inmunidad Avanzada (Linfocitos T y Macrófagos)
cd8a_exp <- as.data.frame(t(log_tpm["925", clinical_final$PATIENT_ID]))
clinical_final$CD8A <- cd8a_exp[,1]

## Correlación con Linfocitos T (CD8A)
cor_cd8 <- cor.test(clinical_final$Risk_Score, clinical_final$CD8A, method = "spearman")
## Correlación con Macrófagos M2 (CD163)
cor_m2 <- cor.test(clinical_final$Risk_Score, clinical_final$CD163, method = "spearman")

print(paste("Rho CD8A (Linfocitos):", round(cor_cd8$estimate, 3), "p-val:", round(cor_cd8$p.value, 4)))
print(paste("Rho CD163 (Macrófagos):", round(cor_m2$estimate, 3), "p-val:", round(cor_m2$p.value, 4)))

#### 8. VALIDACIÓN CRUZADA (10-fold Cross-Validation)####
library(survival)

set.seed(123) # Para que siempre salgan los mismos números
n_folds <- 10
# 1. Crear los grupos de forma aleatoria
indices <- sample(1:nrow(clinical_final))
folds <- cut(1:nrow(clinical_final), breaks = n_folds, labels = FALSE)

cindex_results <- c()

for(i in 1:n_folds){
  test_indices <- indices[which(folds == i)]
  train_data <- clinical_final[-test_indices, ]
  test_data  <- clinical_final[test_indices, ]
  
  modelo_cv <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ Risk_Score, data = train_data)
  
  c_val <- concordance(modelo_cv, newdata = test_data)$concordance
  cindex_results <- c(cindex_results, c_val)
}

cat("Resultado de la Validación Cruzada (C-index):\n")
print(cindex_results)
cat("Media del C-index:", mean(cindex_results), "\n")
cat("Desviación estándar:", sd(cindex_results), "\n")

# 2. Gráfico para ver la estabilidad
boxplot(cindex_results, col = "purple", main = "Estabilidad del C-index (10-fold CV)",
        ylab = "Concordance Index")

#### 9. DECISION CURVE ANALYSIS (DCA) ####
# 1. Definir los Modelos que vamos a comparar
clinical_final$PATH_STAGE <- as.factor(clinical_final$PATH_STAGE2)
## Solo el Estadio Patológico
m1 <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ PATH_STAGE, data = clinical_final)
##Firma de 4 genes + Estadio
m2 <- coxph(Surv(OS_MONTHS, OS_EVENT) ~ Risk_Score + PATH_STAGE, data = clinical_final)

# 2. Función para calcular el Beneficio Neto (Net Benefit)
calcular_nb <- function(modelo, datos, umbral) {
  ## Predecir la probabilidad de fallecer a los 3 años (36 meses)
  pred_surv <- summary(survfit(modelo, newdata = datos), times = 36)$surv
  riesgo <- 1 - pred_surv
  
  n <- nrow(datos)
  ## ¿Quien ha muerto antes de los 36 meses?
  evento_36 <- ifelse(datos$OS_MONTHS <= 36 & datos$OS_EVENT == 1, 1, 0)
  
  tp <- sum(riesgo >= umbral & evento_36 == 1)
  fp <- sum(riesgo >= umbral & evento_36 == 0)
  
  ## Fórmula del Beneficio Neto
  nb <- (tp / n) - (fp / n) * (umbral / (1 - umbral))
  return(nb)
}

# 3. Calcular los resultados para un rango de umbrales (del 5% al 50%)
thresholds <- seq(0.05, 0.50, by = 0.01)

nb_solo_estadio <- sapply(thresholds, function(u) calcular_nb(m1, clinical_final, u))
nb_firma_estadio <- sapply(thresholds, function(u) calcular_nb(m2, clinical_final, u))


# 4. Gráfico
library(ggplot2)
library(tidyr)

df_dca <- data.frame(
  Umbral = thresholds,
  Solo_Estadio = nb_solo_estadio,
  Firma_Mas_Estadio = nb_firma_estadio
)

# Pasamos a formato largo
df_dca_long <- pivot_longer(df_dca, cols = -Umbral, names_to = "Modelo", values_to = "NB")
# Graficar
ggplot(df_dca_long, aes(x = Umbral, y = NB, color = Modelo)) +
  geom_line(size = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(
    values = c("Firma_Mas_Estadio" = "blue", "Solo_Estadio" = "red"),
    labels = c("Firma de 4 genes + Estadio", "Solo Estadio (TNM)")
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Utilidad Clínica: Decision Curve Analysis (DCA)",
    subtitle = "Beneficio neto de añadir la firma genética al diagnóstico estándar",
    x = "Umbral de Riesgo (Probabilidad de Decisión)",
    y = "Beneficio Neto (Net Benefit)",
    color = "Estrategia"
  ) +
  theme(legend.position = "top", plot.title = element_text(face = "bold")) +
  annotate("text", x = 0.4, y = max(nb_firma_estadio)*0.7, 
           label = "Superioridad de la Firma", color = "blue", fontface = "bold.italic")
