######################## Richness and Shannon ##################

##Load libraries
library(vegan)
library(scales)
library(RColorBrewer)
library(VennDiagram)
library(gplots)
library(car)
library(pairwiseAdonis)
library(ggVennDiagram)
library(ggplot2)
library(eulerr)
library(devtools)
library(forcats)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


################# Colapse by Genus
tab_coi <- read.xlsx("Supplementary Data S2.1_COI.xlsx")
tab_18 <- read.xlsx("Supplementary Data S2.2_18S.xlsx")
metadata<-read.xlsx("metadata.xlsx")

tab_coi_genus <- tab_coi %>%
  group_by(genus) %>%
  summarise(across(where(is.numeric), sum))

tab_coi_genus <- tab_coi_genus[c(1:134),]


tab_coi_genus2 <- tab_coi_genus[,c(2:265),]

rownames(tab_coi_genus2)<-tab_coi_genus$genus


# pooledtab2 <- pooledtab2[rowSums(pooledtab2)>1,] ## Fora tots els OTUs amb menys de 2 reads


## Alfa Diversities COI

diver<-diversity(tab_coi_genus2,"shannon",MARGIN=2)

pooledtabt<- as.data.frame(t(tab_coi_genus2))

richsep<-specnumber(pooledtabt)
igu<- diver/log(specnumber(pooledtabt))
simp<-diversity(pooledtabt, index = "simpson")

datos<-data.frame(richsep,diver,simp,igu)
datos$Sample <- rownames(datos) 
datos$SampleID <- str_replace(datos$Sample, "\\.\\d+\\.$", "")
datos <- left_join(datos, metadata, by = "SampleID")

rownames(datos)<-c(datos$Sample)

datos <- datos[,c(1:4,9,7,10)]
colnames(datos) <- c("Richness","Shannon","Simpson", "Igualdad","Age","Core", "Env")


datos_coi <- datos
richness_coi<-datos_coi
richness_coi$Richness<-richness_coi$Richness/10


richness_coi<- richness_coi[,c(1,2,5,7)]
richness_coi$Sample <- rownames(richness_coi)

names(richness_coi)<-c("Richness","Shannon", "Year", "Env", "Sample")


richness_coi <- melt(richness_coi, id= c("Year","Sample", "Env"))
richness_coi$Year<-as.character(richness_coi$Year)



ggplot(richness_coi, mapping=aes(x=Year, y=value, fill=variable)) +
  geom_boxplot()+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, size = 7.5))+
  xlab(NULL)+
  ylab("Shannon diversity - Species richness/100 COI ALL ASVs REPs")+
  facet_grid(factor(.~Env, levels=c("River","Estuary")), scales="free_x")+
  scale_fill_manual(name ="variable", breaks=c("Shannon","Richness"), values=c("orange","#815c97"))+
  facet_wrap(.~Env, ncol=1, scales = 'free_y')+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal() +
  theme( legend.position="bottom", 
         # panel.grid.minor = element_blank(), 
         # panel.grid.major = element_blank(),
         axis.text.x = element_text(angle = 90, size = 9),
         panel.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_discrete()



medians_df <- richness_coi %>%
  group_by(Year, variable, Env) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = 'drop')

ggplot(richness_coi, aes(x = Year, y = value, fill = variable)) +
  geom_boxplot()+
  #geom_point(aes(fill = variable), color = "black", shape = 21, size = 1)+
  geom_line(data = medians_df, aes(x = Year, y = median_value, group = variable, color = variable), alpha=0.5,size = 1) +
  scale_fill_manual(name = "variable", breaks = c("Shannon", "Richness"),
                    values = c("Shannon" = "orange", "Richness" = "#815c97")) +
  scale_color_manual(values = c("Shannon" = "orange", "Richness" = "#815c97")) +
  theme_bw() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, size = 9),
    panel.background = element_rect(fill = "transparent", colour = NA)
  ) +
  xlab(NULL) +
  ylab("Shannon diversity - Species richness/10 COI ASVs") +
  facet_wrap(. ~ Env, ncol=1, scales = 'free_y') +
  scale_x_discrete()




############ 18S ###############

tab_18_genus <- tab_18 %>%
  group_by(genus) %>%
  summarise(across(where(is.numeric), sum))

tab_18_genus <- tab_18_genus[c(1:609),]


tab_18_genus2 <- tab_18_genus[,c(3:265)]

rownames(tab_18_genus2)<-tab_18_genus$genus


# pooledtab2 <- pooledtab2[rowSums(pooledtab2)>1,] ## Fora tots els OTUs amb menys de 2 reads


## Alfa Diversities 18S

diver<-diversity(tab_18_genus2,"shannon",MARGIN=2)

pooledtabt<- as.data.frame(t(tab_18_genus2))

richsep<-specnumber(pooledtabt)
igu<- diver/log(specnumber(pooledtabt))
simp<-diversity(pooledtabt, index = "simpson")

datos<-data.frame(richsep,diver,simp,igu)
datos$Sample <- rownames(datos) 
datos$SampleID <- str_replace(datos$Sample, "\\.\\d+\\.$", "")
datos <- left_join(datos, metadata, by = "SampleID")

rownames(datos)<-c(datos$Sample)

datos <- datos[,c(1:4,9,7,10)]
colnames(datos) <- c("Richness","Shannon","Simpson", "Igualdad","Age","Core", "Env")


datos_18 <- datos
richness_18<-datos_18
richness_18$Richness<-richness_18$Richness/10


richness_18<- richness_18[,c(1,2,5,7)]
richness_18$Sample <- rownames(richness_18)

names(richness_18)<-c("Richness","Shannon", "Year", "Env", "Sample")


richness_18 <- melt(richness_18, id= c("Year","Sample", "Env"))
richness_18$Year<-as.character(richness_18$Year)



ggplot(richness_18, mapping=aes(x=Year, y=value, fill=variable)) +
  geom_boxplot()+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, size = 7.5))+
  xlab(NULL)+
  ylab("Shannon diversity - Species richness/100 COI ALL ASVs REPs")+
  facet_grid(factor(.~Env, levels=c("River","Estuary")))+
  scale_fill_manual(name ="variable", breaks=c("Shannon","Richness"), values=c("orange","#815c97"))+
  facet_grid(.~Env, scales = 'free_x')+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  theme_minimal() +
  theme( legend.position="bottom", 
         # panel.grid.minor = element_blank(), 
         # panel.grid.major = element_blank(),
         axis.text.x = element_text(angle = 90, size = 9),
         panel.background = element_rect(fill = "transparent", colour = NA)) +
  scale_x_discrete()



medians_df <- richness_18 %>%
  group_by(Year, variable, Env) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = 'drop')

ggplot(richness_18, aes(x = Year, y = value, fill = variable)) +
  geom_boxplot()+
  #geom_point(aes(fill = variable), color = "black", shape = 21, size = 1)+
  geom_line(data = medians_df, aes(x = Year, y = median_value, group = variable, color = variable), alpha=0.5,size = 1) +
  scale_fill_manual(name = "variable", breaks = c("Shannon", "Richness"),
                    values = c("Shannon" = "orange", "Richness" = "#815c97")) +
  scale_color_manual(values = c("Shannon" = "orange", "Richness" = "#815c97")) +
  theme_bw() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, size = 9),
    panel.background = element_rect(fill = "transparent", colour = NA)
  ) +
  xlab(NULL) +
  ylab("Shannon diversity - Species richness/10 18S ASVs") +
  facet_grid(. ~ Env, scales = 'free_x') +
  scale_x_discrete()





##### Merge Both Genes Datasets

richness_coi<-datos_coi
richness_coi$Richness<-richness_coi$Richness


richness_coi<- richness_coi[,c(1,2,5,7)]
richness_coi$Sample <- rownames(richness_coi)

names(richness_coi)<-c("Richness","Shannon", "Year", "Env", "Sample")


richness_coi <- melt(richness_coi, id= c("Year","Sample", "Env"))
richness_coi$Year<-as.character(richness_coi$Year)



richness_18<-datos_18
richness_18$Richness<-richness_18$Richness


richness_18<- richness_18[,c(1,2,5,7)]
richness_18$Sample <- rownames(richness_18)

names(richness_18)<-c("Richness","Shannon", "Year", "Env", "Sample")


richness_18 <- melt(richness_18, id= c("Year","Sample", "Env"))
richness_18$Year<-as.character(richness_18$Year)



richness_coi$Gene <- "COI"
richness_18$Gene <- "18S"

richness_all <- rbind(richness_coi, richness_18)

richness_all <- richness_all %>%
  filter(variable=="Richness")

medians_df <- richness_all %>%
  group_by(Year, Env, Gene) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = 'drop')

ggplot(richness_all, aes(x = Year, y = value, fill = Gene)) +
  geom_boxplot()+
  #geom_point(aes(fill = variable), color = "black", shape = 21, size = 1)+
  geom_line(data = medians_df, aes(x = Year, y = median_value, group = Gene, color = Gene), alpha=0.5,size = 1) +
  scale_fill_manual(name = "variable", breaks = c("COI", "18S"),
                    values = c("COI" = "orange", "18S" = "#815c97")) +
  scale_color_manual(values = c("COI" = "orange", "18S" = "#815c97")) +
  theme_bw() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, size = 9),
    panel.background = element_rect(fill = "transparent", colour = NA)
  ) +
  xlab(NULL) +
  ylab("Richness Collapsed by Genus") +
  facet_grid(. ~ Env) +
  scale_x_discrete()

richness_all$Year <- as.numeric(as.character(richness_all$Year))
medians_df$Year <- as.numeric(as.character(medians_df$Year))

ggplot(richness_all, aes(x = as.character(Year), y = value, fill = Gene)) +
  geom_boxplot() +
  geom_line(data = medians_df, aes(x = as.character(Year), y = median_value, group = Gene, color = Gene), alpha = 0.5, size = 1) +
  scale_fill_manual(name = "variable", breaks = c("COI", "18S"),
                    values = c("COI" = "orange", "18S" = "#815c97")) +
  scale_color_manual(values = c("COI" = "orange", "18S" = "#815c97")) +
  theme_bw() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, size = 6),
    panel.background = element_rect(fill = "transparent", colour = NA)
  ) +
  xlab("Year-Gene") +
  ylab("Richness Collapsed by Genus") +
  facet_wrap(~Env, ncol = 1)

ggplot(richness_all, aes(x = Year, y = value, color=Gene, fill = Gene)) +
  geom_boxplot(aes(group = interaction(Year, Gene)),
               position = position_dodge2(width = 0.9, preserve = "single"),
               width = 5) +
  geom_smooth(method = "gam", se = TRUE, formula = y ~ s(x, k = 4), linewidth = 0.5, alpha = 0.3) +
  scale_x_continuous(breaks = pretty(richness_all$Year, n = 30)) +
  scale_y_continuous(breaks = pretty(richness_all$value, n = 10))+
  facet_wrap(~ Env, scales = "free_y", ncol=1) +
  labs(title = "Genus Richness ASV", x = "Age (cal BP)", y = "Richness", color = "Gene") +
  theme_minimal()+
  scale_color_manual(values = c("COI" = "orange", "18S" = "#815c97")) +
  scale_fill_manual(values = c("COI" = "orange", "18S" = "#815c97")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    strip.text = element_text(face = "bold", size=18),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    title = element_text(size=20)
  )




#################### MDS #####################


all_data <- read.xlsx("Supplementary Data S3 Relative Abundance Both Genes.xlsx")

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggforce)  # para círculos con borde

# Paso 1: Pivotear a formato ancho por ASV y sample
asv_matrix <- all_data %>%
  pivot_longer(cols = starts_with("H"), names_to = "SampleID", values_to = "Abundance") %>%
  mutate(SampleGeneID = paste(SampleID, gene, sep = "_")) %>%
  group_by(SampleGeneID, ASV) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = ASV, values_from = Abundance, values_fill = 0)

# 2. Metadatos únicos para cada SampleGeneID
metadata_nmds <- asv_matrix %>%
  select(SampleGeneID) %>%
  left_join(all_data %>%
              pivot_longer(cols = starts_with("H"), names_to = "SampleID", values_to = "Abundance") %>%
              mutate(SampleGeneID = paste(SampleID, gene, sep = "_")) %>%
              select(SampleID, gene, SampleGeneID) %>%
              distinct() %>%
              left_join(metadata, by = "SampleID"),
            by = "SampleGeneID")

# 3. Convertir a matriz y correr NMDS
comm_matrix <- asv_matrix %>%
  column_to_rownames("SampleGeneID") %>%
  as.data.frame()

# Correr NMDS
nmds <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)
##stress
nmds_stress <- nmds$stress
print(nmds_stress)
# 4. Extraer coordenadas
nmds_coords <- as.data.frame(as.matrix(scores(nmds)))
nmds_coords$SampleGeneID <- rownames(nmds_coords)

nmds_plot_data <- left_join(nmds_coords, metadata_nmds, by = "SampleGeneID")

nmds_points <- as.data.frame(nmds$points)
nmds_points$SampleGeneID <- rownames(nmds_points)
nmds_points <- left_join(nmds_points, metadata_nmds, by = "SampleGeneID")


# Obtener edades únicas ordenadas
edades <- sort(unique(nmds_points$Age))

# Crear paleta de gradiente de azul a amarillo
gradiente_azul_amarillo <- colorRampPalette(c("#313695", "yellow"))(length(edades))

# Asignar colores a cada edad
edad_colors <- setNames(gradiente_azul_amarillo, edades)

ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = factor(Age), group=gene)) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = edad_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "nMDS ASVs by Age",
    color = "Age"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggplot(nmds_points, aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = factor(Age), shape = gene), 
    color = "black",     # Borde negro para todos (puedes mapearlo si quieres que cambie por gene)
    size = 4,
    stroke = 1.2,         # Grosor del borde
    alpha = 0.9
  ) +
  scale_fill_manual(values = edad_colors, name = "Age") +
  scale_shape_manual(values = c("COI" = 21, "18S" = 24)) + 
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 21, color = "black", stroke = 1.2))
  ) +
  facet_wrap(~Env) +
  theme_minimal(base_size = 14) +
  labs(
    title = "nMDS ASVs by Age and Gene",
    shape = "Gene"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    #legend.text = element_text(size = 10)
  )



ggplot(nmds_points %>% filter(Env == "Hythe"), aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = factor(Age), shape = gene),
    color = "black",
    size = 4,
    stroke = 1.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = edad_colors, name = "Age") +
  scale_shape_manual(values = c("COI" = 21, "18S" = 24)) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 21, color = "black", stroke = 1.2))
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "nMDS – Hythe ",
    shape = "Gene"
  ) +
  facet_grid(.~gene, scales="free_x")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )


ggplot(nmds_points %>% filter(Env == "Bursledon"), aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = factor(Age), shape = gene),
    color = "black",
    size = 4,
    stroke = 1.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = edad_colors, name = "Age") +
  scale_shape_manual(values = c("COI" = 21, "18S" = 24)) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 21, color = "black", stroke = 1.2))
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "nMDS – Bunny Meadows",
    shape = "Gene"
  ) +
  facet_grid(.~gene, scales="free_x")+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10)
  )

stress_label <- paste0("Stress = ", round(nmds$stress, 3))

ggplot(nmds_points, aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = factor(Age), shape = gene), 
    color = "black",
    size = 4,
    stroke = 1.2,
    alpha = 0.9
  ) +
  scale_fill_manual(values = edad_colors, name = "Age") +
  scale_shape_manual(values = c("COI" = 21, "18S" = 24)) + 
  facet_wrap(~Env) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("nMDS ASVs by Age and Gene (", stress_label, ")"),
    shape = "Gene"
  )
