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

richness_coi <- richness_coi %>%
  mutate(
    Year = as.numeric(Year),
    Env = trimws(Env),
    Env = factor(Env, levels = c("Hythe", "Bursledon"))
  )

medians_df <- medians_df %>%
  mutate(
    Year = as.numeric(Year),
    Env = trimws(Env),
    Env = factor(Env, levels = c("Hythe", "Bursledon"))
  )

ggplot(richness_coi, aes(x = Year, y = value, fill = variable)) +
  geom_boxplot(
    aes(group = interaction(Year, variable)),
    width = 8,
    position = position_dodge(width = 8)
  ) +
  geom_line(
    data = medians_df,
    aes(x = Year, y = median_value, group = variable, color = variable),
    alpha = 0.5,
    linewidth = 1
  ) +
  scale_fill_manual(
    name = "variable",
    breaks = c("Shannon", "Richness"),
    values = c("Shannon" = "orange", "Richness" = "#815c97")
  ) +
  scale_color_manual(values = c("Shannon" = "orange", "Richness" = "#815c97")) +
  facet_wrap(~ Env, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(
      floor(min(richness_coi$Year, na.rm = TRUE) / 10) * 10,
      ceiling(max(richness_coi$Year, na.rm = TRUE) / 10) * 10,
      by = 10
    )
  ) +
  labs(
    title = "COI",
    x = "Age (cal BP)",
    y = "Shannon diversity - Species richness/10 ASVs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = "Shannon&Richness_diversity_COI.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)


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


richness_18 <- richness_18 %>%
  mutate(
    Year = as.numeric(Year),
    Env = trimws(Env),
    Env = factor(Env, levels = c("Hythe", "Bursledon"))
  )
medians_df <- richness_18 %>%
  group_by(Year, variable, Env) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = 'drop')

medians_df <- medians_df %>%
  mutate(
    Year = as.numeric(Year),
    Env = trimws(Env),
    Env = factor(Env, levels = c("Hythe", "Bursledon"))
  )

ggplot(richness_18, aes(x = Year, y = value, fill = variable)) +
  geom_boxplot(
    aes(group = interaction(Year, variable)),
    width = 8,
    position = position_dodge(width = 8)
  ) +
  geom_line(
    data = medians_df,
    aes(x = Year, y = median_value, group = variable, color = variable),
    alpha = 0.5,
    linewidth = 1
  ) +
  scale_fill_manual(
    name = "variable",
    breaks = c("Shannon", "Richness"),
    values = c("Shannon" = "orange", "Richness" = "#815c97")
  ) +
  scale_color_manual(values = c("Shannon" = "orange", "Richness" = "#815c97")) +
  facet_wrap(~ Env, ncol = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(
      floor(min(richness_coi$Year, na.rm = TRUE) / 10) * 10,
      ceiling(max(richness_coi$Year, na.rm = TRUE) / 10) * 10,
      by = 10
    )
  ) +
  labs(
    title = "18S",
    x = "Year",
    y = "Shannon diversity - Species richness/100 ASVs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 24),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = "Shannon&Richness_diversity_18S.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)





##### Merge Both Genes Datasets - SHANNON

richness_coi <- datos_coi
richness_coi <- richness_coi[, c(1, 2, 5, 7)]
richness_coi$Sample <- rownames(richness_coi)

names(richness_coi) <- c("Richness", "Shannon", "Year", "Env", "Sample")

richness_coi <- melt(richness_coi, id = c("Year", "Sample", "Env"))
richness_coi$Year <- as.character(richness_coi$Year)


richness_18 <- datos_18
richness_18 <- richness_18[, c(1, 2, 5, 7)]
richness_18$Sample <- rownames(richness_18)

names(richness_18) <- c("Richness", "Shannon", "Year", "Env", "Sample")

richness_18 <- melt(richness_18, id = c("Year", "Sample", "Env"))
richness_18$Year <- as.character(richness_18$Year)


richness_coi$Gene <- "COI"
richness_18$Gene <- "18S"

diversity_all <- rbind(richness_coi, richness_18)

# Seleccionar Shannon en vez de Richness
shannon_all <- diversity_all %>%
  filter(variable == "Shannon")

shannon_all$Year <- as.numeric(as.character(shannon_all$Year))

medians_df <- shannon_all %>%
  group_by(Year, Env, Gene) %>%
  summarize(median_value = median(value, na.rm = TRUE), .groups = "drop")

shannon_all$Env <- factor(
  shannon_all$Env,
  levels = c("Hythe", "Bursledon")
)

ggplot(shannon_all, aes(x = Year, y = value, color = Gene, fill = Gene)) +
  geom_boxplot(
    aes(group = interaction(Year, Gene)),
    position = position_dodge2(width = 0.9, preserve = "single"),
    width = 5
  ) +
  geom_smooth(
    method = "gam",
    se = TRUE,
    formula = y ~ s(x, k = 4),
    linewidth = 0.5,
    alpha = 0.3
  ) +
  scale_x_continuous(breaks = pretty(shannon_all$Year, n = 30)) +
  scale_y_continuous(breaks = pretty(shannon_all$value, n = 10)) +
  facet_wrap(~ Env, scales = "free_y", ncol = 1) +
  labs(
    title = "Genus Shannon Diversity ASV",
    x = "Age (cal BP)",
    y = "Shannon Index",
    color = "Gene",
    fill = "Gene"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("COI" = "orange", "18S" = "#815c97")) +
  scale_fill_manual(values = c("COI" = "orange", "18S" = "#815c97")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 18),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 20)
  )

library(svglite)

ggsave(
  filename = "Shannon_diversity.svg",
  width = 16,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)
