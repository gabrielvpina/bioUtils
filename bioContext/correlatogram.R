library(ggplot2)
library(dplyr)
library(viridis)
library(Cairo)

df <- read.csv("result.blastn.tsv", header = T, sep = "\t")

# Pacotes necessários
library(ggplot2)
library(dplyr)


df <- df %>%
  mutate(identity_group = ifelse(pident >= 80, "≥ 80%", "< 80%"))

# Plot
imagem <- ggplot(df, aes(x = qseqid, y = sseqid, size = qcovs, color = pident, fill=identity_group)) +
  geom_point(alpha = 0.8, shape = 21) +
  scale_size(range = c(1, 9)) +
  #scale_color_viridis_c(option = "D", end = 0.9) +
  scale_fill_brewer(palette = "Set1", name = "≥ 80% Identity") +
  labs(
    x = "EVEs",
    y = "EVEs",
    size = "Coverage (%)",
    color = "Identity (%)",
    fill = "Upper than 80%"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"))

CairoPDF(paste0("Corrplot_Gemini", ".pdf"), width = 8.06, height = 6.76, family = "arial")
imagem
dev.off()



imagem <- ggplot(df, aes(x = V1, y = V2, size = V6, fill = V8, color = identity_group)) +
  geom_point(alpha = 0.9, shape = 21, stroke = 0.4) +
  scale_size(range = c(1, 9), name = "Coverage (%)") +
  scale_fill_viridis_c(name = "Identity (%)", option = "D", end = 0.9) +
  scale_color_manual(
    name = "≥ 80% Identity",
    values = c("≥ 80%" = "red", "< 80%" = "white")
  ) +
  labs(x = "Viruses", y = "Viruses") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(colour = "black")
  )
imagem
