### XRF Data #########

### Load XRF data normalised
xrf_h<-read.xlsx("Supplementary Data S1_XRF_Dating.xlsx")

xrf_long_h <- xrf_h %>%
  pivot_longer(cols = c(Pb, Zn, Cu, Ca, "Ca/10"), names_to = "Element", values_to = "Concentration")

ggplot(xrf_long_h, aes(x = Year, y = Concentration, color = Element)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = pretty(xrf_long_h$Year, n = 30)) +
  facet_wrap(~ Env, ncol=1, scales="free_y") +
  labs(
    title = "XRF - Elementos Traza",
    x = "Age (cal BP)",
    y = "Concentración (ppm)",
    color = "Elemento"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold")
  )
