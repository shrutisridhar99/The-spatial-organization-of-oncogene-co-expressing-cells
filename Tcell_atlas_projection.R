# Use AddModuleScore instead
tonsil <- AddModuleScore(tonsil, 
                         features = list(CD68_GC_up_top100$X),
                         name = "Dispersed_Score")

tonsil <- AddModuleScore(tonsil,
                         features = list(CD68_GC_down_top100$X),
                         name = "Clustered_Score")

score_summary <- tonsil@meta.data %>%
  mutate(Cluster = Idents(tonsil)) %>%
  group_by(Cluster) %>%
  summarise(
    Dispersed = mean(Dispersed_Score1),
    Clustered = mean(Clustered_Score1)
  ) %>%
  mutate(
    Dispersed_Z = scale(Dispersed)[,1],
    Clustered_Z = scale(Clustered)[,1],
    Diff = Dispersed_Z - Clustered_Z
  )

ggplot(score_summary, aes(x = "Dispersed vs Clustered", y = Cluster,
                          colour = Diff, size = abs(Diff))) +
  geom_point() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                         name = "Dispersed\nvs Clustered") +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  theme_classic() +
  labs(title = "Dispersed vs Clustered enrichment by T cell cluster",
       x = "", y = "Identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tonsil@meta.data$Dispersed_Z <- scale(tonsil@meta.data$Dispersed_Score1)[,1]
tonsil@meta.data$Clustered_Z <- scale(tonsil@meta.data$Clustered_Score1)[,1]
tonsil@meta.data$Diff_Score <- tonsil@meta.data$Dispersed_Z - tonsil@meta.data$Clustered_Z

FeaturePlot(tonsil, features = "Diff_Score") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")),
                         limits = c(-2, 2),
                         oob = scales::squish) +
  labs(title = "Dispersed vs Clustered Score")
