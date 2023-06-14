library(ggplot2)
df_CHG = read.table("/blue/soltis/shan158538/Methylation/OutPut/DMR_heatmap/CHG/Filtered_sorted_all_reformatted.txt", sep="\t", header = TRUE)

# Convert Data Frame Column to Vector by Subsetting Data; and reverse the order
Order <- df_CHG[1:16716,"chr_pos"]
Order <- rev(Order)

# Convert the y variable to a factor with specified order
df_CHG$chr_pos <- factor(df_CHG$chr_pos, levels = Order)

# Plot the heatmap with all settings
CHG_plot <- ggplot(df_CHG, aes(x = comparison, y = chr_pos, fill = diff)) +
  geom_tile() +
  scale_colour_gradient2(low = "deepskyblue3",
                         mid = "grey99",
                         high = "orangered3",
                         midpoint = 0,
                         space = "Lab",
                         breaks = c(95, 10, 0, -10, -95),
                         labels = c("95%", "10%", "", "-10%", "-95%"),
                         name = "Difference",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "fill") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 10),
        plot.margin = margin(2, 115, 2, 115, "mm"),
        text = element_text(family = "Arial")) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))


CHG_plot

