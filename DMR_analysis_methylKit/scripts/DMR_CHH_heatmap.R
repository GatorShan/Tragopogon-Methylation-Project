library(ggplot2)
df_CHH = read.table("/blue/soltis/shan158538/Methylation/OutPut/DMR_heatmap/CHH/Filtered_sorted_all_reformatted.txt", sep="\t", header = TRUE)

# Convert Data Frame Column to Vector by Subsetting Data; and reverse the order
Order <- df_CHH[1:25002,"chr_pos"]
Order <- rev(Order)

# Convert the y variable to a factor with specified order
df_CHH$chr_pos <- factor(df_CHH$chr_pos, levels = Order)

# Plot the heatmap with all settings
CHH_plot <- ggplot(df_CHH, aes(x = comparison, y = chr_pos, fill = diff)) +
  geom_tile() +
  scale_colour_gradient2(low = "deepskyblue3",
                         mid = "grey99",
                         high = "orangered3",
                         midpoint = 0,
                         space = "Lab",
                         breaks = c(70, 10, 0, -10, -70),
                         labels = c("70%", "10%", "", "-10%", "-70%"),
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

CHH_plot

