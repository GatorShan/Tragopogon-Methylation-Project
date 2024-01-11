library(tidyverse)

## CG CONTEXT
df_CG = read.table("/Users/shengchen/OneDrive/Methylation.Trag/Important.output/DMR_parental_legacy_quantitative_analysis/CG_parental_legacy_abs_diff.txt", header=TRUE)

# A density plot is a smoothed, continuous version of a histogram estimated from the data.
# The most common form of estimation is known as kernel density estimation.
# The kernel most often used is a Gaussian (which produces a Gaussian bell curve at each data point).
# The y-axis in a density plot is the probability density function for the kernel density estimation.
# However, we need to be careful to specify this is a probability density and not a probability.
# by default, the values of from and to are cut bandwidths beyond the extremes of the data
# the left and right-most points of the grid at which the density is to be estimated; the defaults are cut * bw outside of range(x).
dens_CG <- density(df_CG$Abs_diff, cut=0)

# creat a dataframe that have three columns: x, y, and additional one.
data_CG <- tibble(x = dens_CG$x, y = dens_CG$y) %>% 
  mutate(variable = case_when(
    (x > 35) ~ "Off",
    TRUE ~ NA_character_))

# plot the CG density plot
ggplot(data_CG, aes(x, y)) + geom_line() +
  geom_area(data = filter(data_CG, variable == 'Off'), fill = 'light blue') +
  geom_vline(xintercept = mean(df_CG$Abs_diff)) +
  geom_vline(xintercept = median(df_CG$Abs_diff), color = "red") +
  scale_x_continuous(limits = c(0, 65), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.065), expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) 

## CHG CONTEXT
df_CHG = read.table("/Users/shengchen/OneDrive/Methylation.Trag/Important.output/DMR_parental_legacy_quantitative_analysis/CHG_parental_legacy_abs_diff.txt", header=TRUE)
dens_CHG <- density(df_CHG$Abs_diff, cut=0)
data_CHG <- tibble(x = dens_CHG$x, y = dens_CHG$y) %>% 
  mutate(variable = case_when(
    (x > 25) ~ "Off",
    TRUE ~ NA_character_))
ggplot(data_CHG, aes(x, y)) + geom_line() +
  geom_area(data = filter(data_CHG, variable == 'Off'), fill = 'light blue') +
  geom_vline(xintercept = mean(df_CHG$Abs_diff)) +
  geom_vline(xintercept = median(df_CHG$Abs_diff), color = "red") +
  scale_x_continuous(limits = c(0, 70), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.07), expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) 

## CHH CONTEXT
df_CHH = read.table("/Users/shengchen/OneDrive/Methylation.Trag/Important.output/DMR_parental_legacy_quantitative_analysis/CHH_parental_legacy_abs_diff.txt", header=TRUE)
dens_CHH <- density(df_CHH$Abs_diff, cut=0)
data_CHH <- tibble(x = dens_CHH$x, y = dens_CHH$y) %>% 
  mutate(variable = case_when(
    (x > 10) ~ "Off",
    TRUE ~ NA_character_))
ggplot(data_CHH, aes(x, y)) + geom_line() +
  geom_area(data = filter(data_CHH, variable == 'Off'), fill = 'light blue') +
  geom_vline(xintercept = mean(df_CHH$Abs_diff)) +
  geom_vline(xintercept = median(df_CHH$Abs_diff), color = "red") +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.18), expand = c(0, 0)) +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) 

