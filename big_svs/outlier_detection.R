library(ggplot2)
library(tidyr)
library(dplyr)
library(nlraa)

# Outlier detection by Wrath combines two approaches: z-scores and modelling of the distribution of barcode sharing by distance from the diagonal
# Outliers are defined as values that fall outside the z-score (absolute) threshold AND outside the prediction bands of the model

zscore_threshold = 10
predition_level = 0.95

#### PART 1: DATA CLEANING ####
#read in the data
args = commandArgs(trailingOnly=TRUE)
data <-read.csv(args[1], sep = ",", header = FALSE)

# assign NAs to anything below and including the diagonal
data_m <- as.matrix(data)
data_m[lower.tri(data_m, diag=TRUE)] <- NA

# matrix distribution taking into account position
data_df <- as.data.frame(data_m)

#get row numbers and save them as a column
data_df$nrow <- seq.int(nrow(data_df))

#get column numbers and save them as a column
data_df <- pivot_longer(data_df, !nrow, names_to = "ncol")

#omit NAs
data_df <- na.omit(data_df)

#modify ncols to make them integers
data_df$ncol <- as.integer(gsub("V", "", data_df$ncol))

# calculate the index by the distance to the matrix
data_df_upper <- group_by(data_df, nrow) %>%
  mutate(index=ncol-nrow)

# create dataframe with x and y values to calculate z-scores from
points <- data_df_upper %>% ungroup() %>% dplyr::select(value, index, ncol, nrow)
colnames(points) <- c("y", "x", "ncol", "nrow")

#### PART 2: Z-SCORE CALCULATION ####
# calculate z scores grouping values by their distance to the diagonal
scaled_full_df <- tibble(z_score=0,y=0, x=0, nrow=0,ncol=0)

for (index in unique(points$x)) {
  group_indices <- which(points$x == index)
  group_values <- points[group_indices,]$y
  scaled_group <- scale(na.omit(group_values))
  scaled_df <- tibble(z_score=scaled_group,y=points[group_indices,]$y, x=index, nrow=points[group_indices,]$nrow, ncol=points[group_indices,]$ncol)
  scaled_full_df <- rbind(scaled_full_df,scaled_df)
}
scaled_full_df <- scaled_full_df[-1,]


#### PART 3: MODEL FIT AND PREDICTION BAND CALCULATION ####
# fit the model
fo <- y ~ exp(a + b * exp(-x*c))
fm <- nls(fo, points, start = list(a = 1, b = 1, c=1))

# calculate and plot the prediction bands
fm.Theoph.prd.bnd <- predict2_nls(fm, interval = "prediction", level = predition_level)

# change name of quantiles so they are not relative to the prediction thresfold 
# then prediction level can be changed and code doesn't break
names_quantiles <- names(fm.Theoph.prd.bnd)[3:4] 
names(fm.Theoph.prd.bnd)[3:4] <- c("Qbottom", "Qtop")
fm.Theoph.prd.bnd.dat <- cbind(points, fm.Theoph.prd.bnd)

# merge dataset containing z scores and model estimates
dataset = left_join(scaled_full_df, fm.Theoph.prd.bnd.dat)
dataset = dataset %>%  mutate(abs_z=abs(z_score))


#### PART 4: PLOT AND OUTLIER LIST OUTPUT ####
## Plot it
plot = dataset %>%
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(x = x, y = Estimate), colour="blue") +
  geom_ribbon(data = dataset,
              aes(x = x, ymin = Qbottom, ymax = Qtop), fill = "purple", alpha=0.3) +
  xlab("Distance from matrix") + ylab("Similarity index") +
  ggtitle(paste(predition_level*100, "% prediction bands", sep = ""))+
  geom_point(data = dataset %>%
               mutate(out=(y>Qtop | y<Qbottom)) %>%
               filter(out==TRUE, abs_z>zscore_threshold), aes(x = x, y = y), colour="#FF6600")+
  theme(legend.position = "none")+
  theme_minimal()

## save the plot
ggsave(paste(args[2], "_plot.png", sep=""), plot, width=6, height=3.5, dpi=300)

##outlier list
# outliers are defines as vales outside the threshold quantile and z-score
outliers = dataset %>% left_join(data_df_upper, by=(c("y"="value", "nrow", "ncol", "x"="index"))) %>%
  mutate(upper=(y>Qtop)) %>%
  mutate(lower=(y<Qbottom)) %>%
  filter(upper==TRUE|lower==TRUE, abs_z>zscore_threshold) %>%
  dplyr::select("nrow", "ncol", "y", "Estimate", "Est.Error", "Qbottom", "Qtop", "upper", "lower", "z_score")

names(outliers) <- c("nrow", "ncol", "value", "Estimate", "Est.Error", names_quantiles[1], names_quantiles[2], "upper", "lower", "z_score")

# write out the outliers and the plot
write.table(outliers, file=paste(args[2], ".csv", sep=""), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)