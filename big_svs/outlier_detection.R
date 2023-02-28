library(ggplot2)
library(tidyr)
library(dplyr)
library(nlraa)

args = commandArgs(trailingOnly=TRUE)
data <-read.csv(args[1], sep = ",", header = FALSE)

# matrix distribution taking into account position
data_df <- as.data.frame(data)

#get row numbers and save them as a column
data_df$nrow <- seq.int(nrow(data_df))

#get column numbers and save them as a column
data_df <- pivot_longer(data_df, !nrow, names_to = "ncol")

#omin nas
data_df <- na.omit(data_df)

#modify ncols to make them integers
data_df$ncol <- as.integer(gsub("V", "", data_df$ncol))

#take only upper triangle values; i>j; ncol>nrow
data_df_upper <- data_df[data_df$ncol>data_df$nrow,]

data_df_upper$nrow <- as.numeric(data_df_upper$nrow)
data_df_upper$ncol <- as.numeric(data_df_upper$ncol)

# calculate the index by the distance to the matrix
data_df_upper <- group_by(data_df_upper, nrow) %>%
  mutate(index=ncol-nrow)

points <- data_df_upper %>% ungroup() %>% dplyr::select(value, index)
colnames(points) <- c("y", "x")

# fit the model
fo <- y ~ exp(a + b * exp(x*-c))
fm <- nls(fo, points, start = list(a = 1, b = 1, c=1))

# calculate and plot the prediction bands
fm.Theoph.prd.bnd <- predict2_nls(fm, interval = "prediction", level = 0.95)
fm.Theoph.prd.bnd.dat <- cbind(points, fm.Theoph.prd.bnd)

## Plot it
plot = fm.Theoph.prd.bnd.dat %>%
  ggplot(aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(x = x, y = Estimate), colour="blue") +
  geom_ribbon(data = fm.Theoph.prd.bnd.dat,
              aes(x = x, ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha=0.3) +
  xlab("Distance from matrix") + ylab("Similarity index") +
  ggtitle("95% prediction bands")+
  #geom_vline(xintercept = 5)+
  geom_point(data = fm.Theoph.prd.bnd.dat %>%
               mutate(out=(y>Q97.5 | y<Q2.5)) %>%
               filter(out==TRUE), aes(x = x, y = y), colour="#FF6600")+
  theme(legend.position = "none")+
  theme_minimal()

##outlier list
outliers = fm.Theoph.prd.bnd.dat %>% cbind(data_df_upper) %>%
  mutate(upper=(y>Q97.5)) %>%
  mutate(lower=(y<Q2.5)) %>%
  filter(upper==TRUE|lower==TRUE,x>5) %>%
  dplyr::select("nrow", "ncol", "value", "Estimate", "Est.Error", "Q2.5", "Q97.5", "upper", "lower")


# write out
write.table(outliers, file=paste(args[2], ".csv", sep=""), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

png(file=paste(args[2], "_plot.png", sep=""),
width=600, height=350)
plot
dev.off()