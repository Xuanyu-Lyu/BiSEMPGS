# test the covariance of the original files and merged files

library(data.table)
testdf <- fread("Data/Paper/Model_r2_8/nfam32000/loop5.rds_32000.txt")
psych::describe(testdf)

x16k <- fread("Data/Paper/Model_latent30/nfam16000/loop3.rds_16000.txt")
x16k2 <- fread("Data/Paper/Model_latent30/nfam16000/loop56.rds_16000.txt")
x48k <- fread("Data/Paper/Model_latent30/nfam64000/loop1.rds_64000.txt")
x32k <- fread("Data/Paper/Model_latent30/nfam64000/loop1.rds_64000.txt")
cor(x16k)
cor(x48k)
cor(x32k)
x16k3 <- rbind(x16k, x16k2)

# visualize the covariance matrix
library(ggplot2)
library(reshape2)

# get the difference between 16k and 32k
diff <- cor(x16k) - cor(x16k3)
diff_var <- cov(x16k) - cov(x16k3)

melted_x16k <- melt(cor(x16k))
melted_x48k <- melt(cor(x48k))
melted_x32k <- melt(cor(x32k))
melted_x16k3 <- melt(cor(x16k3))

melted_diff <- melt(diff)
melted_diff_var <- melt(diff_var)
# create a heatmap with numbers in the cells
ggplot(melted_x16k, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 4)), vjust = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "blue", high = "red")

ggplot(melted_x16k3, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 4)), vjust = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = "blue", high = "red")


ggplot(melted_x48k, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 4)), vjust = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low = "blue", high = "red")

ggplot(melted_x32k, aes(Var1, Var2, fill = value)) +   
    geom_tile() +
    geom_text(aes(label = round(value, 4)), vjust = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low = "blue", high = "red")

print(ggplot(melted_diff, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 4)), vjust = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradient2(low = "blue", high = "red") )
