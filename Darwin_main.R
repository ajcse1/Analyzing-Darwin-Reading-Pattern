# #-----------------------------------------------------------------------
#   @Author - Abhishek Anand Singh
#   aasingh@indiana.edu
#   Code for Project Darwin
# #-----------------------------------------------------------------------

library(entropy)
library(ggplot2)
library(cluster)
library(proxy)
library(fpc)
library(lsa)
library(kknn)

#------------------PCA-----------------

#Applying PCA
Data_PC <- read.csv(file.choose(), header = TRUE) #Select the topic40.csv file
head(Data_PC)
Data_PC[,1] <- NULL #Remove the column containing book names
log_data <- log(Data_PC) #Applying log to make the data continous
head(log_data)
# pc <- prcomp(log_data, center = TRUE, scale. = TRUE,)
pc <- prcomp(Data_PC, center = TRUE, scale. = TRUE,)
plot(pc, type = 'l', main = "PRINCIPAL COMPONENT on T40" )

## ----------------CMD Scale -------------------


Y60 <- read.csv(file.choose(), header = T) # take topic 60 dataset
Y60[,1] <- NULL # remove the first column
YY <- dist(Y60 , method = "euclidean")  # Converting to distance matrix
y_cmd <- cmdscale(YY,k = 4, eig = TRUE) # Applying cmdscale function retrieve the CMDS output
points <- y_cmd$points 
x <- points[,1]
y <- -points[,2]
plot(x,y, xlab = "Dimension 1", ylab = "Dimension 2", main = "CMD Scale on T60")

# Using KL Divergence as distance for CMDS

Y60_KL <- read.csv(file.choose(), header = T) # take topic 60 KL dataset

Y60_KL_data <- Y60_KL[,2:594]
y_cmd_kl <- cmdscale(Y60_KL_data,k = 2, eig = TRUE)
points_KL <- y_cmd_kl$points
x <- points_KL[,1]
y <- -points_KL[,2]
plot(x,y, xlab = "Dimension 1", ylab = "Dimension 2", main = "CMD Scale on T80 KL")
text(x,y, Book_data$Color, cex=0.3, pos=4)



# -----------------ISOMAP-------------

# Gettting configuration using ISOMAP 
Y40_Data <- read.csv(file.choose(), header = T)
head(Y40_Data)
Y40_Data[,1]<-NULL
Y40_Dist <- mds.edm1(Y40_Data)  # Getting distance matrix
Y40_KNN <- graph.knn(Y40_Dist, k= 10) # Matrix KNN for 10 nearest neighbours
head(Y40_KNN)
Y40_adj <- graph.adj(KNN = Y40_KNN) # Calculating adjacency matrix
head(Y40_adj)
Y40_dis <- graph.dis(KNN = Y40_KNN, D = Y40_Dist) # Assigning dissimilarity edge weights based on adjacency
head(Y40_dis)
Y40_short <- graph.short(W = Y40_dis) # Computing all shortest path distances for a weighted graph.
head(Y40_short)
Y40_cmdscale <- cmdscale(Y40_short, k = 2) # Performing CMDS
head(Y40_cmdscale)
x <- Y40_cmdscale[,1]
y <- -Y40_cmdscale[,2]
plot(x,y, xlab = "Dimension 1", ylab = "Dimension 2", main = "ISOMAP on T40")


# ----------------- K-means Clustering with Macqueen -------------------

Y20 <- read.csv(file.choose(),header = T)
Y20[,1]<- NULL
#Y20 <- pr_simil2dist(Y20)
Y20_kmeans <- kmeans(x = Y20, centers = 5, algorithm = "MacQueen")
clusplot(Y20, Y20_kmeans$cluster,col.p = "blue", main = "Angle/Arccos Dissimilarity")
plotcluster(Y20, Y20_kmeans$cluster,xlab = "Coordinate 1", ylab = "Coordinate 2", main = "K-means")


# -----------------Spectral Clustering with knn--------------

Y20_spec <- read.csv(file.choose(),header = T)
Y20_spec[,1]<- NULL
#Y20 <- pr_simil2dist(Y20)
Y20_spec_cluster <- specClust(data = Y20_spec, centers = NULL, nn = 10, method = "symmetric")
clusplot(Y20_spec, Y20_spec_cluster$cluster, main = "Spectral Cluster on Topic 80")
plot(Y20_spec_cluster)
plotcluster(Y20_spec, Y20_spec_cluster$cluster,xlab = "Coordinate 1", ylab = "Coordinate 2", main = "Spectral Cluster on Topic 80")


# ---------- KL Divergence -------------

#Trying KL Divergence -
# KL.divergence(X = Data_KL[2,], Y = Data_KL[1,],k = 19, algorithm = c("kd_tree"))
# X1 <- read.csv(file.choose())
# X1[,1]<- NULL
# X1_1<- X1[1,]
# X1_2 <- X1[2,]

KL_div <- function(x,y) 
  {
    KL.plugin(freqs1 = x, freqs2 = y, unit = 'log2') # Calculates KL divergence between two distributions
  }

# creating my own KL divergence matrix 

KL_mat <- matrix(nrow = 593,ncol = 593) # since the dataset has 593 observations
for(i in 1:593)
{
  for(j in 1:593)
  {
    X1_1 <- X1[i,]
    X1_2 <- X1[j,]
    KL_diverg <- KL.plugin(freqs1 = X1_1, freqs2 = X1_2, unit = 'log2') 
    KL_mat[i,j] <- KL_diverg
  }
}
KL_mat

# ----------JS Divergence-------------

# Trying JS Divergence 

# For 20 Topic model

JS_div <- function(x,y) {
  (0.5 * KL_div(x, (x+y)/2) + 0.5 * KL_div(y, (x+y)/2))  # Formula for JS Divergence
}

# Creating JS Divergence matrix

JD_mat <- matrix(nrow = 593,ncol = 593)
for(i in 500:593)
{
  for(j in 1:593)
  {
    x_jsd <- X1[i,]
    y_jsd <- X1[j,]
    JS_diverg <- JS_div(x = x_jsd, y = y_jsd) 
    JD_mat[i,j] <- JS_diverg
  }
}
head(JD_mat)


write.csv(file = "Desktop/JS_Dis_20.csv", x = JD_mat)

ac_JD_mat_1 <-sqrt(JD_mat)
write.csv(file = "Desktop/JS_Distance_20.csv", x = ac_JD_mat_1)


# For 40 Topic Model 

JD_mat1 <- matrix(nrow = 593,ncol = 593)
for(i in 1:593)
{
  for(j in 1:593)
  {
    x_jsd <- X2[i,]
    y_jsd <- X2[j,]
    JS_diverg <- JS_div(x = x_jsd, y = y_jsd) 
    JD_mat1[i,j] <- JS_diverg
  }
}
head(JD_mat1)

write.csv(file = "Desktop/JS_Dis_40.csv", x = JD_mat1)


# For 60 Topic Model 
JD_mat2 <- matrix(nrow = 593,ncol = 593)
for(i in 1:593)
{
  for(j in 1:593)
  {
    x_jsd <- X3[i,]
    y_jsd <- X3[j,]
    JS_diverg <- JS_div(x = x_jsd, y = y_jsd) 
    JD_mat2[i,j] <- JS_diverg
  }
}
head(JD_mat2)

write.csv(file = "Desktop/JS_Dis_60.csv", x = JD_mat2)



# For 80 Topic Model 
JD_mat3 <- matrix(nrow = 593,ncol = 593)
for(i in 1:593)
{
  for(j in 1:593)
  {
    x_jsd <- X4[i,]
    y_jsd <- X4[j,]
    JS_diverg <- JS_div(x = x_jsd, y = y_jsd) 
    JD_mat3[i,j] <- JS_diverg
  }
}
head(JD_mat3)

write.csv(file = "Desktop/JS_Dis_80.csv", x = JD_mat3)


# -------------COSINE SIMILARITY---------------------

# Creating cosine matrix 
# Topic 20 
COS_mat <- cosine(t(as.matrix(X1))) # Taking the transpose of Topic 20 data set and applying cosine function to get cosine similarity matrix
COS_mat
write.csv(file = "Desktop/Cosine_matrix_T20.csv", x = COS_mat) # Saving matrix as csv file


# -------------ARCOSINE DISSIMILARITY---------------------

# Creating arcosine matrix 
# Topic 20 
ACOS_mat <- acos(COS_mat) # Applying acos function on the cosine similarity matrix
ACOS_mat
write.csv(file = "Desktop/ACosine_matrix_T20.csv", x = ACOS_mat) # Saving matrix as csv file


#----------------Chordal Distance Transformation----------------


Chordal_mat <- sqrt(1-COS_mat) # Applying chordal distance formula to converst cosine similarity matrix to dissimilarity
head(Chordal_mat)
write.csv(file = "Desktop/Chordal_matrix_T20.csv", x = Chordal_mat) # Saving matrix as csv file

# ---------------Creating NULL model--------------

met_data <- read.csv(file.choose())
# topic20_data <-  read.csv(file.choose())
met_data<- met_data[ order(met_data$publication_year), ] # ordering book data based on publishing year
# met_data[,1] <- NULL
# met_data[,2] <- NULL
head(met_data)

met_data_null <- read.csv(file.choose())
met_data_null<- met_data_null[ order(met_data_null$publication_year), ] # ordering book data based on publishing year
head(met_data_null)

# null_data_mat <- matrix(nrow = 593, ncol = 21) 
# 
# for(i in nrow(met_data))
# {
#   for(j in nrow(topic20_data))
#   {
#     if(met_data[i,1] ==  topic20_data[j,1])
#     {
#       null_data_mat[i,2:20]<- topic20_data[j,2:20]
#     }
#     
#   }
# }
# null_data_mat


# ----------------Text To Text-------------------

  t2t_mat <- matrix( nrow = nrow(X1), ncol = 2)
  KL_diverg_cumul = 0
  for(i in 2:nrow(X1))
  {
    book1 <- X1[(i-1),]
    book2 <- X1[i,]
    book1_1 <- met_data[(i-1),]
    book2_1 <- met_data[i,]
    KL_diverg_t2t <- KL_div(x = book1, y = book2) #KL divergence
    KL_diverg_t2t_2 <- KL_div(x = book1_1, y = book2_1) # KL divergence for null model
    KL_diverg_final <- KL_diverg_t2t - KL_diverg_t2t_2
    KL_diverg_cumul <- KL_diverg_cumul + KL_diverg_final
#     j <- i-1
    t2t_mat[i,1] <- i  
    t2t_mat[i,2] <- KL_diverg_cumul
  }

book1=0
book2 = 0
book1_1 = 0
book2_1 = 0

t2t_mat <- t2t_mat[-1,] 
t2t_df <- data.frame(t2t_mat)
ggplot(data = t2t_df, aes(x = t2t_mat[,1], y = t2t_mat[,2], group =1)) + labs(x= "Volumes", y = "Data - Null (bits)", title = "Darwin's Text-to-Text Cumulative Surprise") +geom_line()


# ----------------Past To Text-------------------

p2t_mat <- matrix( nrow = nrow(X1), ncol = 2)
KL_diverg_cumul = 0

for(i in 2:nrow(X1))
{
  book1 = 0
  for(j in 1:(i-1))
  {
    book1 <- book1 + X1[j,]
  }
  book2 <- X1[i,]
  
  book1_1 = 0
  for(j in 1:(i-1))
  {
    book1_1 <- book1_1 + met_data[j,]
  }
  book2_1 <- met_data[i,]
  
  
  KL_diverg_p2t <- KL_div(x = book1, y = book2) # KL divergence between current and previous books
  KL_diverg_p2t_2 <- KL_div(x = book1_1, y = book2_1) # KL divergence between current and previous books with NULL model
  KL_diverg_final <- KL_diverg_p2t - KL_diverg_p2t_2
  KL_diverg_cumul <- KL_diverg_cumul + KL_diverg_final
  #     j <- i-1
  p2t_mat[i,1] <- i  
  p2t_mat[i,2] <- KL_diverg_cumul
}

book1=0
book2 = 0
book1_1 = 0
book2_1 = 0

p2t_mat <- p2t_mat[-1,] 
p2t_df <- data.frame(p2t_mat)

ggplot(data = p2t_df, aes(x = p2t_mat[,1], y = p2t_mat[,2], group =1)) + labs(x= "Volumes", y = "Data - Null (bits)", title = "Darwin's Past-to-Text Cumulative Surprise") + geom_line()


# ----------------Month To Text-------------------

met_data_m2t <- read.csv(file.choose())

m2t_mat <- matrix( nrow = nrow(X1), ncol = 4)
KL_diverg_cumul = 0

for(i in 2:nrow(met_data_m2t))
{
  month <- met_data_m2t[i,3]
  year <- met_data_m2t[i,2]
  book1 = 0
  for(j in 1:nrow(met_data_m2t))
  {
    if(met_data_m2t[j,3] == month && met_data_m2t[j,2] == year) # Filtering books based on month and year
    {
      book1 <- book1 + met_data_m2t[j,4:23]
    }
  }
  
  book2 <- met_data_m2t[i,4:23]
  
  month_2 <- met_data_null[i,24]
  year_2 <- met_data_null[i,23]
  book1_1 = 0
  for(j in 1:nrow(met_data_null))
  {
    if(met_data_null[j,24] == month_2 && met_data_null[j,23] == year_2)
    {
      book1_1 <- book1_1 + met_data_null[j,3:22]
    }
  }
  book2_1 <- met_data_null[i,3:22]
  
  
  KL_diverg_m2t <- KL_div(x = book1, y = book2) # KL divergence between current and previous books 
  KL_diverg_m2t_2 <- KL_div(x = book1_1, y = book2_1) # KL divergence between current and previous books with NULL model
  KL_diverg_final <- KL_diverg_m2t - KL_diverg_m2t_2
  KL_diverg_cumul <- KL_diverg_cumul + KL_diverg_final
  #     j <- i-1
  m2t_mat[i,1] <- i  
  m2t_mat[i,2] <- KL_diverg_cumul
  m2t_mat[i,3] <- month
  m2t_mat[i,4] <- year
}
# book1=0
# book2 = 0
# book1_1 = 0
# book2_1 = 0

m2t_mat <- m2t_mat[-1,] 
m2t_df <- data.frame(m2t_mat)

ggplot(data = m2t_df, aes(x = m2t_mat[,1], y = m2t_mat[,2], group =1)) + labs(x= "Volumes", y = "Data - Null (bits)", title = "Darwin's Month-to-Text Cumulative Surprise")  + geom_line()




# ----------------Year To Text-------------------

met_data_y2t <- read.csv(file.choose())

y2t_mat <- matrix( nrow = nrow(X1), ncol = 3)
KL_diverg_cumul = 0

for(i in 2:nrow(met_data_y2t))
{
  year <- met_data_y2t[i,2]
  book1 = 0
  for(j in 1:nrow(met_data_y2t))
  {
    if(met_data_y2t[j,2] == year)  # Filter only on year
    {
      book1 <- book1 + met_data_y2t[j,4:23]
    }
  }
  
  book2 <- met_data_y2t[i,4:23]
  
  year_2 <- met_data_null[i,23]
  book1_1 = 0
  for(j in 1:nrow(met_data_null))
  {
    if(met_data_null[j,23] == year_2)
    {
      book1_1 <- book1_1 + met_data_null[j,3:22]
    }
  }
  book2_1 <- met_data_null[i,3:22]
  
  
  KL_diverg_y2t <- KL_div(x = book1, y = book2) # KL divergence between current and previous books 
  KL_diverg_y2t_2 <- KL_div(x = book1_1, y = book2_1) # KL divergence between current and previous books with NULL model
  KL_diverg_final <- KL_diverg_y2t - KL_diverg_y2t_2
  KL_diverg_cumul <- KL_diverg_cumul + KL_diverg_final
  #     j <- i-1
  y2t_mat[i,1] <- i  
  y2t_mat[i,2] <- KL_diverg_cumul
  y2t_mat[i,3] <- year
}
# book1=0
# book2 = 0
# book1_1 = 0
# book2_1 = 0

y2t_mat <- y2t_mat[-1,] 
y2t_df <- data.frame(y2t_mat)

ggplot(data = y2t_df, aes(x = y2t_mat[,1], y = y2t_mat[,2], group =1)) + labs(x= "Volumes", y = "Data - Null (bits)", title = "Darwin's Year-to-Text Cumulative Surprise") + geom_line()




# ---------------------------Past-N-to-Text----------------------



n2t_mat <- matrix( nrow = nrow(X1), ncol = 2)
KL_diverg_cumul = 0
n = 100
for(i in 2:nrow(X1))
{
  book1 = 0
  for(j in 1:n)
  {
    book1 <- book1 + X1[j,]
  }
  book2 <- X1[i,]
  
  book1_1 = 0
  for(j in 1:n)
  {
    book1_1 <- book1_1 + met_data[j,]
  }
  book2_1 <- met_data[i,]
  
  
  KL_diverg_n2t <- KL_div(x = book1, y = book2) # KL divergence between current and previous books 
  KL_diverg_n2t_2 <- KL_div(x = book1_1, y = book2_1) # KL divergence between current and previous books with NULL model
  KL_diverg_final <- KL_diverg_n2t - KL_diverg_n2t_2
  KL_diverg_cumul <- KL_diverg_cumul + KL_diverg_final
  #     j <- i-1
  n2t_mat[i,1] <- i  
  n2t_mat[i,2] <- KL_diverg_cumul
}

n2t_mat <- n2t_mat[-1,] 
n2t_df <- data.frame(n2t_mat)

ggplot(data = n2t_df, aes(x = n2t_mat[,1], y = n2t_mat[,2], group =1)) + labs(x= "Volumes", y = "Data - Null (bits)", title = "Darwin's Past-N-to-Text Cumulative Surprise") + geom_line()

# References

# Libraries provided by open source R-community.
# Also used code provided by Prof. Trosset as part of class assignment to calculate graph laplacian and adjacency matrix. 

