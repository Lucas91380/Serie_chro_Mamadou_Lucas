rm(list=ls())
require(astsa)
dataSet <- read.table("Data/mortAlcool",head=TRUE)
head(dataSet)

xt <- ts(dataSet[,3],start=1945) # Consommation d'alcool
yt <- ts(dataSet[,2],start=1945) # Mortalité par cirrhose
head(cbind(yt,xt))

dirName <-"Figures/"
fname <- paste(dirName,"seriesOrg.pdf",sep="")
pdf(fname,width=7, height = 5) #création du fichier pdf
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1) + 0.1)
plot.ts(
  xt,
  las = 1,
  ylab = expression(italic(X[t])),
  xlab = expression(italic(t)),
  main = "Consommation d'alcool",
  frame = FALSE
)
plot.ts(
  yt ,
  las = 1,
  ylab = expression(italic(Y[t])),
  xlab = expression(italic(t)),
  main = "Mortalité liée à la cirrhose",
  frame = FALSE
)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
(corBrute <- cor(xt, yt))