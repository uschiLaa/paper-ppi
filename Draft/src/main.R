## ---- load

library(tourr)
library(tidyverse)
library(reshape2)
library(binostics)
library(gridExtra)
library(wesanderson) #color palette
library(tictoc) #timer
library(mbgraphic) #Katrins package

# function to log transform part of data
asLog <- function(dt){
  dtOut <- dt %>% mutate(BRbsg_ln = log(BRbsg)) %>%
    mutate(BrBsmumu_ln = log(BrBsmumu)) %>%
    mutate(Omegahsq_ln = log(Omegahsq)) %>%
    mutate(Xenonxsec_ln = log(Xenonxsec)) %>%
    mutate(pSDxsec_ln = log(pSDxsec)) %>%
    mutate(gminus2_shiftln = log(gminus2 + 1.8e-9)) %>%
    select(-BRbsg, -BrBsmumu, -Omegahsq, -Xenonxsec, -pSDxsec, -gminus2)
  return(dtOut)
}

#function to plot data and binning of requested projection
plotProj <- function(nProj, fullPath, dMatrix){
  p <- fullPath[nProj][[1]]
  d <- dMatrix %*% p
  s <- scagnostics.default(d[,1],d[,2], bins= 50)
  p1 <- ggplot()+
    geom_point(data = as.data.frame(d), mapping=aes(x=V1, y=V2))
  p2 <- ggplot()+
    geom_point(data = as.data.frame(s$bins), mapping = aes(x=V1, y=V2, color=V3, size=V3))+
    xlim(0,1000)+
    ylim(0,1000)+
    lims(colour = c(1,25), size = c(1,25))+
    scale_color_gradientn(colors=wes_palette(name="Zissou"))+
    scale_size(range = c(1,3))
    grid.arrange(arrangeGrob(p1, widths = 2),arrangeGrob(p2, widths = 7),nrow=1)
    return(s$s)
}

scagIndex <- function(scagType){
  function(mat){
    sR <- scagnostics.default(mat[,1],mat[,2])$s
    return(sR[scagType])
  }
}

splineIndex <- function(){
  function(mat){
    return(splines2d(mat[,1], mat[,2]))
  }
}

dcorIndex <- function(){
  function(mat){
    return(dcor2d(mat[,1], mat[,2]))
  }
}


# get tour path
load("data/fullPath_nObs7.RData")

# read data file
dt <- read.csv("data/obs_Bino.csv")

## ---- scagnostics-summary

n <- 4*100
dfSummary = data.frame( Outlying=rep(0, n), Skewed=rep(0,n), Clumpy=rep(0,n),
                 Sparse=rep(0,n),Striated=rep(0,n), Convex=rep(0,n),
                 Skinny=rep(0,n), Stringy=rep(0,n), Monotonic=rep(0,n), t=rep(0,n), o=rep(0,n), s=rep(0,n))
i <- 1
for(sampleSize in c(500,5000)){ #fix sample size here
  t <- 1
  set.seed(22012018) # all samples include same (smaller) dataset
  dtSample <- asLog(sample_n(dt, sampleSize))
  thisMatrix <- rescale(as.matrix(select(dtSample, -modelName)))
  for(pMatrix in fullPath){
    dProj <- thisMatrix %*% pMatrix
    sgnst <- as.vector(scagnostics.default(dProj[,1],dProj[,2], bins= 50)$s)
    dfSummary[i,] = c(sgnst, t, 1, sampleSize)
    sgnst <- as.vector(scagnostics.default(dProj[,1],dProj[,2], bins= 50, outlierRmv = FALSE)$s)
    dfSummary[i+1,] = c(sgnst, t, 0, sampleSize)
    i = i+2
    t = t+1
    if(t > 100) break # only consider first 100 projections
  }
}

meltR = melt(dfSummary, id = c("t", "o", "s"))

meltR <- mutate(meltR, oR = factor(ifelse(o==1, "On", "Off")))

ggplot() + 
  geom_line(data=subset(meltR, oR=="On"),aes(x=t,y=value, color=as.factor(s),linetype ="On")) +
  geom_line(data=subset(meltR, oR=="Off"),aes(x=t,y=value, color=as.factor(s),linetype = "Off")) +
  geom_vline(xintercept = 19)+
  facet_wrap(~ variable,ncol=1,scales = "free") +
  labs(x = "Projection", y = "Value", color = "Sample Size", linetype = "Outlier Removal")+
  scale_linetype_manual("Outlier Removal",values=c("On"=1,"Off"=2)) #main result for outlier removal = On

## ---- scatter-binning

#reload the 500 point sample
set.seed(22012018)
dt_sample <- asLog(sample_n(dt, 500))
thisMatrix <- rescale(as.matrix(select(dt_sample, -modelName)))

s1 <- plotProj(18, fullPath, thisMatrix)
s2 <- plotProj(19, fullPath, thisMatrix)
s3 <- plotProj(20, fullPath, thisMatrix)

## ---- scatter-table
s1 <- as.data.frame(s1)
s1$s2 <- s2
s1$s3 <- s3
print(s1)
#for some reason kable formatting not working?
#kable(s1, digits=2, format = "markdown")

## ---- scagnositcs-timing
n <- 3*2*100
df = data.frame( t=rep(0, n), nt=rep(0,n), o=rep(0,n), n=rep(0,n))
t <- 1
for(sampleN in c(1,2,3)){
  for(sampleSize in seq(0,10000,by=100)){ 
    set.seed(sampleN*sampleN*10012018) # new seed for each sample selection
    if(sampleSize == 0) sampleSize = 10 # smallest considered sample has 10 points
    dt_sample <- asLog(sample_n(dt, sampleSize))
    thisMatrix <- rescale(as.matrix(select(dt_sample, -modelName))) 
    i <- 1
    tic.clearlog()
    tic() #start timer
    for(pMatrix in fullPath){
      if(i>100) break
      dProj <- thisMatrix %*% pMatrix
      sgnst <- scagnostics.default(dProj[,1],dProj[,2])
      i = i+1
    }
    toc(log=TRUE,quiet=TRUE)
    td <- unlist(tic.log(format=FALSE))["toc.elapsed"]-unlist(tic.log(format=FALSE))["tic.elapsed"]
    df[t,] = c(td, sampleSize, 1, sampleN)
    
    i <- 1
    tic.clearlog()
    tic()
    for(pMatrix in fullPath){
      if(i>100) break
      dProj <- thisMatrix %*% pMatrix
      sgnst <- scagnostics.default(dProj[,1],dProj[,2], outlierRmv = F)
      i = i+1
    }
    toc(log=TRUE,quiet=TRUE)
    td <- unlist(tic.log(format=FALSE))["toc.elapsed"]-unlist(tic.log(format=FALSE))["tic.elapsed"]
    df[t+1,] = c(td, sampleSize, 0, sampleN)
    
    t = t+2
  }
}

ggplot(df, aes(x = nt, y = t, color=as.factor(o))) + 
  geom_line(data=subset(df, n==1),aes(linetype = as.factor(n))) +
  geom_line(data=subset(df, n==2),aes(linetype = as.factor(n))) +
  geom_line(data=subset(df, n==3),aes(linetype = as.factor(n))) +
  labs(x = "Sample Size", y = "Time [s]", color = "Outlier Removal", linetype = "Sample")

## ---- smoothing

dfSmoothing <- dfSummary %>%
  filter(o==1 & s==500) %>%
  select(-o, -s)

ggplot(dfSmoothing) + 
  geom_line(aes(x = t, y = Outlying, colour = "red")) +
  stat_smooth(aes(x = t, y = Outlying), method = lm, formula = y ~ poly(x, 10), se = FALSE)

## ---- guided-tour

#reload the 500 point sample
set.seed(22012018)
dt_sample <- asLog(sample_n(dt, 500))
thisMatrix <- rescale(as.matrix(select(dt_sample, -modelName)))

#record guided tour path based on convex measure
gT <- guided_tour(scagIndex("Convex"))
guidedPath <- save_history(thisMatrix, gT, max_bases = 100)
guidedOriginalOnly <- as.list(guidedPath)
guidedPath <- as.list(interpolate(guidedPath))

## ---- guided-tour-plot

n <- length(guidedPath)
convexDf <- data.frame(Convex=rep(0, n), t=rep(0,n))
i <- 1
for(pMatrix in guidedPath){
  dProj <- thisMatrix %*% pMatrix
  convexDf[i,] <- c(scagnostics.default(dProj[,1],dProj[,2], bins= 50)$s["Convex"],i)
  i = i+1
}

i <- 1
n <- length(guidedOriginalOnly)
convexOriginalOnly <- data.frame(Convex=rep(0, n), t=rep(0,n))
for(pMatrix in guidedOriginalOnly){
  dProj <- thisMatrix %*% pMatrix
  convexOriginalOnly[i,] <- c(scagnostics.default(dProj[,1],dProj[,2], bins= 50)$s["Convex"],i)
  i = i+1
}

ggplot(convexDf, mapping = aes(x=t, y=Convex)) +
  geom_line()+
  ggtitle("Interpolated path")

ggplot(convexOriginalOnly, mapping = aes(x=t, y=Convex)) +
  geom_line()+
  ggtitle("Non-interpolated guided tour planes")

## ---- guided-tour-scatter

s1 <- plotProj(1, guidedPath, thisMatrix)
s2 <- plotProj(13, guidedPath, thisMatrix)
s3 <- plotProj(28, guidedPath, thisMatrix)

## ---- mbgraphics-indices

n <- 200
dfKatrin = data.frame( dcor2d=rep(0, n), splines2d=rep(0,n), t=rep(0,n), s=rep(0,n))
i <- 1
for(sampleSize in c(500,5000)){ #fix sample size here
  t <- 1
  set.seed(22012018) # all samples include same (smaller) dataset
  dtSample <- asLog(sample_n(dt, sampleSize))
  thisMatrix <- rescale(as.matrix(select(dtSample, -modelName)))
  for(pMatrix in fullPath){
    dProj <- thisMatrix %*% pMatrix
    dcor2dV <- dcor2d(dProj[,1],dProj[,2])
    splines2dV <- splines2d(dProj[,1],dProj[,2])
    dfKatrin[i,] = c(dcor2dV, splines2dV, t, sampleSize)
    i = i+1
    t = t+1
    if(t > 100) break # only consider first 100 projections
  }
}

meltR = melt(dfKatrin, id = c("t", "s"))

ggplot(data=meltR) + 
  geom_line(aes(x=t,y=value, color=as.factor(s))) +
  facet_wrap(~ variable,ncol=1,scales = "free") +
  labs(x = "Projection", y = "Value", color = "Sample Size")

## ---- guided-spline2d
#reload the 500 point sample
set.seed(22012018)
dt_sample <- asLog(sample_n(dt, 500))
thisMatrix <- rescale(as.matrix(select(dt_sample, -modelName)))

#record guided tour path based on spline measure
gT <- guided_tour(splineIndex())
guidedPath <- save_history(thisMatrix, gT, max_bases = 100)
guidedOriginalOnly <- as.list(guidedPath)
guidedPath <- as.list(interpolate(guidedPath))

## ---- guided-spline2d-plot
n <- length(guidedPath)
splineDf <- data.frame(spline=rep(0, n), t=rep(0,n))
i <- 1
for(pMatrix in guidedPath){
  dProj <- thisMatrix %*% pMatrix
  splineDf[i,] <- c(splines2d(dProj[,1],dProj[,2]),i)
  i = i+1
}

i <- 1
n <- length(guidedOriginalOnly)
splineOriginalOnly <- data.frame(spline=rep(0, n), t=rep(0,n))
for(pMatrix in guidedOriginalOnly){
  dProj <- thisMatrix %*% pMatrix
  if(i==1) firstProj <- dProj
  splineOriginalOnly[i,] <- c(splines2d(dProj[,1],dProj[,2]),i)
  i = i+1
}

ggplot(splineDf, mapping = aes(x=t, y=spline)) +
  geom_line()+
  ggtitle("Interpolated path")

ggplot(splineOriginalOnly, mapping = aes(x=t, y=spline)) +
  geom_line()+
  ggtitle("Non-interpolated guided tour planes")

p1 <- ggplot()+
  geom_point(data = as.data.frame(firstProj), mapping=aes(x=V1, y=V2))

p2 <- ggplot()+
  geom_point(data = as.data.frame(dProj), mapping=aes(x=V1, y=V2))

grid.arrange(arrangeGrob(p1),arrangeGrob(p2),nrow=1)

## ---- guided-dcor2d
#reload the 500 point sample
set.seed(22012018)
dt_sample <- asLog(sample_n(dt, 500))
thisMatrix <- rescale(as.matrix(select(dt_sample, -modelName)))

#record guided tour path based on spline measure
gT <- guided_tour(dcorIndex())
guidedPath <- save_history(thisMatrix, gT, max_bases = 100)
guidedOriginalOnly <- as.list(guidedPath)
guidedPath <- as.list(interpolate(guidedPath))

## ---- guided-dcor2d-plot
n <- length(guidedPath)
dcorDf <- data.frame(dcor=rep(0, n), t=rep(0,n))
i <- 1
for(pMatrix in guidedPath){
  dProj <- thisMatrix %*% pMatrix
  dcorDf[i,] <- c(dcor2d(dProj[,1],dProj[,2]),i)
  i = i+1
}

i <- 1
n <- length(guidedOriginalOnly)
dcorOriginalOnly <- data.frame(dcor=rep(0, n), t=rep(0,n))
for(pMatrix in guidedOriginalOnly){
  dProj <- thisMatrix %*% pMatrix
  if(i==1) firstProj <- dProj
  dcorOriginalOnly[i,] <- c(dcor2d(dProj[,1],dProj[,2]),i)
  i = i+1
}

ggplot(dcorDf, mapping = aes(x=t, y=dcor)) +
  geom_line()+
  ggtitle("Interpolated path")

ggplot(dcorOriginalOnly, mapping = aes(x=t, y=dcor)) +
  geom_line()+
  ggtitle("Non-interpolated guided tour planes")

p1 <- ggplot()+
  geom_point(data = as.data.frame(firstProj), mapping=aes(x=V1, y=V2))

p2 <- ggplot()+
  geom_point(data = as.data.frame(dProj), mapping=aes(x=V1, y=V2))

grid.arrange(arrangeGrob(p1),arrangeGrob(p2),nrow=1)

print(pMatrix)

