
numbat_circle <- read_csv("Draft/data/numbat.csv")
numbat_only <- numbat_circle %>% filter(group=="A") %>% select(-group)
circle_only <- numbat_circle %>% filter(group=="B") %>% select(-group)

scagIndex2 <- function(proj,dataF,scagType) {
  proj <- matrix(proj, ncol=2, byrow=FALSE)
  proj <- tourr:::orthonormalise(proj)
  if(any(is.nan(proj))) return(0)
  mat <- as.matrix(dataF) %*% proj
  #scagType <- "Monotonic"
  sR <- scagnostics.default(mat[,1],mat[,2])$s
  return(sR[scagType])
}

ix <- length(circle_only)
ntries <- 10

lower <- rep(-1,2*ix)
upper <- rep(1, 2*ix)

idxV <- rep(0,ntries)
prjL <- list()

i <- 1
while(i<=ntries){
  par0 <- runif(2*ix, min = 0, max = 1)
  opt <- hjkb(par0, scagIndex2, lower=lower, upper=upper,
              control = list(maximize=TRUE, tol=0.001, 
                             info=TRUE),dataF=circle_only, scagType="Skinny")
  idxV[i] <- opt$value
  prjL[[i]] <- tourr:::orthonormalise(matrix(opt$par, ncol=2, byrow=FALSE))
  i <- i+1
}

i1 <- order(idxV)

print(idxV[i1])

plannedT <- planned_tour(prjL[i1])
quartz()
animate(numbat_only, plannedT)

