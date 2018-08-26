
# # standard gwa in package
# gwanormal<-GWAS(pheno,geno,plot = F,K=K)
# bkinship<-gwanormal$y
# hist(bkinship, xlab="Inferred allelic selection differentials (b=s)",main="",col="black")
# plot(bkinship ~ genomes$map$physical.pos, xlab="Chromosome position (bp)", pch=16)
# plot(bkinship,s,pch=16)
# 
# bkinap<-(cor(y,X)/abs(cor(y,X))) * sqrt(bkinship) 
# bkinap %>% hist
# plot(bkinap,s)

# library(MASS)
# select(lm.ridge(y ~  (X-1),lambda = seq(0,0.1,0.0001)))
# naive<-lm.ridge(y ~  (X-1),lambda = 1e-04 )
# b<-naive$coef
# plot(b ~ genomes$map$physical.pos, xlab="Chromosome position (bp)")
# plot(b,s)

#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
colnames(M) <- 1:200
geno <- data.frame(marker=1:1000,chrom=rep(1,1000),pos=1:1000,M,check.names=FALSE)
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
u[QTL] <- sample(c(-1,1),size = length(QTL),replace = T)
g <- as.vector(crossprod(M,u))
h2 <- 1
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
pheno <- data.frame(line=1:200,y=y)
scores <- GWAS(pheno,geno,plot=FALSE)

cor(scores$y, u)
plot(scores$y, u)



 
# # standard gwa in package, even worse than kinship and residuals gwa
# gwanormal<-GWAS(pheno,geno,plot = F,K=K)
# bkinship<-gwanormal$y
# hist(bkinship, xlab="Inferred allelic selection differentials (b=s)",main="",col="black")
# plot(bkinship ~ genomes$map$physical.pos, xlab="Chromosome position (bp)", pch=16)
# plot(bkinship,s,pch=16)
# 
# bkinap<-(cor(y,X)/abs(cor(y,X))) * sqrt(bkinship)
# ggdotscolor(c(bkinap),s) %>% moiR::addggregression()





# we can do it also more manually using MCMCglmm 
#Kinship based on Zhou, Carbonetto & Stephens PLOS Genetics 2013
# K=Kcalc(X_)
# dim(K)
# hist(K)

Ainv<-as(ginv(K), "dgCMatrix")
colnames(Ainv)=rownames(Ainv)=genomes$fam$sample.ID

dat=data.frame(y=y, Genotype=genomes$fam$sample.ID) 
head(dat)
Ainv[1:5,1:5]

prior <- list(R = list(V = 1, nu = 0.002), G = list(mu = 0,V = 1e+08))
# R=list(V=diag(2)*(0.002/1.002),nu=1.002),
prior<- list(R = list(V = 1, fix=1), G = list(G1 = list(V = 1, nu = 0.002)))

library(MCMCglmm)
res<-MCMCglmm(data=dat,
              fixed =(y) ~ 1, random=~ Genotype,
              ginverse = list(Genotype=K),
              pl=T,
              pr=T,
              prior=prior
)
prior = prior

y
