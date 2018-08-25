readattachG<-function(genomesfile="data-big/genome.desc", doattach=F)
{
  genomes<-readRDS(genomesfile)
  if(doattach==T){
    genomes<-attachgenomes(genomes)
  }
  return(genomes)
}

attachgenomes<-function (genomes){
    genomes$g <- bigmemory::attach.big.matrix(genomes$genotypes)
    return(genomes)
}
