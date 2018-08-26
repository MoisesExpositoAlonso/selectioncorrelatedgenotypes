leadingSNP<-function(dat, significancecol='significance',window=5e4){

colnames(dat) [which(colnames(dat) == significancecol)] <- 'significance'

dat   <- dat %>%
            mutate(posround=round(pos / window)*window) %>%
            mutate(posround=paste(chr,posround,sep='_'))  %>%
            group_by(posround) %>%
            summarize(SNPtop = head(SNP[significance == min(significance)],n=1)  ) ->
            selected

return(selected$SNPtop)
}

# top_gwa<-function(as){
#
#   res<-lapply(as, FUN = function(a){
#                 message("Reading ", a)
#                 tmp<-.read_assoc(a)
#                 tmp<-.selecttop(tmp)
#                 return(tmp) }
#               ) %>%
#     do.call(rbind,.)
#   return(res)
# }
#
# sub_gwa<-function(as,SNPs){
#   res<-lapply(as, FUN = function(a){
#                 message("Reading ", a)
#                 tmp<-.read_assoc(a)
#                 tmp$SNP<-paste(tmp$chr,tmp$pos,sep='_')
#                 tmp<-dplyr::filter(tmp,SNP %in% SNPs)
#                 }
#               ) %>%
#     do.call(rbind,.)
#   return(res)
# }

all_gwa<-function(as){
  res<-lapply(as, FUN = function(a){
                message("Reading ", a)
                d<-.read_assoc(a)
                  return(d)
                }
              ) %>%
    do.call(rbind,.)
  return(res)
}


all_gwa_onlyeff<-function(as,mycol='effect'){
  data(map)
  res<-lapply(as, FUN = function(a){
                message("Reading ", a)
                d<-.read_assoc(a)
                myname<-d$env[1]
                tmp <- data.frame(d[,'SNP'],d[,mycol])
                colnames(tmp) <- c('SNP',myname)
                tmp2<-merge(map[,'SNP'],tmp, by='SNP',all.x=T)
                return(tmp2)
                }
              ) %>%
    do.call(cbind,.)
  ressnp<-res [,'SNP']
  res<-dplyr::select(data.frame(res),-contains('SNP'))
  res<-data.frame(ressnp,res)

  fin<-merge(map,res, by='SNP',all.x=T)

  return(fin)
}


# wrap_gwa<-function(a,...){
#   gwares=.read_assoc(a)
#   aname=.cleanname(a)
#
#   p<-gws::ggmanhattan(gwares,stat.col = 'p_lrt', stat.type = 'p',subset = 5000,...)
#   es<-gws::ggmanhattan(gwares,stat.col = 'beta', stat.type = 'es',subset = 5000,...)
#
#   comb=plot_grid(p, es,nrow=2)
#
#   save_plot(file.path('figs/',paste0(aname,".pdf") ),plot = comb,base_width = 12,base_height = 4 ,useDingbats = F)
#
#   return(list(p,es))
# }


.cleanname<-function(x){
  sapply(x, function(i){
    strsplit(i, split =  '/', fixed=TRUE) %>% unlist %>% tail(1) %>%
    strsplit(split =  '.', fixed=TRUE) %>% unlist %>%
    head(1)
  })
}

.read_hyperparameters<-function(outfile, hyperparameter=1, MCMCchain=1000,quantiles= c(0.5,0.025,0.975) ){

  hyp<-read.table(outfile,header=TRUE, stringsAsFactors = FALSE) %>% tail(n=MCMCchain)

  hist(hyp[,hyperparameter], main = outfile)
  quantile(hyp[,hyperparameter],p=quantiles)

}

.read_assoc<-function(outfile="output/fecundity_mhi.assoc.txt"){
  require(data.table)
  d <- data.table::fread(outfile,header=TRUE,stringsAsFactors = FALSE)
  # d <- read.table(outfile,header=TRUE,stringsAsFactors = FALSE)
  d<-dplyr::rename(d,pos=ps)
  d<-mutate(d, SNP= paste(chr,pos,sep="_"),
            env= .cleanname(outfile),
            cumpos=1:nrow(d))

  if (grepl('.param.', outfile)){
    d<-mutate(d,effect= alpha + (gamma*beta) )
    d<-mutate(d,sparseeffect= (gamma*beta) )
    d$BAY<- d$gamma !=0

  }else if(grepl('.assoc.', outfile)){
    d<-mutate(d,effect= beta)
    d$FDR<- d$p_lrt < fdr_level(d$p_lrt)
    d$BONF<- d$p_lrt < 1/nrow(d)
  }

  return(d)
}

.find_gemmas<-function(folder='output/', type='heritability',...){
  if(type=='heritability'){
    filelist=list.files(folder,pattern = '.hyp.',...)
  }else if(type=='association_bslmm'){
    filelist=list.files(folder,pattern = '.param.',...)
  }else if(type=='association_lm'){
    filelist=list.files(folder,pattern = '.assoc.',...)
  }else if(type=='association_lmm'){
    filelist=list.files(folder,pattern = '.assoc.',...)
  }else{stop("argument not recognized")}
  return(filelist)
}


read_gemmas<-function(folder, what='heritability'){

if(what=='heritability'){
res=sapply(h2,
       .read_heritability)%>% data.frame() %>%
  `colnames<-`(h2%>% .cleanname())
}
else if(what=='association'){

}

  return(res)
}



#' Generate a fam file
#'
#' @param wannaoverwrite Logical. Default False
makefam515<-function(wannaoverwrite=F){
fam=read.table('data-raw/515g.fam')

if(wannaoverwrite==TRUE) devtools::use_data(fam,overwrite = wannaoverwrite)
else return(fam)

}

################################################################################
#' genplink
#' Generates the plink necessary fam file to run gemma
#'
#' @param data
#' @param id
#' @param phenotype
#' @param naval
#' @export
genplink<-function(data,
                   id='id',
                   phenotype,
                   naval='-9',
                   out='test',
                   path='../gemma/',
                   famfile='515g.fam',
                   bimfile='515g.bim',
                   bedfile='515g.bed',
                   cleardir=F,
                   makerelative=F
                   ){

# load the fam
data(fam)
fam<-data.frame(fam)

idname<-ifelse('id' %in% colnames(fam),'id','sample.ID')

# merge with the dataset
myplink<-merge(fam, data[,c('id',phenotype)] , by.x=idname,by.y='id' ,all.x = T)

myplink<-data.frame(
                    sample.ID= myplink[,idname],
                    family.ID= myplink[,'family.ID'],
                    paternal.ID= myplink[,'paternal.ID'],
                    maternal.ID= myplink[,'maternal.ID'],
                    sex= myplink[,'sex'],
                    affection= myplink[,phenotype]
                    )

# make relative if needed
if(makerelative==TRUE){
    myplink[,'affection'] = myplink[,'affection']/mean(myplink[,'affection'],na.rm=T)
}

# check that the NA are as naval
myplink[,'affection'][is.na(myplink[,'affection'])]<-naval # sixth column is the phenotype

# output name
fampath=file.path(path, out)
famfile=file.path(fampath,'/515g.fam')

# create directory
catchexit=system(paste('mkdir',fampath))

if(catchexit !=0){
message("the directory already exists!")
  if(cleardir==TRUE){
    message('flag cleandir=T, so removing directory')
    clear_gemma_dir(out)
    system(paste('mkdir',fampath))
  }else{
    stop('provide permission to remove directory with cleardir=TRUE')
  }
}
# hard link the genome matrix files for gemma
system(paste('ln', file.path(path, bedfile),fampath))
system(paste('ln', file.path(path, bimfile),fampath))


# out
write.table(myplink,file=famfile,col.names=F,row.names=F,quote=F,sep=" ")
print(head(myplink))

}


clear_gemma_dir<-function(out,gemmapath='../gemma'){
  fampath=file.path(gemmapath,out)
  system(paste('rm -rf',fampath))
}




#' Call GEMMA from within R
#'
#' @param out
#'
#' @export
run_gemma<-function(out,type='bslmm',background=TRUE, gemmapath='../gemma',plinkbase='515g',dryrun=FALSE){
background=ifelse(background==F, " ", " &")

if(type=='bslmm'){
  command= paste0('~/bin/gemma -bfile ', file.path(gemmapath,out, plinkbase) ,'  -bslmm -o ' ,out, background)
}else if(type=='lm'){
  command= paste0(' ~/bin/gemma -bfile ', file.path(gemmapath,out, plinkbase) ,'  -lm 4  -o ' ,out, background)
}else if(type=='lmm'){
  command= paste0(' ~/bin/gemma -bfile ', file.path(gemmapath,out, plinkbase) ,'  -gk 1  -o ' ,out, background)
  command2= paste0(' ~/bin/gemma -bfile ', file.path(gemmapath,out, plinkbase) ,'  -k ', paste0(out,'.cXX.txt') , '  -lmm 4 -o ',out, background)
}
if( !type %in% c('bslmm','lm','lmm')) stop('type of GWA not recognized')

# Run gemma
message(paste('running GEMMA command: ',command))
if(dryrun==FALSE) system(command)
# If lmm, additionally run the lmm command
if(type=='lmm'){
}

}


################################################################################

genplinkX<-function(){

}

'plink --bfile mydata --fisher'





# The other option is to run a gemma gwa with a previous transformation into Logit
# https://en.wikipedia.org/wiki/Logit
