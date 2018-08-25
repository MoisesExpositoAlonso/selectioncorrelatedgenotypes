

readplink<-function (g012file = "data-raw/515g_1000.012", backingfile = "genomes",
    backingpath = "data-big/")
{
    require(bigmemory)
    require(data.table)
    famfile <- sub("\\.012$", ".fam", g012file)
    mapfile <- sub("\\.012$", ".bim", g012file)

    if (!all(c(basename(g012file), basename(mapfile), basename(famfile)) %in%
        list.files(dirname(g012file)))) {
        stop("Not all bim bed and fam files were found!")
    }
    NAMES.MAP <- c("chromosome", "marker.ID", "genetic.dist",
        "physical.pos", "allele1", "allele2")
    NAMES.FAM <- c("family.ID", "sample.ID", "paternal.ID", "maternal.ID",
        "sex", "affection")
    message("Reading the FAM file...")
    fam <- data.table::fread(famfile)
    colnames(fam) <- NAMES.FAM
    message("Reading the BIM/MAP file...")
    bim <- data.table::fread(mapfile)
    colnames(bim) <- NAMES.MAP
    message("Reading the 012 file...")

    geno <- bigmemory::read.big.matrix(filename = g012file, sep = " ",
        # type = "integer", // will dramatically affect implementation
        type = "double",
        backingfile = paste0(backingfile, ".bk"),
        backingpath = backingpath,
        descriptorfile = paste0(backingfile,".desc"),
        skip = 1)

    rda <- paste0(file.path(backingpath, backingfile), ".rda")
    snp_list <- structure(list(genotypes = describe(geno), fam = fam,
        map = bim, savedIn = rda))

    saveRDS(object = snp_list, file = rda)
    message("The genome file is stored in ", rda)
    return(snp_list)
}



read_n_subset_genome<-function(genomefile='databig/genome.rda',selectedsnps){
genomes<-readRDS(genomefile)
map<-genomes$map
map$SNP<-paste0(map$chromosome, "_",map$physical.pos)

Go<-bigmemory::attach.big.matrix(paste0(tools::file_path_sans_ext(genomefile),'.desc'))


sg<- data.frame(Go[,which(map$SNP %in% selectedsnps)])
sg<-sg / 2
head(sg)
colnames(sg) <- selectedsnps

return(sg)
}
