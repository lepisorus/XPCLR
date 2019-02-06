# Usage: Rscript Phys2Genetic.R -i Andes.chr1.snp -o Andes.chr1.snp.geneticPos

library(optparse)

option_list <- list(make_option(c('-i','--in_file'), action='store', type='character', default=NULL, help='Input file'),
make_option(c('-o','--out_file'), action='store', type='character', default=NULL, help='Output file')
                    )
opt <- parse_args(OptionParser(option_list = option_list))



snp <- read.delim(opt$in_file, header=T)
head(snp)
map<-read.table("../NAM_phasedImputed_1cM_CONVERTED_TO_V3.txt") # load map
head(map)
chrLengths <- c(301476924,237917468,232245527,242062272,217959525,169407836,176826311,175377492,157038028,149632204) # these are full lengths for version 3 genome. 
snp$geneticPosition <- NA
for(i in 1:nrow(snp)){
lowerIndex <- which(map$chrom == snp$chr[i] & map$posV3 <= snp$PhysicalPosition[i])
belowPhys <- max(map$posV3[lowerIndex[length(lowerIndex)]],1,na.rm=T)
belowGen <- max(map$cm[lowerIndex[length(lowerIndex)]],map$cm[which(map$chrom==snp$chr[i])][1]-1,na.rm=T)
higherIndex <- which(map$chrom == snp$chr[i] & map$posV3 >= snp$PhysicalPosition[i])
abovePhys <- min(map$posV3[higherIndex[1]],chrLengths[snp$chr[i]],na.rm=T)
aboveGen <- min(map$cm[higherIndex[1]],map$cm[which(map$chrom==snp$chr[i])][length(which(map$chrom==snp$chr[i]))]+1,na.rm=T)
scale <- {snp$PhysicalPosition[i]-belowPhys}/{abovePhys-belowPhys}
newGen <- {aboveGen-belowGen}*scale + belowGen
snp$geneticPosition[i] <- newGen/100
}
head(snp)
snp <- snp[, c(1:2, 6, 3:5)]
head(snp)
write.table(snp, file=opt$out_file, sep="\t", quote=F, row.names=F, col.names=F)

