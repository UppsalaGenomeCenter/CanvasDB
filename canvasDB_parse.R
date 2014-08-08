
## #######################################################
##
## 1. Functions for parsing SNPs/indels from SOLiD data
##
## #######################################################

parseLifescopeSNPs <- function(SNPfileLifescope){

    data <- read.table(SNPfileLifescope, header=TRUE, as.is=TRUE, sep="\t")

    ref <- data[,"reference"]
    genotype <- data[,"genotype"]

    alt <- genotype
    alt[intersect(which(genotype=="R"),which(ref=="A"))] <- "G"
    alt[intersect(which(genotype=="R"),which(ref=="G"))] <- "A"
    alt[intersect(which(genotype=="Y"),which(ref=="T"))] <- "C"
    alt[intersect(which(genotype=="Y"),which(ref=="C"))] <- "T"
    alt[intersect(which(genotype=="S"),which(ref=="G"))] <- "C"
    alt[intersect(which(genotype=="S"),which(ref=="C"))] <- "G"
    alt[intersect(which(genotype=="W"),which(ref=="A"))] <- "T"
    alt[intersect(which(genotype=="W"),which(ref=="T"))] <- "A"
    alt[intersect(which(genotype=="K"),which(ref=="G"))] <- "T"
    alt[intersect(which(genotype=="K"),which(ref=="T"))] <- "G"
    alt[intersect(which(genotype=="M"),which(ref=="C"))] <- "A"
    alt[intersect(which(genotype=="M"),which(ref=="A"))] <- "C"

    data[,"genotype"] <- alt

    ## The following code handles non-refrence heterozygous calls

    idsForHetNonref <- which(!(data[,"genotype"] %in% c("A","C","G","T")))

    dataHetNonref <- NULL

    for(multid in idsForHetNonref){
        if(data[multid,"genotype"] == "R"){
            newEntry <- data[multid,]
            newEntry["het"] <- "1"
            newEntry["genotype"] <- "G"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
            newEntry["genotype"] <- "A"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
        }
        if(data[multid,"genotype"] == "Y"){
            newEntry <- data[multid,]
            newEntry["het"] <- "1"
            newEntry["genotype"] <- "C"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
            newEntry["genotype"] <- "T"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
        }
        if(data[multid,"genotype"] == "S"){
            newEntry <- data[multid,]
            newEntry["het"] <- "1"
            newEntry["genotype"] <- "G"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
            newEntry["genotype"] <- "C"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
        }
        if(data[multid,"genotype"] == "W"){
            newEntry <- data[multid,]
            newEntry["het"] <- "1"
            newEntry["genotype"] <- "A"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
            newEntry["genotype"] <- "T"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
        }
        if(data[multid,"genotype"] == "K"){
            newEntry <- data[multid,]
            newEntry["het"] <- "1"
            newEntry["genotype"] <- "G"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
            newEntry["genotype"] <- "T"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
        }
        if(data[multid,"genotype"] == "M"){
            newEntry <- data[multid,]
            newEntry["het"] <- "1"
            newEntry["genotype"] <- "C"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
            newEntry["genotype"] <- "A"
            dataHetNonref <- rbind(newEntry,dataHetNonref)
        }
    }

    if(length(idsForHetNonref)>0){
        data <- data[-idsForHetNonref,]
    }

    data <- rbind(data, dataHetNonref)

    SNPdata <- data[,c("Seqid","Start","reference","genotype","coverage","refAlleleCounts","refAlleleStarts","refAlleleMeanQV","novelAlleleCounts","novelAlleleStarts","novelAlleleMeanQV","het")]

    colnames(SNPdata) <- c("chr","pos","ref","alt","coverage","ref.counts","ref.starts","ref.mean.qv","alt.counts","alt.starts","alt.mean.qv","zygosity")

    invisible(unique(SNPdata))

}


parseLifescopeIndels <- function(IndelfileLifescope){

    ## Handle Lifescope gff files (from LS >= 2.5.1)
    if(length(grep(".gff3",IndelfileLifescope)) > 0){
        return(parseLifescopeIndelsGff(IndelfileLifescope))

    }

    ## Handle Lifescope gff files (from LS <= 2.5)
    if(length(grep(".txt",IndelfileLifescope)) > 0){
         return(parseLifescopeIndelsTxt(IndelfileLifescope))
     }

    stop(paste("Unknown file extension for indel file:\n",IndelfileLifescope,"\n",sep=""))

}


parseLifescopeIndelsGff <- function(IndelfileLifescope){

    ## Write perl script to parse gff file...
    tmpfile <- "tmpfile_dir/tmp_indelfile.txt"

    system(paste("utils/parse_lifescope_indels_gff.pl ",IndelfileLifescope," > ",tmpfile,sep=""))

    indelData <- read.table(tmpfile, header=FALSE, as.is=TRUE, sep="\t", comment.char="#")

    file.remove(tmpfile)

    colnames(indelData) <- c("chr","start","end","ref","alt","total.counts","ref.counts","alt.counts","heterozygous")

    indelData[which(indelData[,"heterozygous"] == "HOMOZYGOUS-NON-REF"),"heterozygous"] <- 0
    indelData[which(indelData[,"heterozygous"] == "HETEROZYGOUS-NON-REF"),"heterozygous"] <- 1
    indelData[which(indelData[,"heterozygous"] == "HEMIZYGOUS-REF"),"heterozygous"] <- 1
    indelData[which(indelData[,"heterozygous"] == "HEMIZYGOUS-NON-REF"),"heterozygous"] <- 1

    indelData[which(!indelData[,"heterozygous"] %in% c(0,1)),"heterozygous"] <- NA

    return(indelData)

}


parseLifescopeIndelsTxt <- function(IndelfileLifescope){

    data <- read.table(IndelfileLifescope, header=TRUE, as.is=TRUE, sep="\t")

    indelData <- data[,c("chrom","min.chrom.pos","max.chrom.pos","ref.allele","var.allele1","indel.type","unique.indel.size","num.tot.align","var.counts1","zygosity.call")]

    colnames(indelData) <- c("chr","start","end","ref","alt","type","size","total.counts","alt.counts","heterozygous")

    indelData[,"chr"] <- paste("chr",indelData[,"chr"],sep="")

    if("chr23" %in% indelData[,"chr"])
        indelData[which(indelData[,"chr"] == "chr23"), "chr"] <- "chrX"
    if("chr24" %in% indelData[,"chr"])
        indelData[which(indelData[,"chr"] == "chr24"), "chr"] <- "chrY"
    if("chr25" %in% indelData[,"chr"])
        indelData[which(indelData[,"chr"] == "chr25"), "chr"] <- "chrM"


    for(i in 1:nrow(indelData)){

        if(indelData[i,"type"] == "INSERTION"){
            indelData[i,"type"] <- "ins"
            if(!(length(grep(indelData[i,"ref"], indelData[i,"alt"])) == 0)){
                indelData[i,"alt"] <- sub(indelData[i,"ref"],"",indelData[i,"alt"])
            }
            indelData[i,"ref"] <- "-"
        }

        if(indelData[i,"type"] == "DELETION"){
            indelData[i,"type"] <- "del"
            if(!(length(grep(indelData[i,"alt"], indelData[i,"ref"])) == 0)){
                indelData[i,"ref"] <- sub(indelData[i,"alt"],"",indelData[i,"ref"])
            }
            indelData[i,"alt"] <- "-"
            indelData[i,"start"] <- indelData[i,"start"]-1 ## To make it agree with dbSNP deletions
        }
    }

    indelData[which(indelData[,"heterozygous"] == "HOMOZYGOUS-NON-REF"),"heterozygous"] <- 0
    indelData[which(indelData[,"heterozygous"] == "HETEROZYGOUS-NON-REF"),"heterozygous"] <- 1
    indelData[which(indelData[,"heterozygous"] == "HEMIZYGOUS-REF"),"heterozygous"] <- 1
    indelData[which(indelData[,"heterozygous"] == "HEMIZYGOUS-NON-REF"),"heterozygous"] <- 1

    indelData[which(!indelData[,"heterozygous"] %in% c(0,1)),"heterozygous"] <- NA


    indelData <- indelData[, c("chr","start","end","ref","alt","total.counts","total.counts","alt.counts","heterozygous")]
    colnames(indelData) <- c("chr","start","end","ref","alt","total.counts","ref.counts","alt.counts","heterozygous")

    indelData[,"ref.counts"] <- -1 ## We have no info on reference reads for old style lifescope/bioscope output files..

    invisible(unique(indelData))

}



## ##########################################################
##
## 2. Functions for parsing SNPs/indels from IonTorrent data
##
## ##########################################################



## TODO: Create new version for TSS3.6!!!
parseSNPsFromIonTorrentVCFsOLD <- function(VCFfile, variantCaller="TSS"){

    VCFdata <- read.table(VCFfile, sep="\t", as.is=TRUE)
    colnames(VCFdata) <- c("chr","pos","id","ref","alt","qual","filter","info","format","sample")

    chr <- as.character(VCFdata[,"chr"])
    pos <- as.character(VCFdata[,"pos"])
    id <- as.character(VCFdata[,"id"])
    ref <- as.character(VCFdata[,"ref"])
    alt <- as.character(VCFdata[,"alt"])
    qual <- as.character(VCFdata[,"qual"])
    filter <- as.character(VCFdata[,"filter"])

    ## Important that fields are listed in this order in VCF!! Maybe write a check?
    ## 1-GT, 2-GQ, 3-GL, 4-DP, 5-FDP, 6-AD, 7-APSD, 8-AST, 9-ABQV

    sampleRows <- strsplit(as.character(VCFdata[,"sample"]),":")

    genotype <- sapply(sampleRows, function(x){x[1]})
    genotypeQual <- sapply(sampleRows, function(x){x[2]})
    genotypeLikelihood <- sapply(sampleRows, function(x){x[3]})

    readDepth <- sapply(sampleRows, function(x){x[4]})
    filteredReadDepth <- sapply(sampleRows, function(x){x[5]})

    allelicDepths <- lapply(sampleRows, function(x){strsplit(x[6],",")[[1]]})
    coverageRef <- sapply(allelicDepths, function(x){x[1]})
    coverageAlt <- sapply(allelicDepths, function(x){x[2]})

    allelicPlusStrandDepths <- lapply(sampleRows, function(x){strsplit(x[7],",")[[1]]})
    coveragePlusStrandRef <- sapply(allelicPlusStrandDepths, function(x){x[1]})
    coveragePlusStrandAlt <- sapply(allelicPlusStrandDepths, function(x){x[2]})

    allelicUniqueStarts <- lapply(sampleRows, function(x){strsplit(x[8],",")[[1]]})
    uniqueStartsRef <- sapply(allelicUniqueStarts, function(x){x[1]})
    uniqueStartsAlt <- sapply(allelicUniqueStarts, function(x){x[2]})

    allelicAverageBaseQV <- lapply(sampleRows, function(x){strsplit(x[9],",")[[1]]})
    avgBaseQVRef <- sapply(allelicAverageBaseQV, function(x){x[1]})
    avgBaseQVAlt <- sapply(allelicAverageBaseQV, function(x){x[2]})

    dataParsed <- cbind(chr,pos,ref,alt,readDepth,coverageRef,uniqueStartsRef,avgBaseQVRef,coverageAlt,uniqueStartsAlt,avgBaseQVAlt,genotype)
    colnames(dataParsed) <- c("chr","pos","ref","alt","coverage","ref.counts","ref.starts","ref.mean.qv","alt.counts","alt.starts","alt.mean.qv","zygosity")

    dataParsed[which(dataParsed[,"zygosity"] == "1/1"),"zygosity"] <- 0
    dataParsed[which(dataParsed[,"zygosity"] == "0/1"),"zygosity"] <- 1

    ## Only consider homo/heterozygous calls for now. Fix "1/2" and other strange genotypes later..
    dataParsed <- dataParsed[which(dataParsed[,"zygosity"] %in% c(0,1)),]

    ## Make sure there is only one nucleotide for ref and alt alleles
    dataParsed <- dataParsed[intersect(which(nchar(dataParsed[,"ref"])==1),which(nchar(dataParsed[,"alt"])==1)),]

    return(dataParsed)

}




parseSNPsFromIonTorrentVCFs <- function(VCFfile, variantCaller="TSS"){

    VCFdata <- read.table(VCFfile, sep="\t", as.is=TRUE)
    colnames(VCFdata) <- c("chr","pos","id","ref","alt","qual","filter","info","format","sample")

    chr <- as.character(VCFdata[,"chr"])
    pos <- as.character(VCFdata[,"pos"])
    id <- as.character(VCFdata[,"id"])
    ref <- as.character(VCFdata[,"ref"])
    alt <- as.character(VCFdata[,"alt"])
    qual <- as.character(VCFdata[,"qual"])
    filter <- as.character(VCFdata[,"filter"])

    ## Important that fields are listed in this order in VCF!! Maybe write a check?
    ## 1-GT, 2-GQ, 3-GL, 4-DP, 5-FDP, 6-AD, 7-APSD, 8-AST, 9-ABQV

    ##GT:GQ:DP:FDP:RO:FRO:AO:FAO:SAR:SAF:SRF:SRR:FSAR:FSAF:FSRF:FSRR  0/1:15.0764:7:7:3:3:4:4:1:3:1:2:1:3:1:2

    ## 1-GT, 2-GQ, 3-DP, 4-FDP, 5-RO, 6-FRO, 7-AO, 8-FAO, 9-SAR, 10-SAF, 11-SRF, 12-SRR, 13-FSAR, 14-FSAF, 15-FSRF, 16-FSRR

    sampleRows <- strsplit(as.character(VCFdata[,"sample"]),":")

    genotype <- sapply(sampleRows, function(x){x[1]})
    genotypeQual <- sapply(sampleRows, function(x){x[2]})
    readDepth <- sapply(sampleRows, function(x){x[3]})
    refReads <- sapply(sampleRows, function(x){x[5]})
    refReadsFilt <- sapply(sampleRows, function(x){x[6]})
    altReads <- sapply(sampleRows, function(x){x[7]})
    altReadsFilt <- sapply(sampleRows, function(x){x[8]})

    dataParsed <- cbind(chr,pos,ref,alt,readDepth,refReads,refReadsFilt,genotypeQual,altReads,altReadsFilt,genotypeQual,genotype)
    colnames(dataParsed) <- c("chr","pos","ref","alt","coverage","ref.counts","ref.starts","ref.mean.qv","alt.counts","alt.starts","alt.mean.qv","zygosity")

    dataParsed[which(dataParsed[,"zygosity"] == "1/1"),"zygosity"] <- 0
    dataParsed[which(dataParsed[,"zygosity"] == "0/1"),"zygosity"] <- 1

    ## Only consider homo/heterozygous calls for now. Fix "1/2" and other strange genotypes later..
    dataParsed <- dataParsed[which(dataParsed[,"zygosity"] %in% c(0,1)),]

    ## Make sure there is only one nucleotide for ref and alt alleles
    dataParsed <- dataParsed[intersect(which(nchar(dataParsed[,"ref"])==1),which(nchar(dataParsed[,"alt"])==1)),]

    return(dataParsed)

}


## TODO: Create new version for TSS3.6!!!
parseIndelsFromIonTorrentVCFs <- function(VCFfile, variantCaller="TSS"){

    VCFdata <- read.table(VCFfile, sep="\t", as.is=TRUE)
    colnames(VCFdata) <- c("chr","pos","id","ref","alt","qual","filter","info","format","sample")

    chr <- as.character(VCFdata[,"chr"])
    pos <- as.character(VCFdata[,"pos"])
    id <- as.character(VCFdata[,"id"])
    ref <- as.character(VCFdata[,"ref"])
    alt <- as.character(VCFdata[,"alt"])
    qual <- as.character(VCFdata[,"qual"])
    filter <- as.character(VCFdata[,"filter"])

    ## Important that fields are listed in this order in VCF!! Maybe write a check?
    ## 1-GT, 2-AD, 3-DP, 4-FA, 5-GQ, 6-MQ0, 7-PL
q
    sampleRows <- strsplit(as.character(VCFdata[,"sample"]),":")

    genotype <- sapply(sampleRows, function(x){x[1]})

    allelicDepths <- lapply(sampleRows, function(x){strsplit(x[2],",")[[1]]})
    coverageRef <- sapply(allelicDepths, function(x){x[1]})
    coverageAlt <- sapply(allelicDepths, function(x){x[2]})

    readDepth <- sapply(sampleRows, function(x){x[3]})

    genotypeQual <- sapply(sampleRows, function(x){x[5]})

    uniqueStartsRef <- NA
    uniqueStartsAlt <- NA
    avgBaseQVRef <- NA
    avgBaseQVAlt <- NA

    ## Check that first letter of ref allele is always same as first letter of alt
    if(!length(substring(ref,1,1) == substring(alt,1,1)) == nrow(VCFdata)){
        stop("Error in parsing indel VCF file. First nucleotide of ref and alt alleles doesn't agree")
    }

    ## Remove first nucleotide of ref and alt alleles
    ref <- substring(ref,2)
    alt <- substring(alt,2)

    starts <- pos
    ends <- as.numeric(pos)+nchar(ref)

    dataParsed <- cbind(chr,starts, ends, ref, alt, readDepth, coverageRef, coverageAlt, genotype)
    colnames(dataParsed) <- c("chr","start","end","ref","alt","coverage", "ref.counts","alt.counts","heterozygous")

    dataParsed[which(dataParsed[,"heterozygous"] == "1/1"),"heterozygous"] <- 0
    dataParsed[which(dataParsed[,"heterozygous"] == "0/1"),"heterozygous"] <- 1

    dataParsed[which(!(dataParsed[,"heterozygous"] %in% c(0,1))),] <- NA

    dataParsed[which(dataParsed[,"ref"]==""),"ref"] <- "-"
    dataParsed[which(dataParsed[,"alt"]==""),"alt"] <- "-"

    return(dataParsed)

}



## ##########################################################
##
## 3. Functions for parsing SNPs/indels from Illumina data
##
## ##########################################################



parseSNPsFromIlluminaVCFs <- function(VCFfile, variantCaller="GATK"){

    VCFdata <- read.table(VCFfile, sep="\t", as.is=TRUE)
    colnames(VCFdata) <- c("chr","pos","id","ref","alt","qual","filter","info","format","sample")

    VCFdata <- VCFdata[which(nchar(VCFdata[,"ref"])==1 & nchar(VCFdata[,"alt"])==1),]
    VCFdata <- VCFdata[-1*grep("GL",VCFdata[,1]),]

    chr <- paste("chr",as.character(VCFdata[,"chr"],sep=""))
    chr <- sub("chr ","chr",chr)
    pos <- as.character(VCFdata[,"pos"])
    id <- as.character(VCFdata[,"id"])
    ref <- as.character(VCFdata[,"ref"])
    alt <- as.character(VCFdata[,"alt"])
    qual <- as.character(VCFdata[,"qual"])
    filter <- as.character(VCFdata[,"filter"])

    ## Important that fields are listed in this order in VCF!! Maybe write a check?
    ## 1-GT, 2-AD, 3-DP, 4-GQ, 5-PL

    sampleRows <- strsplit(as.character(VCFdata[,"sample"]),":")

    genotype <- sapply(sampleRows, function(x){x[1]})

    allelicDepths <- lapply(sampleRows, function(x){strsplit(x[2],",")[[1]]})
    coverageRef <- sapply(allelicDepths, function(x){x[1]})
    coverageAlt <- sapply(allelicDepths, function(x){x[2]})

    readDepth <- sapply(sampleRows, function(x){x[3]})
    genotypeQual <- sapply(sampleRows, function(x){x[4]})
    phredScore <- sapply(sampleRows, function(x){x[5]})

    avgBaseQVRef <- NA
    avgBaseQVAlt <- NA
    uniqueStartsRef <- NA
    uniqueStartsAlt <- NA

    dataParsed <- cbind(chr,pos,ref,alt,readDepth,coverageRef,uniqueStartsRef,avgBaseQVRef,coverageAlt,uniqueStartsAlt,avgBaseQVAlt,genotype)
    colnames(dataParsed) <- c("chr","pos","ref","alt","coverage","ref.counts","ref.starts","ref.mean.qv","alt.counts","alt.starts","alt.mean.qv","zygosity")

    dataParsed[which(dataParsed[,"zygosity"] == "1/1"),"zygosity"] <- 0
    dataParsed[which(dataParsed[,"zygosity"] == "0/1"),"zygosity"] <- 1

    ## Only consider homo/heterozygous calls for now. Fix "1/2" and other strange genotypes later..
    dataParsed <- dataParsed[which(dataParsed[,"zygosity"] %in% c(0,1)),]

    return(dataParsed)

}




parseIndelsFromIlluminaVCFs <- function(VCFfile, variantCaller="GATK"){

    VCFdata <- read.table(VCFfile, sep="\t", as.is=TRUE)
    colnames(VCFdata) <- c("chr","pos","id","ref","alt","qual","filter","info","format","sample")

    VCFdata <- VCFdata[which(!(nchar(VCFdata[,"ref"])==1 & nchar(VCFdata[,"alt"])==1)),]
    VCFdata <- VCFdata[-1*grep("GL",VCFdata[,1]),]

    chr <- paste("chr",as.character(VCFdata[,"chr"]),sep="")
    chr <- sub("chr ","chr",chr)
    pos <- as.character(VCFdata[,"pos"])
    id <- as.character(VCFdata[,"id"])
    ref <- as.character(VCFdata[,"ref"])
    alt <- as.character(VCFdata[,"alt"])
    qual <- as.character(VCFdata[,"qual"])
    filter <- as.character(VCFdata[,"filter"])

    ## Important that fields are listed in this order in VCF!! Maybe write a check?
    ## 1-GT, 2-AD, 3-DP, 4-GQ, 5-PL

    sampleRows <- strsplit(as.character(VCFdata[,"sample"]),":")

    genotype <- sapply(sampleRows, function(x){x[1]})

    allelicDepths <- lapply(sampleRows, function(x){strsplit(x[2],",")[[1]]})
    coverageRef <- sapply(allelicDepths, function(x){x[1]})
    coverageAlt <- sapply(allelicDepths, function(x){x[2]})

    readDepth <- sapply(sampleRows, function(x){x[3]})
    genotypeQual <- sapply(sampleRows, function(x){x[4]})
    phredScore <- sapply(sampleRows, function(x){x[5]})

    avgBaseQVRef <- NA
    avgBaseQVAlt <- NA
    uniqueStartsRef <- NA
    uniqueStartsAlt <- NA

    ## Check that first letter of ref allele is always same as first letter of alt
    if(!length(substring(ref,1,1) == substring(alt,1,1)) == nrow(VCFdata)){
        stop("Error in parsing indel VCF file. First nucleotide of ref and alt alleles doesn't agree")
    }

    ## Remove first nucleotide of ref and alt alleles
    ref <- substring(ref,2)
    alt <- substring(alt,2)

    starts <- pos
    ends <- as.numeric(pos)+nchar(ref)

    dataParsed <- cbind(chr,starts, ends, ref, alt, readDepth, coverageRef, coverageAlt, genotype)
    colnames(dataParsed) <- c("chr","start","end","ref","alt","coverage", "ref.counts","alt.counts","heterozygous")

    dataParsed[which(dataParsed[,"heterozygous"] == "1/1"),"heterozygous"] <- 0
    dataParsed[which(dataParsed[,"heterozygous"] == "0/1"),"heterozygous"] <- 1

    dataParsed[which(!(dataParsed[,"heterozygous"] %in% c(0,1))),] <- NA

    dataParsed[which(dataParsed[,"ref"]==""),"ref"] <- "-"
    dataParsed[which(dataParsed[,"alt"]==""),"alt"] <- "-"

    return(dataParsed)

}
