#################################################################
##
## File: canvasDB_annotate.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Functions for annotating SNPs and small indels
##
#################################################################


## Import dbSNP data, SIFT scores and PolyPhen2 scores into local databases to enable quick annotations
setupAnnotationDBtables <- function(tmpAnnotationDir=tmpfileDir, dbSNPversion=137, TALK=FALSE){

    ## Import dbSNPcommon rsIds if not already in database...
    dbSNPcommon.table <- annotTables[["dbsnpCommon"]]
    con <- connectToInhouseDB()
    query <- paste("SELECT * FROM ",dbSNPcommon.table,";",sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    if(nrow(tmp) == 0){

        cat("Imporing dbSNP common ids into database...\n")

        dbsnpCommonFile <- paste(tmpAnnotationDir,"hg19_snp",dbSNPversion,"common_parsed.txt",sep="")

	if(!file.exists(dbsnpCommonFile)){
            dbsnpCommonRawFile <- paste(ANNOVARpathDB,"hg19_snp",dbSNPversion,"Common.txt",sep="")

            if(!file.exists(dbsnpCommonRawFile)){
                stop("dbsnpcommon file missing! Download using ANNOVAR.\n")
            }

            parseCmd <- paste("utils/prepare_annot_dbsnp_tables.pl ",dbsnpCommonRawFile," > ",dbsnpCommonFile,sep="")
            system(parseCmd)
        }

        query <- paste("CREATE TABLE ",dbSNPcommon.table,
                       "(SNP_id varchar(20), ",
                       "name varchar(20), ",
                       "PRIMARY KEY(SNP_id)",
                       ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",dbsnpCommonFile,"' REPLACE INTO TABLE ",dbSNPcommon.table,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)
    }


    ## Import dbSNPcommon indel rsIds if not already in database...
    dbSNPcommonIndel.table <- annotTables[["dbsnpCommonIndels"]]
    query <- paste("SHOW TABLES LIKE '",dbSNPcommonIndel.table,"';",sep="")
    con <- connectToInhouseDB()
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    if(nrow(tmp) == 0){

        cat("Imporing dbSNP common indel ids into database...\n")

        dbsnpCommonIndelFile <- paste(tmpAnnotationDir,"hg19_snp",dbSNPversion,"common_indels_parsed.txt",sep="")

        if(!file.exists(dbsnpCommonIndelFile)){
            dbsnpCommonRawFile <- paste(ANNOVARpathDB,"hg19_snp",dbSNPversion,"Common.txt",sep="")

            if(!file.exists(dbsnpCommonRawFile)){
                stop("dbsnpcommon file missing! Download using ANNOVAR.\n")
            }

            parseCmd <- paste("utils/prepare_annot_dbsnp_tables_indels.pl ",dbsnpCommonRawFile," > ",dbsnpCommonIndelFile,sep="")
            system(parseCmd)
        }

        query <- paste("CREATE TABLE ",dbSNPcommonIndel.table,
                       "(indel_id varchar(500), ",
                       "name varchar(20), ",
                       "PRIMARY KEY(indel_id)",
                       ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",dbsnpCommonIndelFile,"' REPLACE INTO TABLE ",dbSNPcommonIndel.table,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)
    }


    ## Import dbSNP rsIds if not already in database...
    dbSNP.table <- annotTables[["dbsnp"]]
    query <- paste("SHOW TABLES LIKE '",dbSNP.table,"';",sep="")
    con <- connectToInhouseDB()
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    if(nrow(tmp) == 0){

        cat("Imporing dbSNP ids into database...\n")

        dbsnpFile <- paste(tmpAnnotationDir,"hg19_snp",dbSNPversion,"_parsed.txt",sep="")

       	if(!file.exists(dbsnpFile)){
            dbsnpRawFile <- paste(ANNOVARpathDB,"hg19_snp",dbSNPversion,".txt",sep="")

            if(!file.exists(dbsnpRawFile)){
                stop("dbsnp file missing! Download using ANNOVAR.\n")
            }

            parseCmd <- paste("utils/prepare_annot_dbsnp_tables.pl ",dbsnpRawFile," > ",dbsnpFile,sep="")
            system(parseCmd)
        }

        query <- paste("CREATE TABLE ",dbSNP.table,
                       "(SNP_id varchar(20), ",
                       "name varchar(20), ",
                       "PRIMARY KEY(SNP_id)",
                       ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",dbsnpFile,"' REPLACE INTO TABLE ",dbSNP.table,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)
    }

    ## Import dbSNP indel rsIds if not already in database...
    dbSNP.tableIndels <- annotTables[["dbsnpIndels"]]
    query <- paste("SHOW TABLES LIKE '",dbSNP.tableIndels,"';",sep="")
    con <- connectToInhouseDB()
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    if(nrow(tmp) == 0){

        cat("Imporing dbSNP indel ids into database...\n")

        dbsnpFileIndels <- paste(tmpAnnotationDir,"hg19_snp",dbSNPversion,"_indels_parsed.txt",sep="")

        if(!file.exists(dbsnpFileIndels)){
            dbsnpRawFile <- paste(ANNOVARpathDB,"hg19_snp",dbSNPversion,".txt",sep="")

            if(!file.exists(dbsnpRawFile)){
                stop("dbsnp file missing! Download using ANNOVAR.\n")
            }

            parseCmd <- paste("utils/prepare_annot_dbsnp_tables_indels.pl ",dbsnpRawFile," > ",dbsnpFileIndels,sep="")
            system(parseCmd)
        }

        query <- paste("CREATE TABLE ",dbSNP.tableIndels,
                       "(indel_id varchar(500), ",
                       "name varchar(20), ",
                       "PRIMARY KEY(indel_id)",
                       ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",dbsnpFileIndels,"' REPLACE INTO TABLE ",dbSNP.tableIndels,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)
    }


    ## Import SNP scores if not already in database...
    score.table <- annotTables[["score"]]
    query <- paste("SHOW TABLES LIKE '",score.table,"';",sep="")
    con <- connectToInhouseDB()
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    if(nrow(tmp) == 0){

        cat("Imporing scores into database...\n")

        scoreFile <- paste(tmpAnnotationDir,"hg19_ljb_all_parsed.txt",sep="")

	if(!file.exists(scoreFile)){
            scoreRawFile <- paste(ANNOVARpathDB,"hg19_ljb_all.txt",sep="")

            if(!file.exists(scoreRawFile)){
                stop("ljb_all file missing! Download using ANNOVAR.\n")
            }

            parseCmd <- paste("utils/prepare_annot_score_tables.pl ",scoreRawFile," > ",scoreFile,sep="")
            system(parseCmd)
        }


        query <- paste("CREATE TABLE ",score.table,
                       "(SNP_id varchar(20), ",
                       "sift double(5,2), ",
                       "polyphen double(5,2), ",
                       "phylop double(5,2), ",
                       "lrt double(5,2), ",
                       "mut_taster double(5,2), ",
                       "gerp double(5,2), ",
                       "PRIMARY KEY(SNP_id)",
                       ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",scoreFile,"' REPLACE INTO TABLE ",score.table,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)
    }

}





## Remove temporary files created by annovar
cleanUpTmpfiles <- function(tmpAnnotationDir){

    tmpfiles <- paste(tmpAnnotationDir,
                      c("tmp_exonicVariantsToBeAnnotated.txt",
                        "tmp_exonicVariantsToBeAnnotated.txt.log",
                        "tmp_variantsToBeAnnotated.txt",
                        "tmp_variantsToBeAnnotated.txt.exonic_variant_function",
                        "tmp_variantsToBeAnnotated.txt.invalid_input",
                        "tmp_variantsToBeAnnotated.txt.log",
                        "tmp_variantsToBeAnnotated.txt.variant_function"),sep="")

    for(tmpfile in tmpfiles){
        if(file.exists(tmpfile)){
            file.remove(tmpfile)
        }
    }
}



## Annotate SNPs with the following information
##
## 1. dbSNP and dbSNPCommon - Assign rsId if available
## 2. refgene - Determine effect on transcript level
## 3. Scores - Amino acid substitution effects
annotateSNPs <- function(inputSNPs, tmpAnnotationDir=".", dbSNPversion=137, TALK=FALSE){

    ps <- proc.time()[3]

    SNPseverity <- array()
    SNPseverity["stopgain"] <- 5
    SNPseverity["nonsynonymous"] <- 4
    SNPseverity["stoploss"] <- 4
    SNPseverity["splicing"] <- 4
    SNPseverity["exonic;splicing"] <- 4
    SNPseverity["ncRNA_exonic"] <- 2
    SNPseverity["ncRNA_intronic"] <- 1
    SNPseverity["intergenic"] <- 1
    SNPseverity["intronic"] <- 1
    SNPseverity["synonymous"] <- 1
    SNPseverity["UTR5"] <- 1
    SNPseverity["UTR3"] <- 1
    SNPseverity["downstream"] <- 1
    SNPseverity["upstream"] <- 1
    SNPseverity["ncRNA_UTR3"] <- 1
    SNPseverity["upstream;downstream"] <- 1
    SNPseverity["ncRNA_splicing"] <- 1
    SNPseverity["exonic"] <- 2
    SNPseverity["ncRNA_UTR5"] <- 1
    SNPseverity["UTR5;UTR3"] <- 1

    inputSNPids <- as.character(inputSNPs[,"SNP_id"])
    inputSNPidsSplit <- strsplit(inputSNPids,"\\|")
    inputChr <- sapply(inputSNPidsSplit,function(x){x[[1]]})
    inputPos <- sapply(inputSNPidsSplit,function(x){x[[2]]})
    inputRef <- sapply(inputSNPidsSplit,function(x){x[[3]]})
    inputAlt <- sapply(inputSNPidsSplit,function(x){x[[4]]})
    inputNrSamples <- as.numeric(inputSNPs[,"nr_samples"])
    inputSampleStr <- as.character(inputSNPs[,"sample_str"])

    colnames <- c("chr","pos","pos","ref","alt","SNP_id")

    tmpDatafile <- paste(tmpAnnotationDir,"/tmp_variantsToBeAnnotated.txt",sep="")
    SNPdata <- unique(cbind(inputChr,inputPos,inputPos,inputRef,inputAlt,inputSNPids,inputNrSamples,inputSampleStr))
    colnames(SNPdata) <- c("chr","pos","pos","ref","alt","SNP_id","nr_samples","samples")
    dataToBeAnnotated <- SNPdata[,colnames,drop=FALSE]
    dataToBeAnnotated[,"chr"] <- sub("chr","",dataToBeAnnotated[,"chr"])
    write.table(dataToBeAnnotated, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    con <- connectToInhouseDB()

    tmpSNP.table <- "tmp_SNP_to_be_annotated"

    query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",tmpSNP.table,
                   "(chr varchar(10), ",
                   "start integer(11), ",
                   "end integer(11), ",
                   "ref varchar(1), ",
                   "alt varchar(1), ",
                   "SNP_id varchar(20), ",
                   "PRIMARY KEY(SNP_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",tmpSNP.table,";",sep="")

    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    dbDisconnect(con)

    dataToBeAdded <- matrix("", nrow=nrow(SNPdata), ncol=19)
    colnames(dataToBeAdded) <- c("SNP_id","chr","pos","ref","alt","nr_samples","samples","dbSNP","dbSNPcommon","class","severity","gene","details","sift","polyphen","phylop","lrt","mut_taster","gerp")
    dataToBeAdded[, c("SNP_id","chr","pos","ref","alt","nr_samples","samples")] <- as.matrix(SNPdata[,c("SNP_id","chr","pos","ref","alt","nr_samples","samples")])


    ## Step 1. Filter against dbSNP and dbSNPCommon
    cat("  - Annotating SNPs with dbSNP rs-ids...",sep="")
    con <- connectToInhouseDB()
    query <- paste("SELECT t1.SNP_id,t2.name FROM ",tmpSNP.table," as t1, ",annotTables[["dbsnp"]]," as t2 WHERE t1.SNP_id=t2.SNP_id;", sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)
    if(nrow(tmp)>0){
    dataToBeAdded[match(tmp[,"SNP_id"],dataToBeAdded[,"SNP_id"]),"dbSNP"] <- tmp[,"name"]
    }

    con <- connectToInhouseDB()
    query <- paste("SELECT t1.SNP_id,t2.name FROM ",tmpSNP.table," as t1, ",annotTables[["dbsnpCommon"]]," as t2 WHERE t1.SNP_id=t2.SNP_id;", sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)
    if(nrow(tmp)>0){
        dataToBeAdded[match(tmp[,"SNP_id"],dataToBeAdded[,"SNP_id"]),"dbSNP"] <- tmp[,"name"]
        dataToBeAdded[match(tmp[,"SNP_id"],dataToBeAdded[,"SNP_id"]),"dbSNPcommon"] <- tmp[,"name"]
    }

    cat(proc.time()[3] - ps,"s\n");


    ## Step 2. Annotate exonic variants against score table
    cat("  - Annotating SNPs with SNP scores...",sep="")
    con <- connectToInhouseDB()
    query <- paste("SELECT t1.SNP_id,t2.sift,t2.polyphen,t2.phylop,t2.lrt,t2.mut_taster,t2.gerp FROM ",tmpSNP.table," as t1, ",annotTables[["score"]]," as t2 WHERE t1.SNP_id=t2.SNP_id;", sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)
    if(nrow(tmp)>0){
        idsToBeAdded <- match(tmp[,"SNP_id"],dataToBeAdded[,"SNP_id"])
        dataToBeAdded[idsToBeAdded,"sift"] <- tmp[,"sift"]
        dataToBeAdded[idsToBeAdded,"polyphen"] <- tmp[,"polyphen"]
        dataToBeAdded[idsToBeAdded,"phylop"] <- tmp[,"phylop"]
        dataToBeAdded[idsToBeAdded,"lrt"] <- tmp[,"lrt"]
        dataToBeAdded[idsToBeAdded,"mut_taster"] <- tmp[,"mut_taster"]
        dataToBeAdded[idsToBeAdded,"gerp"] <- tmp[,"gerp"]
    }
    cat(proc.time()[3] - ps,"s\n");

    con <- connectToInhouseDB()
    query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)

    ## Step 3. Annotate against refSeq genes
    cat("  - Annotating SNPs against refGene...",sep="")
    ANNOVARcmd <- paste(ANNOVARpath," -geneanno -buildver hg19 -dbtype refgene ",tmpDatafile," ",ANNOVARpathDB,sep="")
    system(ANNOVARcmd,ignore.stderr=TRUE)

    variantFunctionFile <- paste(tmpDatafile,".variant_function",sep="")
    exonicVariantFunctionFile <- paste(tmpDatafile,".exonic_variant_function",sep="")

    if(file.info(variantFunctionFile)$size > 0){
        variantFunctions <- read.table(variantFunctionFile, as.is=TRUE)
        variantCategory <- as.character(variantFunctions[,1])
        names(variantCategory) <- as.character(variantFunctions[,8])
        variantGene <- as.character(variantFunctions[,2])
        names(variantGene) <- as.character(variantFunctions[,8])

        dataToBeAdded[match(names(variantGene),dataToBeAdded[,"SNP_id"]),"gene"] <- variantGene
        dataToBeAdded[match(names(variantCategory),dataToBeAdded[,"SNP_id"]),"class"] <- variantCategory
    }

    if(file.info(exonicVariantFunctionFile)$size > 0){
        exonicVariantFunctions <- read.table(exonicVariantFunctionFile, sep="\t")
        exonicVariantFunctions <- exonicVariantFunctions[which(exonicVariantFunctions[,2] != "unknown"),]
        exonicVariantCategory <- as.character(exonicVariantFunctions[,2])
        names(exonicVariantCategory) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])
        exonicVariantDetail <- as.character(exonicVariantFunctions[,3])
        names(exonicVariantDetail) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])

        dataToBeAdded[match(names(exonicVariantCategory),dataToBeAdded[,"SNP_id"]),"class"] <- exonicVariantCategory
        dataToBeAdded[match(names(exonicVariantDetail),dataToBeAdded[,"SNP_id"]),"details"] <- exonicVariantDetail
    }

    dataToBeAdded[,"class"] <- gsub(" SNV","",dataToBeAdded[,"class"])
    dataToBeAdded[, "severity"] <- SNPseverity[dataToBeAdded[, "class"]]

    cat(proc.time()[3] - ps,"s\n");

    con <- connectToInhouseDB()
    query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)

    cleanUpTmpfiles(tmpAnnotationDir)



    return(dataToBeAdded)

}



##
## Annotate indels with the following information
##
## 1. dbSNP and dbSNPCommon - Assign rsId if available
## 2. refgene - Determine effect on transcript level
annotateIndels <- function(inputIndels, tmpAnnotationDir=".", dbSNPversion=137, TALK=FALSE){

    ps <- proc.time()[3]

    indelSeverity <- array()
    indelSeverity["frameshift"] <- 5
    indelSeverity["frameshift deletion"] <- 5
    indelSeverity["frameshift insertion"] <- 5
    indelSeverity["stopgain"] <- 5
    indelSeverity["stopgain SNV"] <- 5
    indelSeverity["stoploss SNV"] <- 5
    indelSeverity["nonframeshift"] <- 3
    indelSeverity["nonframeshift deletion"] <- 3
    indelSeverity["nonframeshift insertion"] <- 3
    indelSeverity["exonic"] <- 3
    indelSeverity["exonic;splicing"] <- 3
    indelSeverity["splicing"] <- 4
    indelSeverity["UTR3"] <- 2
    indelSeverity["UTR5"] <- 2
    indelSeverity["ncRNA_exonic"] <- 2
    indelSeverity["ncRNA_splicing"] <- 2
    indelSeverity["ncRNA_UTR5"] <- 1
    indelSeverity["ncRNA_UTR3"] <- 1
    indelSeverity["ncRNA_intronic"] <- 1
    indelSeverity["intronic"] <- 1
    indelSeverity["intergenic"] <- 1
    indelSeverity["upstream"] <- 1
    indelSeverity["upstream;downstream"] <- 1
    indelSeverity["downstream"] <- 1
    indelSeverity["unknown"] <- 1

    inputIndelIds <- as.character(inputIndels[,"indel_id"])
    inputIndelIdsSplit <- strsplit(inputIndelIds,"\\|")
    inputChr <- sapply(inputIndelIdsSplit,function(x){x[[1]]})
    inputStart <- sapply(inputIndelIdsSplit,function(x){x[[2]]})
    inputEnd <- sapply(inputIndelIdsSplit,function(x){x[[3]]})
    inputRef <- sapply(inputIndelIdsSplit,function(x){x[[4]]})
    inputAlt <- sapply(inputIndelIdsSplit,function(x){x[[5]]})
    inputNrSamples <- as.numeric(inputIndels[,"nr_samples"])
    inputSampleStr <- as.character(inputIndels[,"sample_str"])
    inputType <- rep(NA,length(inputIndelIds))
    inputSize <- rep(NA,length(inputIndelIds))

    colnames <- c("chr","start","end","ref","alt","indel_id")

    tmpDatafile <- paste(tmpAnnotationDir,"/tmp_variantsToBeAnnotated.txt",sep="")
    indelData <- unique(cbind(inputChr,inputStart,inputEnd,inputRef,inputAlt,inputIndelIds,inputNrSamples,inputSampleStr,inputType,inputSize))
    colnames(indelData) <- c("chr","start","end","ref","alt","indel_id","nr_samples","samples","type","size")

    indelDataTmp <- indelData

    ## Create synthetic DNA strings for ref/alt of correct length. Should be ok since ANNOVAR anyway doesn't use the actual sequence..
    insertIds <- which(indelDataTmp[,"ref"] == "-")
    deletionIds <- which(indelDataTmp[,"alt"] == "-")

    if(length(insertIds)>0){
       indelData[insertIds,"type"] <- "ins"
       indelData[insertIds,"size"] <- sapply(insertIds, function(i){nchar(indelDataTmp[i,"alt"])})
    }

    if(length(deletionIds)>0){
       indelData[deletionIds,"type"] <- "del"
       indelData[deletionIds,"size"] <- sapply(deletionIds, function(i){nchar(indelDataTmp[i,"ref"])})
       indelDataTmp[deletionIds,"start"] <- as.numeric(indelDataTmp[deletionIds,"start"])+1 ## Required for ANNOVAR annotation of deletions.
    }

    otherIds <- which(!(indelData[,"type"] %in% c("ins","del")))
    indelData[otherIds,"type"] <- "indel"

    dataToBeAnnotated <- indelDataTmp[,colnames,drop=FALSE]
    dataToBeAnnotated[,"chr"] <- sub("chr","",dataToBeAnnotated[,"chr"])

    write.table(dataToBeAnnotated, file=tmpDatafile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    con <- connectToInhouseDB()

    tmpIndel.table <- "tmp_indels_to_be_annotated"

    query <- paste("DROP TABLE IF EXISTS ",tmpIndel.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",tmpIndel.table,
                   "(chr varchar(10), ",
                   "start integer(11), ",
                   "end integer(11), ",
                   "ref varchar(1), ",
                   "alt varchar(1), ",
                   "indel_id varchar(500), ",
                   "PRIMARY KEY(indel_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",tmpIndel.table,";",sep="")

    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    dbDisconnect(con)

    dataToBeAdded <- matrix("", nrow=nrow(dataToBeAnnotated), ncol=16)
    colnames(dataToBeAdded) <- c("indel_id","chr","start","end","ref","alt","type","size","nr_samples","samples","dbSNP","dbSNPcommon","class","severity","gene","details")
    dataToBeAdded[, c("indel_id","chr","start","end","ref","alt")] <- as.matrix(dataToBeAnnotated[,c("indel_id","chr","start","end","ref","alt")])

    ## Step 1. Filter against dbSNP and dbSNPCommon
    cat("  - Annotating indels with dbSNP rs-ids...",sep="")
    con <- connectToInhouseDB()
    query <- paste("SELECT t1.indel_id,t2.name FROM ",tmpIndel.table," as t1, ",annotTables[["dbsnpIndels"]]," as t2 WHERE t1.indel_id=t2.indel_id;", sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)
    if(nrow(tmp)>0){
    dataToBeAdded[match(tmp[,"indel_id"],dataToBeAdded[,"indel_id"]),"dbSNP"] <- tmp[,"name"]
    }

    con <- connectToInhouseDB()
    query <- paste("SELECT t1.indel_id,t2.name FROM ",tmpIndel.table," as t1, ",annotTables[["dbsnpCommonIndels"]]," as t2 WHERE t1.indel_id=t2.indel_id;", sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)
    if(nrow(tmp)>0){
        dataToBeAdded[match(tmp[,"indel_id"],dataToBeAdded[,"indel_id"]),"dbSNP"] <- tmp[,"name"]
        dataToBeAdded[match(tmp[,"indel_id"],dataToBeAdded[,"indel_id"]),"dbSNPcommon"] <- tmp[,"name"]
    }

    cat(proc.time()[3] - ps,"s\n");

    ## Step 2. Annotate against refSeq genes

    cat("  - Annotating indels against refGene...",sep="")
    ANNOVARcmd <- paste(ANNOVARpath," -geneanno -buildver hg19 -dbtype refgene ",tmpDatafile," ",ANNOVARpathDB,sep="")
    system(ANNOVARcmd,ignore.stderr=TRUE)
    cat(proc.time()[3] - ps,"s\n");

    variantFunctionFile <- paste(tmpDatafile,".variant_function",sep="")
    exonicVariantFunctionFile <- paste(tmpDatafile,".exonic_variant_function",sep="")

    ## #################################################
    ##
    ## Combine all annotations into a single output file
    ##
    ## #################################################

    cat("  - Combining annotations and preparing output...",sep="")

    dataToBeAdded[, c("indel_id","chr","start","end","ref","alt","type","size","nr_samples","samples")] <- as.matrix(indelData[,c("indel_id","chr","start","end","ref","alt","type","size","nr_samples","samples")])

    ## Add refSeq annotations
    if(file.info(variantFunctionFile)$size > 0){
        variantFunctions <- read.table(variantFunctionFile, as.is=TRUE)
        variantCategory <- as.character(variantFunctions[,1])
        names(variantCategory) <- as.character(variantFunctions[,8])
        variantGene <- as.character(variantFunctions[,2])
        names(variantGene) <- as.character(variantFunctions[,8])

        dataToBeAdded[match(names(variantGene),dataToBeAdded[,"indel_id"]),"gene"] <- variantGene
        dataToBeAdded[match(names(variantCategory),dataToBeAdded[,"indel_id"]),"class"] <- variantCategory
    }

    if(file.info(exonicVariantFunctionFile)$size > 0){
        exonicVariantFunctions <- read.table(exonicVariantFunctionFile, sep="\t")
        exonicVariantFunctions <- exonicVariantFunctions[which(exonicVariantFunctions[,2] != "unknown"),]
        exonicVariantCategory <- as.character(exonicVariantFunctions[,2])
        names(exonicVariantCategory) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])
        exonicVariantDetail <- as.character(exonicVariantFunctions[,3])
        names(exonicVariantDetail) <- as.character(exonicVariantFunctions[,ncol(exonicVariantFunctions)])

        dataToBeAdded[match(names(exonicVariantCategory),dataToBeAdded[,"indel_id"]),"class"] <- exonicVariantCategory
        dataToBeAdded[match(names(exonicVariantDetail),dataToBeAdded[,"indel_id"]),"details"] <- exonicVariantDetail
    }

    dataToBeAdded[, "severity"] <- indelSeverity[dataToBeAdded[, "class"]]

    cat(proc.time()[3] - ps,"s\n");

    con <- connectToInhouseDB()
    query <- paste("DROP TABLE IF EXISTS ",tmpIndel.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)
    dbDisconnect(con)

    cleanUpTmpfiles(tmpAnnotationDir)

    return(dataToBeAdded)

}
