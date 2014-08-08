################################################################
##
## File: canvasDB_filter.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Filtering and analysis functions for canvasDB
##
################################################################


## SNP filtering function.
##  - 'filterSamples' is either 'all.other' (default) or an array with sample names
##  - 'discardSamples' is either 'none' (default) or an array with sample names
## The two parameters 'filterSamples' and 'discardSamples' can not be used at the same time
filterSNPs <- function(inSamples, filterSamples="all.other", discardSamples="none", minIn=NA, maxOthers=0, minSeverity=3, noCache=FALSE, dbSNPfilterCommon=FALSE, summaryOnly=FALSE, TALK=FALSE){

    ps <- proc.time()[3]

    con <- connectToInhouseDB()

    if(is.na(minIn)){
        minIn <- length(inSamples)
    }

    SNPsummary.table <- dbTables[["SNP.summary"]]
    sample.table <- dbTables[["sample"]]

    noCacheStr <- ""

    if(noCache){
        "SQL_NO_CACHE"
    }

    ## step 1. Extract sample ids for selected samples
    q1 <-  paste("SELECT ",noCacheStr," sample_id,canvas_id FROM ",sample.table,";",sep="")
    tmp <- dbGetQuery_E(con,q1,TALK=TALK)
    dbDisconnect(con)

    allSampleIds <- tmp[,"sample_id"]
    names(allSampleIds) <- tmp[,"canvas_id"]

    samplesNotInDB <- setdiff(inSamples, names(allSampleIds))

    if(length(samplesNotInDB)>0){
        stop(" Sample(s) not in database: ",paste(samplesNotInDB,collapse=", "),"!")
    }

    selectedSampleIds <- allSampleIds[inSamples]
    filterSampleIds <- NA

    if(length(filterSamples) == 1 && filterSamples=="all.other"){
        filterSampleIds <- allSampleIds[!(names(allSampleIds) %in% names(selectedSampleIds))]
    }
    else{
        if(length(which(!(filterSamples %in% names(allSampleIds))))>0){
            stop("Filter sample(s) not available in sample db table!")
        }
        if(length(which((filterSamples %in% names(selectedSampleIds))))>0){
            stop("Filter sample(s) overlapping with selected sample(s)!")
        }
        filterSampleIds <- allSampleIds[filterSamples]
    }

    if(length(filterSamples) == 0){
        if(is.na(filterSampleIds)){
            stop("Error in filter sample argument!")
        }
    }

    discardSampleIds <- c()

    if(!(length(discardSamples) == 1 && discardSamples=="none")){
        if(filterSamples=="all.other"){
            if(length(which(!(discardSamples %in% names(allSampleIds))))>0){
                stop("Discard sample(s) not available in sample db table!")
            }
            discardSampleIds <- allSampleIds[discardSamples]
            filterSampleIds <- filterSampleIds[!(names(filterSampleIds) %in% names(discardSampleIds))]
        }
        else{
            stop("Error. The discard sample parameter may only be used when filterSamples is set to 'all.other'!")
        }
    }

    #cat(proc.time()[3] - ps,"s\n",sep="");

    ## step 3. Extract SNPids from summary table. Filter on dbSNPcommon if specified..
    maxIn <- length(selectedSampleIds)
    maxTot <- maxIn+maxOthers+length(discardSampleIds)

    filteredSNPdata <- NULL
    samplesForCandidates <- c()

    q3opt1 <- ""
    q3opt2 <- ""

    if(dbSNPfilterCommon){
        q3opt1 <- paste(" AND snp",dbSNPversion,"common=''",sep="")
    }
    if(minSeverity>0){
        q3opt2 <- paste(" AND severity>=",minSeverity,sep="")
    }

    for(i in minIn:maxTot){
        #cat(i,"\n")
        q3 <- paste("SELECT ",noCacheStr," * FROM ",SNPsummary.table,
                    " WHERE nr_samples=",i,q3opt1,q3opt2,";",sep="")
        con <- connectToInhouseDB()
        SNPdata <- dbGetQuery_E(con,q3,TALK=TALK)
        dbDisconnect(con)

        if(nrow(SNPdata) < 1)
            next


        ##cat(proc.time()[3] - ps,"s\n",sep="");
        ## step 4. Filter SNPids based on selected and filter samples. Store results in tmp db table

        sampleStringsUnique <- unique(SNPdata[,"samples"], drop=FALSE)
        sampleMatrix <- matrix(0,nrow=length(sampleStringsUnique),ncol=length(allSampleIds))

        for(j in 1:length(sampleStringsUnique)){
            sampleStr <- sampleStringsUnique[j]
            sampleVec <- eval(parse(text=paste("c(",sampleStr,")",sep="")))
            sampleMatrix[j,sampleVec] <- 1
        }

        nrInSamples <- rowSums(sampleMatrix[,selectedSampleIds, drop=FALSE])
        nrDiscardSamples <- rowSums(sampleMatrix[,discardSampleIds, drop=FALSE])
        nrFilterSamples <- i-(nrInSamples+nrDiscardSamples)

        okSampleStrings <- sampleStringsUnique[(nrInSamples>=minIn) & (nrFilterSamples<=maxOthers)]

        if(length(okSampleStrings)>0){
            for(okSampleString in okSampleStrings){
                SNPidsForSampleStr <- which(SNPdata[,"samples"] == okSampleString)
                okSampleVec <- eval(parse(text=paste("c(",okSampleString,")",sep="")))
                nr_in <- length(intersect(okSampleVec,selectedSampleIds))
                nr_discard <- length(intersect(okSampleVec,discardSampleIds))
                nr_filter <- i-(nr_in+nr_discard)

                samplesForCandidates <- unique(c(okSampleVec,samplesForCandidates))
                okSNPdata <- cbind(SNPdata[SNPidsForSampleStr,,drop=FALSE],nr_in,nr_discard,nr_filter)
                filteredSNPdata <- rbind(filteredSNPdata,okSNPdata)
            }
        }

    }

    if(summaryOnly)
        return(filteredSNPdata)


    if(is.null(filteredSNPdata)){
        return(filteredSNPdata)
    }


    ## Create tmp SNP table for filrered SNP data
    tmpSNP.table <- "tmp_filteredSNPs"

    con <- connectToInhouseDB()

    query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",tmpSNP.table,
                   "(SNP_id varchar(20), ",
                   "chr varchar(10), ",
                   "pos integer(11), ",
                   "ref varchar(1), ",
                   "alt varchar(1), ",
                   "snp",dbSNPversion," varchar(20), ",
                   "class varchar(40), ",
                   "gene varchar(100), ",
                   "details varchar(200), ",
                   "sift varchar(10), ",
                   "polyphen varchar(10), ",
                   "phylop varchar(10), ",
                   "lrt varchar(10), ",
                   "mut_taster varchar(10), ",
                   "gerp varchar(10), ",
                   "PRIMARY KEY(SNP_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    tmpDatafile <- "tmp_filteredSNP.txt"

    dataCols <- c("SNP_id","chr","pos","ref","alt",paste("snp",dbSNPversion,sep=""),"class","gene","details","sift","polyphen","phylop","lrt","mut_taster","gerp")

    write.table(filteredSNPdata[,dataCols], file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

    query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",tmpSNP.table,";",sep="")
    tmp <- dbGetQuery_E(con ,query, TALK=TALK)
    file.remove(tmpDatafile)

    SNPresults <- NULL

    ## step 5. Loop through samples and extract final SNP table..
    for(sample in samplesForCandidates){
        SNPdata.table <- paste("snp_data_",sample,sep="")
        sample_name <- names(allSampleIds)[match(sample, allSampleIds)]

        q5 <- paste("SELECT ",noCacheStr," t1.SNP_id, t1.class,t1.chr,t1.pos,t1.ref,t1.alt,snp",dbSNPversion,",t1.gene, t2.heterozygous, t1.sift,t1.polyphen,t1.phylop,t1.lrt,t1.mut_taster,t1.gerp,t2.coverage,t2.ref_counts,t2.ref_starts,t2.alt_counts,t2.alt_starts,t1.details FROM ",tmpSNP.table," as t1, ",SNPdata.table," as t2 ",
                    " WHERE t1.SNP_id=t2.SNP_id;",sep="")

        SNPresultsForSample <- dbGetQuery_E(con, q5, TALK=TALK)

        if(nrow(SNPresultsForSample)>0){
            SNPresults <- rbind(SNPresults,cbind(sample_name, SNPresultsForSample))
        }

    }

    query <- paste("DROP TABLE IF EXISTS ",tmpSNP.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    dbDisconnect(con)

    totalTime <- proc.time()[3] - ps

    cat("Total time for filtering: ",totalTime,"s\n",sep="");

    SNPresults <- SNPresults[order(SNPresults[,"pos"]),]
    SNPresults <- SNPresults[order(SNPresults[,"chr"]),]

    return(SNPresults)

}


## Indel filtering function.
##  - 'filterSamples' is either 'all.other' (default) or an array with sample names
##  - 'discardSamples' is either 'none' (default) or an array with sample names
## The two parameters 'filterSamples' and 'discardSamples' can not be used at the same time
filterIndels <- function(inSamples, filterSamples="all.other", discardSamples="none", minIn=NA, maxOthers=0, minSeverity=3, summaryOnly=FALSE,  dbSNPfilterCommon=FALSE, TALK=FALSE){

    ps <- proc.time()[3]

    con <- connectToInhouseDB()

    if(is.na(minIn)){
        minIn <- length(inSamples)
    }

    indelSummary.table <- dbTables[["indel.summary"]]
    sample.table <- dbTables[["sample"]]

    ## step 1. Extract sample ids for selected samples
    q1 <-  paste("SELECT sample_id,canvas_id FROM ",sample.table,";",sep="")
    tmp <- dbGetQuery_E(con,q1,TALK=TALK)
    dbDisconnect(con)

    allSampleIds <- tmp[,"sample_id"]
    names(allSampleIds) <- tmp[,"canvas_id"]

    samplesNotInDB <- setdiff(inSamples, names(allSampleIds))

    if(length(samplesNotInDB)>0){
        stop(" Sample(s) not in database: ",paste(samplesNotInDB,collapse=", "),"!")
    }

    selectedSampleIds <- allSampleIds[inSamples]

    filterSampleIds <- NA

    if(length(filterSamples) == 1 && filterSamples=="all.other"){
        filterSampleIds <- allSampleIds[!(names(allSampleIds) %in% names(selectedSampleIds))]
    }
    else{
        if(length(which(!(filterSamples %in% names(allSampleIds))))>0){
            stop("Filter sample(s) not available in sample db table!")
        }
        if(length(which((filterSamples %in% names(selectedSampleIds))))>0){
            stop("Filter sample(s) overlapping with selected sample(s)!")
        }
        filterSampleIds <- allSampleIds[filterSamples]
    }

    if(length(filterSamples) == 0){
        if(is.na(filterSampleIds)){
            stop("Error in filter sample argument!")
        }
    }

    discardSampleIds <- c()

    if(!(length(discardSamples) == 1 && discardSamples=="none")){
        if(filterSamples=="all.other"){
            if(length(which(!(discardSamples %in% names(allSampleIds))))>0){
                stop("Discard sample(s) not available in sample db table!")
            }
            discardSampleIds <- allSampleIds[discardSamples]
            filterSampleIds <- filterSampleIds[!(names(filterSampleIds) %in% names(discardSampleIds))]
        }
        else{
            stop("Error. The discard sample parameter may only be used when filterSamples is set to 'all.other'!")
        }
    }

    ## step 3. Extract SNPids from summary table. Filter on dbSNPcommon if specified..
    maxIn <- length(selectedSampleIds)
    maxTot <- maxIn+maxOthers+length(discardSampleIds)

    filteredIndelData <- NULL
    samplesForCandidates <- c()

    q3opt1 <- ""
    q3opt2 <- ""

    if(dbSNPfilterCommon){
        q3opt1 <- paste(" AND snp",dbSNPversion,"common=''",sep="")
    }
    if(minSeverity>0){
        q3opt2 <- paste(" AND severity>=",minSeverity,sep="")
    }

    for(i in minIn:maxIn){
        cat(i,"\n")
        q3 <- paste("SELECT  * FROM ",indelSummary.table,
                    " WHERE nr_samples=",i,q3opt1,q3opt2,";",sep="")
        con <- connectToInhouseDB()
        indelData <- dbGetQuery_E(con,q3,TALK=TALK)
        dbDisconnect(con)

        if(nrow(indelData) < 1)
            next


        cat(proc.time()[3] - ps,"s\n",sep="");
        ## step 4. Filter SNPids based on selected and filter samples. Store results in tmp db table

        sampleStringsUnique <- unique(indelData[,"samples"])
        sampleMatrix <- matrix(0,nrow=length(sampleStringsUnique),ncol=length(allSampleIds))

        for(j in 1:length(sampleStringsUnique)){
            sampleStr <- sampleStringsUnique[j]
            sampleVec <- eval(parse(text=paste("c(",sampleStr,")",sep="")))
            sampleMatrix[j,sampleVec] <- 1
        }

        nrInSamples <- rowSums(sampleMatrix[,selectedSampleIds,drop=FALSE])
        nrDiscardSamples <- rowSums(sampleMatrix[,discardSampleIds,drop=FALSE])
        nrFilterSamples <- i-(nrInSamples+nrDiscardSamples)

        okSampleStrings <- sampleStringsUnique[(nrInSamples>=minIn) & (nrFilterSamples<=maxOthers)]

        if(length(okSampleStrings)>0){
            for(okSampleString in okSampleStrings){
                indelIdsForSampleStr <- which(indelData[,"samples"] == okSampleString)
                okSampleVec <- eval(parse(text=paste("c(",okSampleString,")",sep="")))
                nr_in <- length(intersect(okSampleVec,selectedSampleIds))
                nr_discard <- length(intersect(okSampleVec,discardSampleIds))
                nr_filter <- i-(nr_in+nr_discard)

                samplesForCandidates <- unique(c(okSampleVec,samplesForCandidates))
                okIndelData <- cbind(indelData[indelIdsForSampleStr,,drop=FALSE],nr_in,nr_discard,nr_filter)
                filteredIndelData <- rbind(filteredIndelData,okIndelData)
            }
        }

    }

    if(summaryOnly)
        return(filteredIndelData)

    if(is.null(filteredIndelData)){
        return(filteredIndelData)
    }

    ## Create tmp SNP table for filrered SNP data
    tmpIndel.table <- "tmp_filteredIndels"

    con <- connectToInhouseDB()

    query <- paste("DROP TABLE IF EXISTS ",tmpIndel.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",tmpIndel.table,
                   "(indel_id varchar(500), ",
                   "chr varchar(10), ",
                   "start integer(11), ",
                   "end integer(11), ",
                   "ref varchar(500), ",
                   "alt varchar(500), ",
                   "snp",dbSNPversion," varchar(20), ",
                   "type varchar(100), ",
                   "size integer(10), ",
                   "class varchar(40), ",
                   "gene varchar(100), ",
                   "details varchar(200), ",
                   "PRIMARY KEY(indel_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    tmpDatafile <- "tmp_filteredIndels.txt"

    dataCols <- c("indel_id","chr","start","end","ref","alt",paste("snp",dbSNPversion,sep=""),"type","size","class","gene","details")

    write.table(filteredIndelData[,dataCols], file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

    query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",tmpIndel.table,";",sep="")
    tmp <- dbGetQuery_E(con ,query, TALK=TALK)
    file.remove(tmpDatafile)

    indelResults <- NULL

    ## step 5. Loop through samples and extract final SNP table..
    for(sample in samplesForCandidates){

        indelData.table <- paste("indel_data_",sample,sep="")
        sample_name <- names(allSampleIds)[match(sample, allSampleIds)]

        ## step 5. Extract final Indel table
        q5 <- paste("SELECT t1.indel_id,t1.class,t1.chr,t1.start,t1.end,t1.ref,t1.alt,t1.snp",dbSNPversion,",t1.type,t1.gene,t2.heterozygous,t2.coverage,t2.alt_counts, t1.details FROM ",tmpIndel.table," as t1, ",indelData.table," as t2 ",
                    " WHERE t1.indel_id=t2.indel_id;",sep="")

        indelResultsForSample <- dbGetQuery_E(con, q5, TALK=TALK)

        if(nrow(indelResultsForSample)>0){
            indelResults <- rbind(indelResults,cbind(sample_name, indelResultsForSample))
        }

    }


    query <- paste("DROP TABLE IF EXISTS ",tmpIndel.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    dbDisconnect(con)

    indelResults <- indelResults[order(indelResults[,"start"]),]
    indelResults <- indelResults[order(indelResults[,"chr"]),]


    cat(proc.time()[3] - ps,"s\n");

    return(indelResults)

}




getStatsForSample <- function(sample, TALK=FALSE){

    SNPsummary.table <- dbTables[["SNP.summary"]]

    run.table <- dbTables[["run"]]

    con <- connectToInhouseDB()

    ## step 1. Extract sample ids for selected samples
    q1 <-  paste("SELECT sample_id,canvas_id,total_reads,reads_on_target FROM ",run.table,
                 " WHERE canvas_id='",sample,"';",sep="")

    runData <- dbGetQuery_E(con,q1,TALK=TALK)
    SNPdata.table <- paste("snp_data_",runData[1,"sample_id"],sep="")


    sampleId <- runData[1,"sample_id"]

    q2 <- paste("SELECT t2.* FROM ",SNPdata.table," AS t1, ",SNPsummary.table," AS t2",
                " WHERE t1.SNP_id=t2.SNP_id;",sep="")

    SNPdata <- dbGetQuery_E(con,q2,TALK=TALK)

    dbDisconnect(con)

    totalSNPs <- nrow(SNPdata)
    SNPsInDBSNP <- length(which(SNPdata[,paste("snp",dbSNPversion,sep="")]!=""))

    res <- c(runData[1,"sample_id"], runData[1,"canvas_id"], runData[1,"total_reads"], runData[1,"reads_on_target"], totalSNPs, SNPsInDBSNP)
    names(res) <- c("sample_id", "canvas_id", "total_reads", "reads_on_target", "snps_total", "snps_in_dbsnp")

    return(res)

}



getSNPsAndIndelsForSample <- function(sample, getIndels=TRUE, TALK=FALSE){

    SNPsummary.table <- dbTables[["SNP.summary"]]
    indelSummary.table <- dbTables[["indel.summary"]]

    run.table <- dbTables[["run"]]

    con <- connectToInhouseDB()

    ## step 1. Extract sample ids for selected samples
    q1 <-  paste("SELECT sample_id,canvas_id,total_reads,reads_on_target FROM ",run.table,
                 " WHERE canvas_id='",sample,"';",sep="")

    runData <- dbGetQuery_E(con,q1,TALK=TALK)
    SNPdata.table <- paste("snp_data_",runData[1,"sample_id"],sep="")

    sampleId <- runData[1,"sample_id"]

    q2 <- paste("SELECT t1.*,t2.snp137,t2.class,t2.gene,t2.sift,t2.polyphen,t2.phylop,t2.lrt,t2.mut_taster,t2.gerp FROM ",SNPdata.table," AS t1, ",SNPsummary.table," AS t2",
                " WHERE t1.SNP_id=t2.SNP_id;",sep="")

    SNPdata <- dbGetQuery_E(con,q2,TALK=TALK)

    indelData <- NULL

    if(getIndels){

       indelData.table <- paste("indel_data_",runData[1,"sample_id"],sep="")

       q2 <- paste("SELECT t1.*,t2.snp137,t2.class,t2.gene FROM ",indelData.table," AS t1, ",indelSummary.table," AS t2",
                " WHERE t1.indel_id=t2.indel_id;",sep="")

       indelData <- dbGetQuery_E(con,q2,TALK=TALK)

       dbDisconnect(con)
    }


    invisible(list("SNPs"=SNPdata, "indels"=indelData))
}


getUgcIdsByGeographicLocation <- function(geographicLocation, TALK=FALSE){

    sample.table <- dbTables[["sample"]]

    con <- connectToInhouseDB()

    q1 <-  paste("SELECT canvas_id FROM ",sample.table,
                  " WHERE geographic_location='",geographicLocation,"';",sep="")

    canvasIds <- dbGetQuery_E(con,q1,TALK=TALK)

    dbDisconnect(con)

    return(canvasIds[,"canvas_id"])

}

getUgcIdsWithMappingStats <- function(TALK=FALSE){

    run.table <- dbTables[["run"]]

    con <- connectToInhouseDB()

    q1 <-  paste("SELECT canvas_id FROM ",run.table,
                  " WHERE total_reads>0;",sep="")

    canvasIds <- dbGetQuery_E(con,q1,TALK=TALK)

    dbDisconnect(con)

    return(canvasIds[,"canvas_id"])

}


## Function for filtering for recessive diseases: homozygous & compound heterozygous variants
filterRecessive <- function(inSamples, minSeverity=3, dbSNPfilterCommon=FALSE, maxFreq=0.01, TALK=FALSE, includeCompoundHeteroz=TRUE, outfile=NA){

    maxOthers <- 0

    if(maxFreq>0){
        sample.table <- dbTables[["sample"]]
        con <- connectToInhouseDB()
        q <-  paste("SELECT count(*) canvas_id FROM ",sample.table,";",sep="")
        totalNrSamples <- dbGetQuery_E(con,q,TALK=TALK)
        dbDisconnect(con)

        maxOthers <- as.numeric(round((totalNrSamples-length(inSamples))*maxFreq))
    }

    data <- filterSNPs(inSamples, minSeverity=minSeverity, dbSNPfilterCommon=dbSNPfilterCommon, maxOthers=maxOthers)

    nrSamples <- length(inSamples)

    ## 1. Get shared homozygous variants
    candSNPs <- NULL
    tmp <- data[data[,"heterozygous"] == 0,]
    if(!is.null(tmp)){
        if(nrow(tmp)>0){
            tmp2 <- tmp[sapply(tmp[,"SNP_id"], function(x){ length(which(tmp[,"SNP_id"] %in% x)) == nrSamples }),]
            if(nrow(tmp2)>0){
                tmp3 <- tmp2[which(tmp2[,"sample_name"] %in% inSamples),]
                if(nrow(tmp3)>0){
                    candSNPs <- tmp3[sapply(tmp3[,"SNP_id"], function(x){ length(which(tmp3[,"SNP_id"] %in% x)) == nrSamples }),]
                }
            }
        }
    }

    ## 2. Get shared indels
    dataIndels <- filterIndels(inSamples, minSeverity=minSeverity, dbSNPfilterCommon=dbSNPfilterCommon, maxOthers=maxOthers)
    candIndels <- NULL
    tmp <- dataIndels[dataIndels[,"heterozygous"] == 0,]
    if(!is.null(tmp)){
        if(nrow(tmp)>0){
            tmp2 <- tmp[sapply(tmp[,"indel_id"], function(x){ length(which(tmp[,"indel_id"] %in% x)) == nrSamples }),]
            if(nrow(tmp2)>0){
                tmp3 <- tmp2[which(tmp2[,"sample_name"] %in% inSamples),]
                if(nrow(tmp3)>0){
                    candIndels <- tmp3[sapply(tmp3[,"indel_id"], function(x){ length(which(tmp3[,"indel_id"] %in% x)) == nrSamples }),]
                }
            }
        }
    }


    ## 3. Get shared compound heterozygous variants (could be both SNPs and indels)
    candCompHeteroz <- list()

    if(length(inSamples) == 1){
        includeCompoundHeteroz <- FALSE
    }

    if(includeCompoundHeteroz){
        tmpSNPs <- data[data[,"heterozygous"] == 1,]
        tmpSNPs <- tmpSNPs[which(tmpSNPs[,"sample_name"] %in% inSamples),]
        tmpIndels <- dataIndels[dataIndels[,"heterozygous"] == 1,]
        tmpIndels <- tmpIndels[which(tmpIndels[,"sample_name"] %in% inSamples),]
        candGenes <- unique(c(tmpSNPs[,"gene"],tmpIndels[,"gene"]))

        for(gene in candGenes){
            nrMutsInGene <- length(which(gene == tmpSNPs[,"gene"])) + length(which(gene == tmpIndels[,"gene"]))
            if(nrMutsInGene >= nrSamples*2){
                candCompHeteroz[[gene]] <- list()
                if(length(which(tmpSNPs[,"gene"] %in% gene)>0))
                    candCompHeteroz[[gene]][["snps"]] <- tmpSNPs[which(tmpSNPs[,"gene"] %in% gene),]
                if(length(which(tmpIndels[,"gene"] %in% gene)>0))
              candCompHeteroz[[gene]][["indels"]] <- tmpIndels[which(tmpIndels[,"gene"] %in% gene),]
            }
        }
    }

    cand <- list()

    if(!is.null(candSNPs)){
        if(nrow(candSNPs) == 0)
            candSNPs <- NULL
    }
    if(!is.null(candIndels)){
        if(nrow(candIndels) == 0)
            candIndels <- NULL
    }

    cand[["homoz_snps"]] <- candSNPs
    cand[["homoz_indels"]] <- candIndels
    cand[["comp_heteroz"]] <- candCompHeteroz

    if(!is.na(outfile)){
        cat("*********** CanvasDB filtering analysis report ***********\n", file=outfile)
        cat("**  Analysis of recessive variants in samples: ",paste(inSamples, collapse=",",sep=""),"\n", file=outfile, append=TRUE)
        cat("**  Variants allowed to be detected in at most ",maxOthers," (",100*maxFreq,"%) of other samples\n", sep="", file=outfile, append=TRUE)
        cat("**********************************************************\n\n", file=outfile, append=TRUE)
        cat("-------------------  Homozygous SNPs ---------------------\n", file=outfile, append=TRUE)
        if(!(is.null(cand[["homoz_snps"]]))){
            cat(colnames(cand[["homoz_snps"]]), sep="\t", file=outfile, append=TRUE)
            cat("\n", file=outfile, append=TRUE)
            write.table(cand[["homoz_snps"]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
        }
        else{
            cat("no shared homozygous SNPs found", file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
        cat("-------------------  Homozygous Indels ---------------------\n", file=outfile, append=TRUE)
        if(!(is.null(cand[["homoz_indels"]]))){
            cat(colnames(cand[["homoz_indels"]]), sep="\t", file=outfile, append=TRUE)
            cat("\n", file=outfile, append=TRUE)
            write.table(cand[["homoz_indels"]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
        }
        else{
            cat("no shared homozygous indels found", file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
        cat("-------------------  Compound heterozygous variants --------\n", file=outfile, append=TRUE)
        if(length(cand[["comp_heteroz"]])>0){
            for(i in 1:length(cand[["comp_heteroz"]])){
                cat(names(cand[["comp_heteroz"]])[i],":\n", file=outfile, append=TRUE)
                for(type in names(cand[["comp_heteroz"]][[i]])){
                    write.table(cand[["comp_heteroz"]][[i]][[type]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
                }
            }
        }
        else{
            cat("no shared compound heterozygous variants found", file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
    }

    return(cand)

}


## Function for filtering for recessive diseases: homozygous & compound heterozygous variants
filterDominant <- function(inSamples, dbSNPfilterCommon=FALSE, maxFreq=0, minSeverity=3, TALK=FALSE, outfile=NA){

    maxOthers <- 0

    if(maxFreq>0){
        sample.table <- dbTables[["sample"]]
        con <- connectToInhouseDB()
        q <-  paste("SELECT count(*) canvas_id FROM ",sample.table,";",sep="")
        totalNrSamples <- dbGetQuery_E(con,q,TALK=TALK)
        dbDisconnect(con)

        maxOthers <- as.numeric(round((totalNrSamples-length(inSamples))*maxFreq))
    }

    data <- filterSNPs(inSamples, minSeverity=minSeverity, dbSNPfilterCommon=dbSNPfilterCommon, maxOthers=maxOthers)

    nrSamples <- length(inSamples)

    ## 1. Get shared heterozygous SNPs
    candSNPs <- NULL
    tmp <- data[data[,"heterozygous"] == 1,]
    if(!is.null(tmp)){
        if(nrow(tmp)>0){
            tmp2 <- tmp
            ##tmp2 <- tmp[sapply(tmp[,"SNP_id"], function(x){ length(which(tmp[,"SNP_id"] %in% x)) >= nrSamples }),]
            if(nrow(tmp2)>0){
                tmp3 <- tmp2[which(tmp2[,"sample_name"] %in% inSamples),]
                if(nrow(tmp3)>0){
                    candSNPs <- tmp3[sapply(tmp3[,"SNP_id"], function(x){ length(which(tmp3[,"SNP_id"] %in% x)) == nrSamples }),]
                }
            }
        }
    }

    ## 2. Get shared heterozygous indels
    dataIndels <- filterIndels(inSamples, minSeverity=minSeverity, dbSNPfilterCommon=dbSNPfilterCommon, maxOthers=maxOthers)
    candIndels <- NULL
    tmp <- dataIndels
    tmp <- dataIndels[dataIndels[,"heterozygous"] == 1,]
    if(!is.null(tmp)){
        if(nrow(tmp)>0){
            tmp2 <- tmp
            ##tmp2 <- tmp[sapply(tmp[,"indel_id"], function(x){ length(which(tmp[,"indel_id"] %in% x)) >= nrSamples }),]
            if(nrow(tmp2)>0){
                tmp3 <- tmp2[which(tmp2[,"sample_name"] %in% inSamples),]
                if(nrow(tmp3)>0){
                    candIndels <- tmp3[sapply(tmp3[,"indel_id"], function(x){ length(which(tmp3[,"indel_id"] %in% x)) == nrSamples }),]
                }
            }
        }
    }

    cand <- list()

    if(!is.null(candSNPs)){
        if(nrow(candSNPs) == 0)
            candSNPs <- NULL
    }
    if(!is.null(candIndels)){
        if(nrow(candIndels) == 0)
            candIndels <- NULL
    }

    cand[["hetero_snps"]] <- candSNPs
    cand[["hetero_indels"]] <- candIndels

    if(!is.na(outfile)){
        cat("*********** CanvasDB filtering analysis report ***********\n", file=outfile)
        cat("**  Analysis of dominant variants in samples: ",paste(inSamples, collapse=",",sep=""),"\n", file=outfile, append=TRUE)
        cat("**  Variants allowed to be detected in at most ",maxOthers," (",100*maxFreq,"%) of other samples\n", sep="", file=outfile, append=TRUE)
        cat("**********************************************************\n\n", file=outfile, append=TRUE)
        cat("-------------------  Heterozygous SNPs ---------------------\n", file=outfile, append=TRUE)
        if(!(is.null(cand[["hetero_snps"]]))){
            cat(colnames(cand[["hetero_snps"]]), sep="\t", file=outfile, append=TRUE)
            cat("\n", file=outfile, append=TRUE)
            write.table(cand[["hetero_snps"]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
        }
        else{
            cat("no shared heterozygous SNPs found", file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
        cat("-------------------  Heterozygous Indels ---------------------\n", file=outfile, append=TRUE)
        if(!(is.null(cand[["hetero_indels"]]))){
            cat(colnames(cand[["hetero_indels"]]), sep="\t", file=outfile, append=TRUE)
            cat("\n", file=outfile, append=TRUE)
            write.table(cand[["hetero_indels"]], quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
        }
        else{
            cat("no shared heterozygous indels found", file=outfile, append=TRUE)
        }
        cat("\n\n", file=outfile, append=TRUE)
    }

    return(cand)

}


## Function for filtering for de novo variants: Detects variants present in one sample only
filterDenovo <- function(sample, discardSamples=NULL, minAltReads=3, minSeverity=3, dbSNPfilterCommon=TRUE,  outfile=NA){

    candSNPs <- filterSNPs(sample, discardSamples=discardSamples, minSeverity=minSeverity, dbSNPfilterCommon=dbSNPfilterCommon)
    candSNPs <- candSNPs[candSNPs[,"alt_counts"]>=minAltReads,]

    candIndels <- filterIndels(sample, discardSamples=discardSamples, minSeverity=minSeverity, dbSNPfilterCommon=dbSNPfilterCommon)
    candIndels <- candIndels[as.numeric(candIndels[,"coverage"])>=minAltReads,]

    cat("*********** CanvasDB filtering analysis report ***********\n", file=outfile)
    cat("**  Analysis of de novo variants in sample: ",sample,"\n", file=outfile, append=TRUE)
    if(!is.null(discardSamples))
        cat("**  Samples excluded from filtering: ",paste(discardSamples, collapse=",",sep=""),"\n", file=outfile, append=TRUE)
    cat("**  Minimum number of alt reads: ",minAltReads,"\n", file=outfile, append=TRUE)
    cat("**********************************************************\n\n", file=outfile, append=TRUE)
    cat("-------------------  SNPs ---------------------\n", file=outfile, append=TRUE)
    if(!(is.null(candSNPs))){
        cat(colnames(candSNPs), sep="\t", file=outfile, append=TRUE)
        cat("\n", file=outfile, append=TRUE)
        write.table(candSNPs, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
    }
    else{
        cat("no candidate SNPs found", file=outfile, append=TRUE)
    }
    cat("\n\n", file=outfile, append=TRUE)
    cat("-------------------  Indels ---------------------\n", file=outfile, append=TRUE)
    if(!(is.null(candIndels))){
        cat(colnames(candIndels), sep="\t", file=outfile, append=TRUE)
        cat("\n", file=outfile, append=TRUE)
        write.table(candIndels, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", file=outfile, append=TRUE)
    }
    else{
        cat("no candidate indels found", file=outfile, append=TRUE)
    }
    cat("\n\n", file=outfile, append=TRUE)
}
