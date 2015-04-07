############################################################################
##
## File: canvasDB.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Functions for building and populating a MySQL database with
## SNP and indel data from high-throughput sequencing experinents.
##
############################################################################

rootDir <- "/Volumes/Data/canvasDB/" # Path to this directory, where R code is located

## Load functions for conftructing an canvasDB
source(paste(rootDir,"canvasDB_config.R",sep=""))
source(paste(rootDir,"canvasDB_parse.R",sep=""))
source(paste(rootDir,"canvasDB_annotate.R",sep=""))
source(paste(rootDir,"canvasDB_setup.R",sep=""))
source(paste(rootDir,"canvasDB_filter.R",sep=""))

## Load plugins (if available)
if(file.exists(paste(rootDir,"plugins/canvasDB_homozmap.R",sep="")))
    source(paste(rootDir,"plugins/canvasDB_homozmap.R",sep=""))

## Sub-function to addVariantData. This function is called with an object containing SNPs to be updated
updateSNPsummaryTable <- function(SNPsToBeUpdated, SNPsToBeAdded, TALK=FALSE){

    SNPsummary.table <- dbTables[["SNP.summary"]]

    ## ##############################################
    ##
    ## Step1: Update sample info in SNPsummary table
    ##
    ## ##############################################

    if(nrow(SNPsToBeUpdated)>0){

        cat(" Updating SNP summary table\n")

        cat(" Step 1: Updating SNPs already in summary table...")
        ps <- proc.time()[3]

        tmpDatafile <- paste(tmpfileDir,"tmp_SNPsToBeUpdated.txt",sep="")
        SNPidsToBeUpdated <- SNPsToBeUpdated[,"SNP_id"]
        write.table(SNPidsToBeUpdated, file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

        SNPsToBeUpdated.table <- "tmp_SNPsToBeUpdated"
        con <- connectToInhouseDB()

        query <- paste("DROP TABLE IF EXISTS ",SNPsToBeUpdated.table, sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        query <- paste("CREATE TABLE ",SNPsToBeUpdated.table,
                       "(SNP_id varchar(20), ",
                       "PRIMARY KEY(SNP_id)",
                       ") ENGINE=MyISAM DEFAULT CHARSET=latin1",sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",SNPsToBeUpdated.table,";",sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        ## Query for fetching data from summary data for SNPs to be updated
        query <- paste("SELECT * FROM ",SNPsummary.table," AS t1, ",SNPsToBeUpdated.table," AS t2 WHERE t1.SNP_id=t2.SNP_id;")

        SNPsToUpdate <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)

        file.remove(tmpDatafile)

        if(nrow(SNPsToUpdate)>0){
            SNPsamplesToBeUpdated <- SNPsToBeUpdated[,"sample_str"]
            names(SNPsamplesToBeUpdated) <- SNPidsToBeUpdated

            SNPnrSamplesToBeUpdated <- as.numeric(SNPsToBeUpdated[,"nr_samples"])
            names(SNPnrSamplesToBeUpdated) <- SNPidsToBeUpdated

            SNPsToUpdate[,"nr_samples"] <- as.numeric(SNPsToUpdate[,"nr_samples"])+SNPnrSamplesToBeUpdated[SNPsToUpdate[,"SNP_id"]]
            SNPsToUpdate[,"samples"] <- paste(SNPsToUpdate[,"samples"],",",SNPsamplesToBeUpdated[SNPsToUpdate[,"SNP_id"]],sep="")

            tmpDatafile <- paste(tmpfileDir,"tmp_SNPsToBeReplacedInSummary.txt",sep="")

            write.table(SNPsToUpdate, file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

            con <- connectToInhouseDB()

            #query <- paste("DELETE t1 FROM ",SNPsummary.table," as t1 INNER JOIN ",SNPsToBeUpdated.table," AS t2 ON t1.SNP_id=t2.SNP_id;",sep="")
            #tmp <- dbGetQuery_E(con, query, TALK=TALK)

            #query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",SNPsummary.table,";",sep="")
            #tmp <- dbGetQuery_E(con, query, TALK=TALK)

            query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' REPLACE INTO TABLE ",SNPsummary.table,";",sep="")
            tmp <- dbGetQuery_E(con, query, TALK=TALK)

            dbDisconnect(con)

            file.remove(tmpDatafile)
        }

        cat(proc.time()[3] - ps,"s\n");

    }



    ## ############################################
    ##
    ##  Step2: Add new SNPs to SNP summary table
    ##
    ## ############################################

    cat(" Step 2: Adding new SNPs to summary table..\n")

    nrSNPsToBeAdded <- nrow(SNPsToBeAdded)

    if(nrSNPsToBeAdded>0){
        cat(" ",nrow(SNPsToBeAdded)," new entries to be added...\n",sep="")

        dataToBeAdded <- annotateSNPs(SNPsToBeAdded, tmpAnnotationDir=tmpfileDir, dbSNPversion=dbSNPversion)

        cat("  - Loading annotated SNPs into summary table...")
        ps <- proc.time()[3]

        tmpDatafile <- paste(tmpfileDir,"tmp_SNPsToBeAdded.txt",sep="")

        write.table(dataToBeAdded, file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

        con <- connectToInhouseDB()

        query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",SNPsummary.table,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        file.remove(tmpDatafile)

        cat(proc.time()[3] - ps,"s\n");

        dbDisconnect(con)
    }

}



## Sub-function to addVariantData. This function is called with an object containing indels to be updated
updateIndelSummaryTable <- function(indelsToBeUpdated, indelsToBeAdded, TALK=FALSE){

    indelSummary.table <- dbTables[["indel.summary"]]

    ## ##############################################
    ##
    ## Step1: Update sample info in indelSummary table
    ##
    ## ##############################################

    if(nrow(indelsToBeUpdated)>0){

        cat(" Updating indel summary table\n")

        cat(" Step 1: Updating indels already in summary table...")
        ps <- proc.time()[3]

        tmpDatafile <- paste(tmpfileDir,"tmp_indelsToBeUpdated.txt",sep="")
        indelIdsToBeUpdated <- indelsToBeUpdated[,"indel_id"]
        write.table(indelIdsToBeUpdated, file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

        indelsToBeUpdated.table <- "tmp_indelsToBeUpdated"
        con <- connectToInhouseDB()

        query <- paste("DROP TABLE IF EXISTS ",indelsToBeUpdated.table, sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        query <- paste("CREATE TABLE ",indelsToBeUpdated.table,
                       "(indel_id varchar(500), ",
                       "PRIMARY KEY(indel_id)",
                       ") ENGINE=MyISAM DEFAULT CHARSET=latin1",sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",indelsToBeUpdated.table,";",sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        ## Query for fetching data from summary data for indels to be updated
        query <- paste("SELECT * FROM ",indelSummary.table," AS t1, ",indelsToBeUpdated.table," AS t2 WHERE t1.indel_id=t2.indel_id;")

        indelsToUpdate <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)

        file.remove(tmpDatafile)

        if(nrow(indelsToUpdate)>0){
            indelSamplesToBeUpdated <- indelsToBeUpdated[,"sample_str"]
            names(indelSamplesToBeUpdated) <- indelIdsToBeUpdated

            indelNrSamplesToBeUpdated <- as.numeric(indelsToBeUpdated[,"nr_samples"])
            names(indelNrSamplesToBeUpdated) <- indelIdsToBeUpdated

            indelsToUpdate[,"nr_samples"] <- as.numeric(indelsToUpdate[,"nr_samples"])+indelNrSamplesToBeUpdated[indelsToUpdate[,"indel_id"]]
            indelsToUpdate[,"samples"] <- paste(indelsToUpdate[,"samples"],",",indelSamplesToBeUpdated[indelsToUpdate[,"indel_id"]],sep="")

            tmpDatafile <- paste(tmpfileDir,"tmp_indelsToBeReplacedInSummary.txt",sep="")

            write.table(indelsToUpdate, file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

            con <- connectToInhouseDB()

            #query <- paste("DELETE t1 FROM ",indelSummary.table," as t1 INNER JOIN ",indelsToBeUpdated.table," AS t2 ON t1.indel_id=t2.indel_id;",sep="")
            #tmp <- dbGetQuery_E(con, query, TALK=TALK)

            #query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",indelSummary.table,";",sep="")
            #tmp <- dbGetQuery_E(con, query, TALK=TALK)

            query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' REPLACE INTO TABLE ",indelSummary.table,";",sep="")
            tmp <- dbGetQuery_E(con, query, TALK=TALK)

            dbDisconnect(con)

            file.remove(tmpDatafile)
        }

        cat(proc.time()[3] - ps,"s\n");
    }



    ## ############################################
    ##
    ##  Step2: Add new Indels to indel summary table
    ##
    ## ############################################

    cat(" Step 2: Adding new indels to summary table..\n")

    nrIndelsToBeAdded <- nrow(indelsToBeAdded)

    if(nrIndelsToBeAdded>0){
        cat(" ",nrow(indelsToBeAdded)," new entries to be added...\n",sep="")

        dataToBeAdded <- annotateIndels(indelsToBeAdded, tmpAnnotationDir=tmpfileDir, dbSNPversion=dbSNPversion)

        cat("  - Loading annotated indels into summary table...")
        ps <- proc.time()[3]

        tmpDatafile <- paste(tmpfileDir,"tmp_indelsToBeAdded.txt",sep="")

        write.table(dataToBeAdded, file=tmpDatafile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

        con <- connectToInhouseDB()

        query <- paste("LOAD DATA LOCAL INFILE '",tmpDatafile,"' INTO TABLE ",indelSummary.table,";",sep="")

        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        dbDisconnect(con)

        file.remove(tmpDatafile)

        cat(proc.time()[3] - ps,"s\n");
    }

}



## Function for adding SNP and indel data to the canvasDB. Samples are added one by one.
addVariantData <- function(SNP.file, canvasId, seq.platform, library.type, read.type, fileFormat, indel.file="", capture.method="", instrument.name="unknown", sampleName="", gender="", geographic.location="", phenotypes="", comments="", principal.investigator="", totalReads="", readsOnTarget="", date="", TALK=FALSE){

    cat("Processing sample: ",canvasId,"\n",sep="")

    ## Supported platforms
    if(!(seq.platform %in% c("SOLiD4", "SOLiD5500", "SOLiD5500W", "Proton", "Illumina", ""))){
        stop(paste("Unsupported sequencing platform:",seq.platform))
    }

    ## Supported library types
    if(!(library.type %in% c("WholeExome","WholeTranscriptome","WholeGenome"))){
        stop(paste("Unsupported library type:",library.type))
    }

    ## Supported SNP calling formats
    if(!(fileFormat %in% c("Lifescope", "Bioscope", "TorrentSuite", "CanvasDB", "GATK"))){
        stop(paste("Unsupported SNP format:",fileFormat))
    }


    ###########################
    ##
    ## Add run and sample data
    ##
    ###########################

    con <- connectToInhouseDB()

    ## Add sample data
    sample.table <- dbTables[["sample"]]

    ## Check if sample already exists...
    query <- paste("SELECT sample_id FROM ",sample.table," WHERE canvas_id='",canvasId,"';",sep="")
    sampleId <- dbGetQuery_E(con,query,TALK=TALK)

    if(nrow(sampleId) > 0){
        warning(paste("Sample '",canvasId,"' already exists in db. Run/sample data cannot be added.\n",sep=""))
        dbDisconnect(con)
        stop()
    }

    query <- paste("INSERT INTO ",sample.table," (canvas_id, gender, sample_name, geographic_location, phenotypes, comments, principal_investigator) VALUES ('",canvasId,"','",gender,"','",sampleName,"','",geographic.location,"','",phenotypes,"','",comments,"','",principal.investigator,"');",sep="")

    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    ## Get sample id
    query <- paste("SELECT sample_id FROM ",sample.table," WHERE canvas_id='",canvasId,"';",sep="")
    sampleId <- dbGetQuery_E(con,query,TALK=TALK)

    ## Add run data
    run.table <- dbTables[["run"]]
    query <- paste("REPLACE INTO ",run.table," (sample_id, canvas_id, seq_platform, library_type, read_type, capture_method, instrument_name, file_format, total_reads, reads_on_target, run_start) VALUES ('",sampleId,"','",canvasId,"','",seq.platform,"','",library.type,"','",read.type,"','",capture.method,"','",instrument.name,"','",fileFormat,"','",totalReads,"','",readsOnTarget,"','",date,"');",sep="")

    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    dbDisconnect(con)

    variantsForSample <- list()
    variantsForSample[["sample_id"]] <- sampleId
    variantsForSample[["SNPs"]] <- NULL
    variantsForSample[["indels"]] <- NULL

    ## #####################################
    ##
    ## Create SNP db table and add SNP data
    ##
    ## #####################################

    SNPdata <- NULL

    ## As default SNP.file is assumed to be on indb SNP format. Other supported formats are parsed below.
    indbFormatedSNPfile <- SNP.file

    ## Lifescope/Bioscope format
    if(fileFormat %in% c("Lifescope","Bioscope")){
        cat("Reading/parsing SNP data file...")
        ps <- proc.time()[3]

        SNPdata <- parseLifescopeSNPs(SNP.file)
        if(fileFormat == "Bioscope"){
            SNPdata[,"chr"] <- paste("chr",SNPdata[,"chr"],sep="")
        }

        SNPids <- paste(SNPdata[,"chr"],SNPdata[,"pos"],SNPdata[,"ref"],SNPdata[,"alt"],sep="|")

        tmp <- cbind(SNPids, SNPdata)

        indbFormatedSNPfile <- paste(tmpfileDir,"tmp_indb_snps.txt",sep="")

        write.table(tmp,file=indbFormatedSNPfile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

        cat(proc.time()[3] - ps,"s\n");
    }

    ## Torrent Suite format (for PGM/Proton)
    if(fileFormat %in% c("TorrentSuite")){
        cat("Reading/parsing SNP data file...\n")
        ps <- proc.time()[3]

        SNPdata <- parseSNPsFromIonTorrentVCFs(SNP.file)

        SNPids <- paste(SNPdata[,"chr"],SNPdata[,"pos"],SNPdata[,"ref"],SNPdata[,"alt"],sep="|")

        tmp <- cbind(SNPids, SNPdata)

        indbFormatedSNPfile <- paste(tmpfileDir,"tmp_indb_snps.txt",sep="")

        write.table(tmp,file=indbFormatedSNPfile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

        cat(proc.time()[3] - ps,"s\n");
    }

    ## GATK format (for Illumina)
    if(fileFormat %in% c("GATK")){
        cat("Reading/parsing SNP data file...\n")
        ps <- proc.time()[3]

        SNPdata <- parseSNPsFromIlluminaVCFs(SNP.file)

        SNPids <- paste(SNPdata[,"chr"],SNPdata[,"pos"],SNPdata[,"ref"],SNPdata[,"alt"],sep="|")

        tmp <- cbind(SNPids, SNPdata)

        indbFormatedSNPfile <- paste(tmpfileDir,"tmp_indb_snps.txt",sep="")

        write.table(tmp,file=indbFormatedSNPfile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

        cat(proc.time()[3] - ps,"s\n");
    }


    cat("Loading SNP data into database...")
    ps <- proc.time()[3]

    con <- connectToInhouseDB()
    SNPdata.table <- paste("snp_data_",sampleId,sep="")

    query <- paste("DROP TABLE IF EXISTS ",SNPdata.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",SNPdata.table,
                   "(SNP_id varchar(20), ",
                   "chr varchar(10), ",
                   "pos integer(11), ",
                   "ref varchar(1), ",
                   "alt varchar(1), ",
                   "coverage integer(10), ",
                   "ref_counts integer(10), ",
                   "ref_starts integer(10), ",
                   "ref_mean_qv integer(3), ",
                   "alt_counts integer(10), ",
                   "alt_starts integer(10), ",
                   "alt_mean_qv integer(3), ",
                   "heterozygous integer(1), ",
                   "PRIMARY KEY(SNP_id)",
                   ") ENGINE=MyISAM DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    colsFromFile <- "(SNP_id,chr,pos,ref,alt,coverage,ref_counts,ref_starts,ref_mean_qv,alt_counts,alt_starts,alt_mean_qv,heterozygous)"
    query <- paste("LOAD DATA LOCAL INFILE '",indbFormatedSNPfile,"' INTO TABLE ",SNPdata.table," ",colsFromFile,";",sep="")
    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    ## Remove temporary files...
    if(fileFormat %in% c("Lifescope","Bioscope","TorrentSuite")){
      file.remove(indbFormatedSNPfile)
    }

    query <- paste("SELECT SNP_id FROM ",SNPdata.table,";",sep="")
    snpIds <- dbGetQuery_E(con,query,TALK=TALK)

    variantsForSample[["SNPs"]] <- "tjena"
    variantsForSample[["SNPs"]] <- snpIds[,1]

    dbDisconnect(con)

    cat(proc.time()[3] - ps,"s\n");


    ## ###############################
    ##
    ## Add indel data (if available)
    ##
    ## ###############################


    if(!(indel.file=="")){

        indbFormatedIndelFile <- indel.file

        indelData <- NULL

        ## Lifescope format.
        if(fileFormat %in% c("Lifescope")){
            cat("Reading/parsing indel data file...")
            ps <- proc.time()[3]

            indelData <- parseLifescopeIndels(indel.file)

            indelIds <- array(NA, nrow(indelData))

            for(i in 1:nrow(indelData)){
                chr <- as.character(indelData[i,"chr"])
                start <- as.character(as.numeric(indelData[i,"start"]))
                end <- as.character(as.numeric(indelData[i,"end"]))
                ref <- as.character(indelData[i,"ref"])
                alt <- as.character(indelData[i,"alt"])

                indelIds[i] <- paste(chr,start,end,ref,alt,sep="|")
            }

            tmp <- cbind(indelIds, indelData)

            indbFormatedIndelFile <- paste(tmpfileDir,"tmp_indb_indels.txt",sep="")

            write.table(tmp,file=indbFormatedIndelFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

            cat(proc.time()[3] - ps,"s\n");

        }

        if(fileFormat %in% c("TorrentSuite")){
            cat("Reading/parsing Indel data file...")
            ps <- proc.time()[3]

            indelData <- parseIndelsFromIonTorrentVCFs(indel.file)

            indelIds <- paste(indelData[,"chr"],indelData[,"start"],indelData[,"end"],indelData[,"ref"],indelData[,"alt"],sep="|")

            tmp <- cbind(indelIds, indelData)

            indbFormatedIndelFile <- paste(tmpfileDir,"tmp_indb_indels.txt",sep="")

            write.table(tmp,file=indbFormatedIndelFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

            cat(proc.time()[3] - ps,"s\n");
        }

        if(fileFormat %in% c("GATK")){
            cat("Reading/parsing Indel data file...")
            ps <- proc.time()[3]

            indelData <- parseIndelsFromIlluminaVCFs(indel.file)

            indelIds <- paste(indelData[,"chr"],indelData[,"start"],indelData[,"end"],indelData[,"ref"],indelData[,"alt"],sep="|")

            tmp <- cbind(indelIds, indelData)

            indbFormatedIndelFile <- paste(tmpfileDir,"tmp_indb_indels.txt",sep="")

            write.table(tmp,file=indbFormatedIndelFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

            cat(proc.time()[3] - ps,"s\n");
        }


        cat("Loading indel data into database...")
        ps <- proc.time()[3]

        con <- connectToInhouseDB()

        indelData.table <- paste("indel_data_",sampleId,sep="")

        query <- paste("DROP TABLE IF EXISTS ",indelData.table, sep="")
        tmp <- dbGetQuery_E(con, query, TALK=TALK)

        query <- paste("CREATE TABLE ",indelData.table,
                       "(indel_id varchar(500), ",
                       "chr varchar(10), ",
                       "start integer(11), ",
                       "end integer(11), ",
                       "ref varchar(100), ",
                       "alt varchar(100), ",
                       "coverage integer(10), ",
                       "ref_counts integer(10), ",
                       "alt_counts integer(10), ",
                       "heterozygous integer(1), ",
                       "PRIMARY KEY(indel_id)",
                       ") ENGINE=MyISAM DEFAULT CHARSET=latin1",sep="")

        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        colsFromFile <- "(indel_id,chr,start,end,ref,alt,coverage,ref_counts,alt_counts,heterozygous)"
        query <- paste("LOAD DATA LOCAL INFILE '",indbFormatedIndelFile,"' INTO TABLE ",indelData.table," ",colsFromFile,";",sep="")
        tmp <- dbGetQuery_E(con,query,TALK=TALK)

        ## Remove temporary files...
        if(fileFormat %in% c("Lifescope","Bioscope","TorrentSuite")){
            file.remove(indbFormatedIndelFile)
        }

        query <- paste("SELECT indel_id FROM ",indelData.table,";",sep="")
        indelIds <- dbGetQuery_E(con,query,TALK=TALK)

        variantsForSample[["indels"]] <- indelIds[,1]

        dbDisconnect(con)

        cat(proc.time()[3] - ps,"s\n");

    }


    return(variantsForSample)

}



## Function for adding a set of samples listed in 'infile' to the canvasDB
batchImport <- function(infile="dummy_exomes.txt", nlines=NA, TALK=FALSE){

    ## Data structures for keeping track of SNPs and indels to be updated
    SNPsToBeInserted <- list()
    SNPsToBeInserted[["nrSamples"]] <- array(NA,0)
    SNPsToBeInserted[["sampleStr"]] <- array(NA,0)

    indelsToBeInserted <- list()
    indelsToBeInserted[["nrSamples"]] <- array(NA,0)
    indelsToBeInserted[["sampleStr"]] <- array(NA,0)

    ## Data structures for keeping track of SNPs and indels in summary tables
    SNPsummary.table <- dbTables[["SNP.summary"]]
    indelSummary.table <- dbTables[["indel.summary"]]

    SNPidsInSummary <- NULL
    indelIdsInSummary <- NULL

    con <- connectToInhouseDB()

    query <- paste("SELECT SNP_id FROM ",SNPsummary.table,";")
    currSNPids <- dbGetQuery_E(con,query,TALK=TALK)

    if(nrow(currSNPids)>0){
        SNPidsInSummary <- currSNPids[,"SNP_id"]
    }

    query <- paste("SELECT indel_id FROM ",indelSummary.table,";")
    currIndelIds <- dbGetQuery_E(con,query,TALK=TALK)
    if(nrow(currIndelIds)>0){
        indelIdsInSummary <- currIndelIds[,"indel_id"]
    }

    dbDisconnect(con)

    ## Get canvas_ids for samples already in database
    con <- connectToInhouseDB()

    sample.table <- dbTables[["sample"]]

    ## Query for fetching sample_ids that are already in summary
    query <- paste("SELECT canvas_id FROM ",sample.table,";",sep="")

    res <- dbGetQuery_E(con,query,TALK=TALK)

    dbDisconnect(con)

    canvasIds <- NULL

    if(nrow(res)>0){
        canvasIds <- as.character(res[,1])
    }

    cat("*** Batch import of variant data from: ",infile," ***\n",sep="")

    data <- read.table(infile,sep="|",header=TRUE, as.is=TRUE)
	## make sure that there are valid numerical numbers in the read count columns. 
	data[,"reads.total"][is.na(data[,"reads.total"])]  <- -1  
	data[,"reads.on.target"][is.na(data[,"reads.on.target"])]  <- -1  	  


    ## Don't add any sample_ids already in the database
    if(!is.null(canvasIds)){
        data <- data[which(!(data[,"canvasId"] %in% canvasIds)),]
    }

    if(is.na(nlines))
        nlines <- nrow(data)

    if(!(nrow(data)>0)){
        stop("No new samples to be added!")
    }

    updateList <- seq(1,nlines,100)
    chunkSize <- 500000

    ## Go through each sample and import data
    for(i in 1:nlines){

        cat("Processing sample ",i," of ",nlines,"...\n",sep="")

        data[i,][is.na(data[i,])] <- ""

        canvasId <- data[i,"canvasId"]
        sampleName <- data[i,"sampleName"]
        seqPlatform <- data[i,"seq.platform"]
        libraryType <- data[i,"library.type"]
        readType <- data[i,"read.type"]
        captureMethod <- data[i,"capture.method"]
        date <- data[i,"date"]
        captureMethod <- data[i,"capture.method"]
        pi <- data[i,"principal.investigator"]
        instrumentName <- data[i,"instrument.name"]
        gender <- data[i,"gender"]
        geoLocation <- data[i,"geographic.location"]
        phenotypes <- data[i,"phenotypes"]
        comments <- data[i,"comments"]
        SNPfile <- data[i,"SNP.file"]
        indelFile <- data[i,"indel.file"]
        fileFormat <- data[i,"file.format"]
        totalReads <- data[i,"reads.total"]
        readsOnTarget <- data[i,"reads.on.target"]

        ## Add data only for samples not already in database..
        if(canvasId %in% canvasIds){
            cat("Failed to import data for ",canvasId,". Already in database!!\n",sep="")
        }
        else{
            if(file.exists(SNPfile)){
                if((indelFile=="") || file.exists(indelFile)){
                    variantsForSample <- addVariantData(SNPfile, canvasId, seqPlatform, libraryType, readType, fileFormat, indel.file=indelFile, capture.method=captureMethod, instrument.name=instrumentName, sampleName=sampleName, gender=gender, geographic.location=geoLocation, phenotypes=phenotypes, comments=comments, principal.investigator=pi, totalReads=totalReads, readsOnTarget=readsOnTarget, date=date, TALK=TALK)

                    cat("Updating SNP/indel datastructures...")
                    ps <- proc.time()[3]

                    sampleId <- variantsForSample[["sample_id"]]
                    SNPids <- variantsForSample[["SNPs"]]
                    indelIds <- variantsForSample[["indels"]]

                    sampleStr <- paste(",",sampleId,sep="")

                    ## For already seen SNPs, update datastructure with latest sample
                    seenSNPs <- (names(SNPsToBeInserted[["nrSamples"]]) %in% SNPids)
                    SNPsToBeInserted[["sampleStr"]][seenSNPs] <- paste(SNPsToBeInserted[["sampleStr"]][seenSNPs],sampleStr,sep="")
                    SNPsToBeInserted[["nrSamples"]][seenSNPs] <- SNPsToBeInserted[["nrSamples"]][seenSNPs]+1

                    ## Add new SNPs to datastructure
                    newSNPids <- SNPids[!(SNPids %in% names(SNPsToBeInserted[["nrSamples"]]))]

                    sampleStrsToBeAdded <- rep(as.character(sampleId), length(newSNPids))
                    names(sampleStrsToBeAdded) <- newSNPids
                    nrSamplesToBeAdded <- rep(1, length(newSNPids))
                    names(nrSamplesToBeAdded) <- newSNPids
                    SNPsToBeInserted[["sampleStr"]] <- c(SNPsToBeInserted[["sampleStr"]], sampleStrsToBeAdded)
                    SNPsToBeInserted[["nrSamples"]] <- c(SNPsToBeInserted[["nrSamples"]], nrSamplesToBeAdded)

                    ## For already seen indels, update datastructure with latest sample
                    seenIndels <- (names(indelsToBeInserted[["nrSamples"]]) %in% indelIds)
                    indelsToBeInserted[["sampleStr"]][seenIndels] <- paste(indelsToBeInserted[["sampleStr"]][seenIndels],sampleStr,sep="")
                    indelsToBeInserted[["nrSamples"]][seenIndels] <- indelsToBeInserted[["nrSamples"]][seenIndels]+1

                    ## Add new indels to datastructure
                    newIndelIds <- indelIds[!(indelIds %in% names(indelsToBeInserted[["nrSamples"]]))]

                    sampleStrsToBeAdded <- rep(as.character(sampleId), length(newIndelIds))
                    names(sampleStrsToBeAdded) <- newIndelIds
                    nrSamplesToBeAdded <- rep(1, length(newIndelIds))
                    names(nrSamplesToBeAdded) <- newIndelIds
                    indelsToBeInserted[["sampleStr"]] <- c(indelsToBeInserted[["sampleStr"]], sampleStrsToBeAdded)
                    indelsToBeInserted[["nrSamples"]] <- c(indelsToBeInserted[["nrSamples"]], nrSamplesToBeAdded)

                    SNPobjectsize <- object.size(SNPsToBeInserted)
                    indelObjectsize <- object.size(indelsToBeInserted)

                    cat(proc.time()[3] - ps,"s\n");
                    cat("Object sizes, SNPs:",SNPobjectsize," indels:",indelObjectsize,"\n")

                    if((i %in% updateList) || (i == nlines)){ ## This should be done at some intervals..
                        cat(" Preparing for SNP summary update...")
                        ps <- proc.time()[3]

                        SNPsToBeInsertedTable <- cbind(names(SNPsToBeInserted[["nrSamples"]]), SNPsToBeInserted[["nrSamples"]], SNPsToBeInserted[["sampleStr"]])
                        colnames(SNPsToBeInsertedTable) <- c("SNP_id","nr_samples","sample_str")

                        totalNrToInsert <- nrow(SNPsToBeInsertedTable)

                        ## Divide large data into smaller chunks, to make sure MySQL doesn't crash
                        chunkStart <- 1
                        chunkEnd <- min(chunkSize,totalNrToInsert)
                        moreChunksToInsert <- TRUE

                        while(moreChunksToInsert){

                            SNPsToBeInsertedTableChunk <- SNPsToBeInsertedTable[chunkStart:chunkEnd,]

                            alreadyInSummary <- (SNPsToBeInsertedTableChunk[,"SNP_id"] %in% SNPidsInSummary)
                            SNPsToBeUpdated <- SNPsToBeInsertedTableChunk[alreadyInSummary,,drop=FALSE]
                            SNPsToBeAdded <- SNPsToBeInsertedTableChunk[!alreadyInSummary,,drop=FALSE]
                            SNPidsInSummary <- c(SNPidsInSummary,SNPsToBeAdded[,"SNP_id"])

                            cat(" chunk:",chunkStart,"-",chunkEnd," total:",nrow(SNPsToBeInsertedTable)," to update:",nrow(SNPsToBeUpdated)," to add:",nrow(SNPsToBeAdded),proc.time()[3] - ps,"s\n");
                            updateSNPsummaryTable(SNPsToBeUpdated, SNPsToBeAdded)

                            if(chunkEnd<totalNrToInsert){
                                chunkStart <- chunkEnd+1
                                chunkEnd <-  min(chunkEnd+chunkSize,totalNrToInsert)
                            }
                            else{
                                moreChunksToInsert <- FALSE
                            }
                        }

                        ## Reset SNPs to be updated. Need to release memory??
                        SNPsToBeInserted <- list()
                        SNPsToBeInserted[["nrSamples"]] <- array(NA,0)
                        SNPsToBeInserted[["sampleStr"]] <- array(NA,0)

                        if(length(indelsToBeInserted$nrSamples)>0){
                            cat(" Preparing for indel summary update...")
                            ps <- proc.time()[3]

                            indelsToBeInsertedTable <- cbind(names(indelsToBeInserted[["nrSamples"]]), indelsToBeInserted[["nrSamples"]], indelsToBeInserted[["sampleStr"]])
                            colnames(indelsToBeInsertedTable) <- c("indel_id","nr_samples","sample_str")

                            totalNrToInsert <- nrow(indelsToBeInsertedTable)

                            ## Divide large data into smaller chunks, to make sure MySQL doesn't crash
                            chunkStart <- 1
                            chunkEnd <- min(chunkSize,totalNrToInsert)
                            moreChunksToInsert <- TRUE

                            while(moreChunksToInsert){

                                 indelsToBeInsertedTableChunk <- indelsToBeInsertedTable[chunkStart:chunkEnd,]

                                 alreadyInSummary <- (indelsToBeInsertedTableChunk[,"indel_id"] %in% indelIdsInSummary)
                                 indelsToBeUpdated <- indelsToBeInsertedTableChunk[alreadyInSummary,,drop=FALSE]
                                 indelsToBeAdded <- indelsToBeInsertedTableChunk[!alreadyInSummary,,drop=FALSE]
                                 indelIdsInSummary <- c(indelIdsInSummary,indelsToBeAdded[,"indel_id"])

                                 cat(proc.time()[3] - ps,"s\n");
                                 updateIndelSummaryTable(indelsToBeUpdated, indelsToBeAdded)

                                 if(chunkEnd<totalNrToInsert){
                                     chunkStart <- chunkEnd+1
                                     chunkEnd <-  min(chunkEnd+chunkSize,totalNrToInsert)
                                 }
                                 else{
                                     moreChunksToInsert <- FALSE
                                 }

                             }

                             ## Reset indels to be updated. Need to release memory??
                             indelsToBeInserted <- list()
                             indelsToBeInserted[["nrSamples"]] <- array(NA,0)
                             indelsToBeInserted[["sampleStr"]] <- array(NA,0)
                         }

                    }
                }
                else{
                    cat("Failed to import data. Indel data file does not exist: ",indelFile,"\n",sep="")
                }
            }
            else{
                cat("Failed to import data. SNP data file does not exist: ",SNPfile,"\n",sep="")
            }
        }

    }
}




































