############################################################################
##
## File: canvasDB_setup.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Functions for building and populating a MySQL database with
## SNP and indel data from high-throughput sequencing experinents.
##
############################################################################

## Easy way to set up the MySQL connection.
## file - a .my.cnf file
## group - an appropriate group within the .my.cnf file
getDBCon_file <- function(file, group)
{
  require(RMySQL)
  drv <- dbDriver("MySQL")
  dbConnect(drv, default.file = file, group = group)
}


## Error controlled version of dbGetQuery
dbGetQuery_E <- function(con,q,TALK = FALSE)
{
    if(TALK)
        cat(q,"\n")
    err <- try(tmp <- dbGetQuery(con, q) ,silent = FALSE)
    if(inherits(err,"try-error"))
    {
        dbDisconnect(con)
        stop("Failed to execute SQL-query")
    }
    if(is.null(tmp))
        tmp <- data.frame()
    return(tmp)
}


## Returns a db connection
connectToInhouseDB <- function(db.name=indbDatabaseName){
    mysqlfile <- paste(rootDir,mysqlConfigFile,sep="")
    con <- getDBCon_file(mysqlfile, db.name)
    return(con)
}


## Creates an empty db
initDB <- function(TALK=FALSE){

    con <- connectToInhouseDB()

    ## Create run info table
    run.table <- dbTables[["run"]]

    query <- paste("DROP TABLE IF EXISTS ",run.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",run.table,
                   "(sample_id int(10), ",
                   "canvas_id varchar(20), ",
                   "library_type varchar(20), ",
                   "read_type varchar(20), ",
                   "capture_method varchar(20), ",
                   "seq_platform varchar(20), ",
                   "instrument_name varchar(20), ",
                   "file_format varchar(20), ",
                   "total_reads int(20), ",
                   "reads_on_target int(20), ",
                   "run_start date, ",
                   "PRIMARY KEY (sample_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con,query,TALK=TALK)


    ## Create sample table
    sample.table <- dbTables[["sample"]]

    query <- paste("DROP TABLE IF EXISTS ",sample.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",sample.table,
                   "(sample_id int(10) NOT NULL AUTO_INCREMENT , ",
                   "canvas_id varchar(20), ",
                   "sample_name varchar(50), ",
                   "gender varchar(10), ",
                   "geographic_location varchar(50), ",
                   "phenotypes varchar(200), ",
                   "comments varchar(50), ",
                   "principal_investigator varchar(10), ",
                   "PRIMARY KEY (sample_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    ## Create SNP summary table
    SNPsummary.table <- dbTables[["SNP.summary"]]
    dbSNP.table <- dbTables[["dbSNP"]]
    dbSNPcommon.table <- dbTables[["dbSNPcommon"]]

    query <- paste("DROP TABLE IF EXISTS ",SNPsummary.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",SNPsummary.table,
                   "(SNP_id varchar(20), ",
                   "chr varchar(10), ",
                   "pos integer(11), ",
                   "ref varchar(1), ",
                   "alt varchar(1), ",
                   "nr_samples integer(6),",
                   "samples blob,",
                   "snp",dbSNPversion," varchar(20), ",
                   "snp",dbSNPversion,"common varchar(20), ",
                   "class varchar(40), ",
                   "severity integer(2), ",
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

    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    query <- paste("CREATE INDEX summary_idx ON ",SNPsummary.table," (nr_samples, snp",dbSNPversion,"common,  snp",dbSNPversion,", severity);",sep="")
    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    query <- paste("CREATE INDEX chrpos ON ",SNPsummary.table," (chr, pos);",sep="")
    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    indelSummary.table <- dbTables[["indel.summary"]]

    query <- paste("DROP TABLE IF EXISTS ",indelSummary.table, sep="")
    tmp <- dbGetQuery_E(con, query, TALK=TALK)

    query <- paste("CREATE TABLE ",indelSummary.table,
                   "(indel_id varchar(500), ",
                   "chr varchar(10), ",
                   "start integer(11), ",
                   "end integer(11), ",
                   "ref varchar(200), ",
                   "alt varchar(200), ",
                   "type varchar(100), ",
                   "size integer(10), ",
                   "nr_samples integer(6),",
                   "samples blob,",
                   "snp",dbSNPversion," varchar(20), ",
                   "snp",dbSNPversion,"common varchar(20), ",
                   "class varchar(40), ",
                   "severity integer(2), ",
                   "gene varchar(100), ",
                   "details varchar(200), ",
                   "PRIMARY KEY(indel_id)",
                   ") ENGINE=",mysqlEngine," DEFAULT CHARSET=latin1",sep="")

    tmp <- dbGetQuery_E(con,query,TALK=TALK)


    query <- paste("CREATE INDEX summary_idx ON ",indelSummary.table," (nr_samples, snp",dbSNPversion,"common,  snp",dbSNPversion,", severity);",sep="")
    tmp <- dbGetQuery_E(con,query,TALK=TALK)

    query <- paste("CREATE INDEX chrpos ON ",indelSummary.table," (chr, start, end);",sep="")
    tmp <- dbGetQuery_E(con,query,TALK=TALK)


    ## Drop SNP data tables
    SNPtablesExists <- TRUE

    while(SNPtablesExists){
        query <- paste("SELECT GROUP_CONCAT(table_name) FROM information_schema.tables WHERE table_schema = '",indbDatabaseName,"' AND table_name like 'snp_data_%';",sep="")
        tmp <- dbGetQuery_E(con,query,TALK=TALK)
        SNPdataTables <- strsplit(tmp[,1],",")[[1]]

        if(length(SNPdataTables)==1){
            if(is.na(SNPdataTables)){
                SNPtablesExists <- FALSE
            }
        }

        if(SNPtablesExists){
            for(SNPdataTable in SNPdataTables){
                print(SNPdataTable)
                query <- paste("DROP TABLE IF EXISTS ",SNPdataTable,sep="")
                tmp <- dbGetQuery_E(con,query, TALK=TALK)
            }
        }
    }

    ## Drop indel data tables
    indelTablesExists <- TRUE

    while(indelTablesExists){
        query <- paste("SELECT GROUP_CONCAT(table_name) FROM information_schema.tables WHERE table_schema = '",indbDatabaseName,"' AND table_name like 'indel_data_%';",sep="")
        tmp <- dbGetQuery_E(con,query,TALK=TALK)
        indelDataTables <- strsplit(tmp[,1],",")[[1]]

        if(length(indelDataTables)==1){
            if(is.na(indelDataTables)){
                indelTablesExists <- FALSE
            }
        }

        if(indelTablesExists){
            for(indelDataTable in indelDataTables){
                print(indelDataTable)
                query <- paste("DROP TABLE IF EXISTS ",indelDataTable,sep="")
                tmp <- dbGetQuery_E(con,query, TALK=TALK)
            }
        }
    }


    dbDisconnect(con)


}
