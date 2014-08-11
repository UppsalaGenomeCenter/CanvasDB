############################################################################
##
## File: canvasDB_config.R
##
## Author: Adam Ameur, Uppsala Genome Center
##
## Description: Global variables and configurations needed for the canvasDB
##
############################################################################

tmpfileDir <- paste(rootDir,"tmpfile_dir/",sep="") # Directory for temporary files.

if(!file.exists(tmpfileDir)){
    dir.create(tmpfileDir)
}

indbDatabaseName <- "canvasdb_test"
annotationDatabaseName <- "canvasdb_annot_hg19"

ANNOVARpath <- "/Volumes/Data/tools/annovar/annotate_variation.pl" ## Add path to ANNOVAR executable
ANNOVARpathDB <- "/Volumes/Data/tools/annovar/humandb/" ## Add path to ANNOVAR annotation database

dbSNPversion <- 137

mysqlConfigFile <- ".my.cnf"  # File containing MySQL user info

mysqlEngine <- "MyISAM"

## MySQL database tables
annotTables <- list()
annotTables[["score"]] <- paste(annotationDatabaseName,".annot_score",sep="")
annotTables[["dbsnp"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_single",sep="")
annotTables[["dbsnpCommon"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_single_common",sep="")
annotTables[["dbsnpIndels"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_indels",sep="")
annotTables[["dbsnpCommonIndels"]] <- paste(annotationDatabaseName,".annot_snp",dbSNPversion,"_indels_common",sep="")


dbTables <- list()
dbTables[["run"]] <- "runs"
dbTables[["sample"]] <- "samples"
dbTables[["SNP.summary"]] <- "snp_summary"
dbTables[["indel.summary"]] <- "indel_summary"

## This should stop all forms of "Scientific writing, ie 1e+05"
options(scipen = 999)
