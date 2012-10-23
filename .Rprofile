upVers <- function(path,update="snapshot",date=TRUE,simplify=TRUE)
{
  # This function updates the description file from package
  # in path (assumed work directory by default, as typical
  # with projects in RStudio using GitHub).

  # Usage:
    # path: path to contents of a package
    # update: What to update? "version", "major", "minor", "snapshot"
    # date: Update date as well?
    # simplfy: omit trailing zeros?

  # Assumes following numbering system:
  # version.major.minor-snapshot

  uplist <- c("version","major","minor","snapshot")

  if (missing(path)) path <- getwd()
  DESCfile <- paste0(path,"/DESCRIPTION")
  if (!file.exists(DESCfile)) stop("DESCRIPTION does not exist. Is this the folder of a package?")

  DESC <- readLines(DESCfile)

  ### Update date:
  if (date)
  {
    DESC <- gsub("(?<=Date: )\\d{4}-\\d{2}-\\d{2}",Sys.Date(),DESC,perl=TRUE)
  }

  ### Update version:
  Vers <- regmatches(DESC,regexpr("(?<=Version: )\\d+\\.?\\d*\\.?\\d*\\-?\\d*",DESC,perl=TRUE))
  Vers <- as.numeric(unlist(strsplit(Vers,split="\\.|\\-")))
  Vers <- c(Vers,rep(0,length=4-length(Vers)))
  Vers[grep(update,uplist,ignore.case=TRUE)] <- Vers[grep(update,uplist,ignore.case=TRUE)] + 1
  Vers[1:4>grep(update,uplist,ignore.case=TRUE)] <- 0

  # Combine and replace:
  Vers <- paste(paste(Vers[1:3],collapse="."),Vers[4],sep="-")
  if (simplify)
  {
    Vers <- gsub("\\.?0?\\.?0?\\-?0?$","",Vers)
  }
  DESC <- gsub("(?<=Version: )\\d+\\.?\\d*\\.?\\d*\\-?\\d*",Vers,DESC,perl=TRUE)

  # Write Description:
  writeLines(DESC,DESCfile)
}

up <- utils:::menu(c("Yes","No"),TRUE,"Update snapshot?")
if (up==1)
{
	upVers()
	cat("Snapshot updated!")
} else cat("Snapshot NOT updated!")