# create a FRANz input file for a single population from 
#   - a standard set of PLINK files  
#   - a file with id - sex - birthyear (assuming 1=female, 2=male, 3=unknown)
# and run franz from within R on a windows machine with cygwin installed,
# and FRANz installed according to the user manual
# requires package plyr


#==========================================

MakeFranzIN <- function(PlinkFile,  # without extension
                        bfile = FALSE, 
                        LHfile, 
                        OutFile)  # without extension
{
  if (bfile) {
    system("cmd", input = paste("plink --bfile", PlinkFile, "--recodeA --out",  PlinkFile))
  } else {
    system("cmd", input = paste("plink --file", PlinkFile, "--recodeA --out",  PlinkFile))
  }
  system("cmd", input = paste0('cat ', PlinkFile, ".raw",
                              ' | sed -e "s/NA/-9/g"> ', PlinkFile, '.txt'))  # TO BE TESTED !!!
  GenoTmp <- readLines(paste0(PlinkFile, ".txt"))
  TmpL    <- strsplit(GenoTmp[-1], split=" ")  # skip header row
  Geno <- plyr::ldply(TmpL)
  IDs <- Geno[,2]
  Geno <- as.matrix(Geno[, -(1:6)])
  GenoF <- Geno  # same dimensions
  dc <- c("0" = "1/1", "1"="1/2", "2" = "2/2", "-9" = "?/?")
  for (i in 1:nrow(Geno)) {
    GenoF[i,] <- sapply(Geno[i,], function(z) dc[z])
  }
  
  LH <- read.table(LHfile, header=TRUE, stringsAsFactors=FALSE)
  tmp <- unique(merge(as.data.frame(IDs), LH, by.x="IDs", by.y=names(LH)[1],
                      all.x=TRUE, sort=FALSE))
  LHF <- data.frame(ID = tmp$IDs,
                    nClones = 1,
                    BY = tmp$BY,
                    DY = "?",  # death year - assume unknown
                    Sex = c("F", "M", "?")[tmp$Sex])  
  GenoF <- cbind(LHF, GenoF)
  GenoF[,1] <- sprintf("%-10s", IDs)
  
  cat(1, ncol(Geno), "/", OutFile, "\n", file=paste0(OutFile, ".dat"), append=FALSE)
  cat(nrow(Geno), "Pop \n", file=paste0(OutFile, ".dat"), append=TRUE)
  write.table(GenoF, paste0(OutFile, ".dat"), quote=FALSE, col.names=FALSE,
              row.names=FALSE, sep=" ", append=TRUE, na="?")
}


#==========================================
# Run FRANz:

repro <- "1:20"   # for example
CygwinPath <- "c:\\cygwin64\\bin\\env"  # escaped windows-style backward slash
FRANzPath <- "/cygdrive/E/Manuscripts/Sequoia/FRANz/FRANz-2.0.0/src/FRANz.exe"  # unix-style forward slash
system(paste('cmd.exe /c', CygwinPath, FRANzPath,
             '--Nmax 500 --femrepro', repro, '--malerepro', repro, '--pedigreeoutformat 3',
             paste0(OutFile,".dat")), 
       TRUE)


# optionally: replace * by NA for non-assigned parents
system("cmd", input = paste('cat', "pedigree.txt",
                            '| sed -e "s/*/NA/g">',
                            "PedigreeFranz.txt"))
