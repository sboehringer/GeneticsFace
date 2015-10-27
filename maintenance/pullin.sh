# pullin.sh
# refresh externally maintained files

export CWD=`pwd`
( cd ~/src/Rprivate ; cp Rdata.R Rgraphics.R Rsystem.R $CWD/../GeneticsFace/R )
