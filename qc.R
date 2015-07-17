setwd("~/projekty/easd_nh")

plink_path <- '/usr/bin/p-link'

temp_dir = tempdir()
data_dir = 'data'

plink <- function(..., infile, outfile = NULL){
  if (is.null(outfile)) {
    outfile = tempfile(tmpdir = tmp)
  }
  system(paste(plink_path, '--noweb', '--bfile', infile, '--make-bed', '--out', outfile, ...))
  return(outfile)
}

#' # Update GenomeStudio result with data from phentoype.csv
phenotype <- read.csv('data/phenotype.csv', sep='\t')

#' Convert file exported from Genome studio to binary file
gs_export <- tempfile(tmpdir = temp_dir)
system(paste(plink_path, '--noweb', '--file', paste(c(data_dir, 'easd_nh', sep='/')), '--make-bed', '--out', gs_export))

#' Update FIDs
update_ids_file <- tempfile(tmpdir = temp_dir)
write.table(subset(phenotype, select=c(GenomeStudio_FID, OMICRON_ID, FAMILY_ID, OMICRON_ID)), file=update_ids_file, col.names = F)
updated_ids <- plink('--update-ids', update_ids_file, infile=gs_export)

#' Update Affected and Sex fields

#' # Filtering
#' Drop individuals excluded individuals, eg. NEWBORN-*

#' MAF filter
#' CR filter
#' HWE filter

#' Remove tmp dir when done?!?