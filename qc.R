#' Configuration
#+ echo=FALSE
data_dir <- '/data/easd_nh'
phenotype_dir <- paste(data_dir, 'phenotype', sep='/')
genotype_dir <- paste(data_dir, 'genotype', sep='/')
raw_genotypes <- paste(data_dir, 'genotype', 'easd_nh', sep='/')
temp_dir = tempdir()

#' PLINK wrapper
#+ echo=FALSE
plink_path <- '/usr/bin/plink'
plink <- function(..., infile, outfile = NULL){
  if (is.null(outfile)) {
    outfile = tempfile(tmpdir = temp_dir)
  }
  print(system(paste(plink_path, '--noweb', '--bfile', infile, '--make-bed', '--out', outfile, ...), intern = T))
  return(outfile)
}

to_file <- function(x) {
  file_name <- tempfile(tmpdir = temp_dir)
  write.table(x, file=file_name, col.names = F, row.names = F, quote = F)
  return(file_name)
}

#' # Update GenomeStudio result with data from phentoype.csv
phenotype <- read.csv(paste(data_dir, 'phenotype.csv', sep='/'), sep='\t')

#' Convert file exported from Genome studio to binary file
if (!file.exists(paste(raw_genotypes, 'bed', sep='.'))) {
  system(paste(plink_path, '--noweb', '--file', paste(genotype_dir, 'PLINK_090715_0107', 'easd_nh', sep='/'), '--make-bed', '--out', raw_genotypes))
}

#' Update FIDs
genotypes <- plink('--update-ids',
                   to_file(subset(phenotype, select=c(GENOMESTUDIO_FID, OMICRON_ID, FAMILY_ID, OMICRON_ID))),
                   infile=raw_genotypes)

#' Update Affected and Sex fields
genotypes <- plink('--update-sex',
                   to_file(subset(phenotype, select=c(FAMILY_ID, OMICRON_ID, SEX_PLINK))),
                   infile=genotypes)

#' # Filtering
#' Drop excluded individuals (not included in phentoype file)
genotypes <- plink('--keep',
                   to_file(subset(phenotype, select=c(FAMILY_ID, OMICRON_ID))),
                   infile=genotypes)

#' MAF fitering
genotypes <- plink('--maf', '0.01', infile=genotypes)

#' Missingnes per marker
genotypes <- plink('--geno', '0.1', infile=genotypes)

#' HWE filtering
# TODO record - do not remove
genotypes <- plink('--hwe', '1e-5', infile=genotypes)

#' Missingnes per subject
genotypes <- plink('--mind', '0.02', infile=genotypes)

#' Sex check
genotypes <- plink('--check-sex', infile=genotypes)

# TODO Remove tmp dir when done?!?