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
  system(paste(plink_path, '--noweb', '--file', paste(genotype_dir, 'PLINK_250715_1022', 'easd_nh', sep='/'), '--make-bed', '--out', raw_genotypes))
}

#' Update FIDs
genotypes <- plink('--update-ids',
                   to_file(subset(phenotype, select=c(GENOMESTUDIO_FID, OMICRON_ID, FAMILY_ID, OMICRON_ID))),
                   infile=raw_genotypes)

#' Update Sex field
genotypes <- plink('--update-sex',
                   to_file(subset(phenotype, select=c(FAMILY_ID, OMICRON_ID, SEX_PLINK))),
                   infile=genotypes)

#' Drop excluded individuals (not included in phentoype file)
updated_genotypes <- plink('--keep',
                   to_file(subset(phenotype, select=c(FAMILY_ID, OMICRON_ID))),
                   infile=genotypes)

#' Plot missinges
library(ggplot2)

missingnes <- function (genotypes) {
  missing_per_sample <- read.table(paste(genotypes, 'imiss', sep='.'), header = T)
  missing_per_sample <- merge(missing_per_sample, phenotype, by.x='IID', by.y='OMICRON_ID', all.x=T)
  print(qplot(factor(ARRAY_ID), F_MISS, fill=factor(PLATE), data=missing_per_sample, geom = 'boxplot') + coord_flip()+ ggtitle('Fraction missing per sample by array'))
  print(qplot(ARRAY_COL, ARRAY_ROW, fill=F_MISS, facets = ~ARRAY_ID, data=missing_per_sample, geom='tile') + scale_fill_gradient(low="white", high="red") + ggtitle('Fraction missing per sample vs location by array'))
  print(qplot(WELL_COL, WELL_ROW, fill=F_MISS, facets = ~PLATE, data=missing_per_sample, geom='tile') + scale_fill_gradient(low="white", high="red") + ggtitle('Fraction missing per sample vs location by plate'))
  missing_per_marker <- read.table(paste(genotypes, 'lmiss', sep='.'), header = T)
  print(qplot(factor(CHR), F_MISS, color=factor(CHR), data = missing_per_marker, geom='boxplot') + ggtitle('Missingnes per marker'))
  return(missing_per_sample)
}

missing_per_sample <- missingnes(plink('--missing', infile=updated_genotypes))

#' Remove replicates by missingnes
library(plyr)
replicates_to_remove <- ddply(subset(missing_per_sample, grepl('.*_.*', missing_per_sample$IID)), .(SAMPLE_ID), function (df){ df[order(df$F_MISS)[2:nrow(df)], c('FAMILY_ID', 'IID')] })[,2:3]
genotypes <- plink('--remove', to_file(replicates_to_remove), infile=updated_genotypes)

#' Drop failed arrays
mean_missingnes <- ddply(missing_per_sample, .(ARRAY_ID), summarize, mean=mean(F_MISS))
prefiltered_genotypes <- plink('--remove', to_file(subset(phenotype, ARRAY_ID %in% mean_missingnes[mean_missingnes$mean > 0.5, 1], select=c('FAMILY_ID', 'OMICRON_ID'))), infile=genotypes)

#' Plot missinges after removing replicates
invisible(missingnes(plink('--missing', infile=prefiltered_genotypes)))

#' Sex check
genotypes <- plink('--merge-x', infile=prefiltered_genotypes)
genotypes <- plink('--split-x', 'hg19', infile=genotypes)
genotypes <- plink('--check-sex', infile=genotypes)
sex_check <- read.table(paste(genotypes, 'sexcheck', sep='.'), header = T)
sex_check <- merge(sex_check, phenotype, by.x='IID', by.y='OMICRON_ID', all.x=T)
table(sex_check$STATUS)
table(sex_check$COUNTRY, sex_check$STATUS)
qplot(seq_along(F), F, color=factor(STATUS), data=sex_check) + geom_hline(y=0.2)

qplot(ARRAY_COL, ARRAY_ROW, fill=factor(STATUS), facets = ~ARRAY_ID, data=sex_check, geom='tile') + ggtitle('Sex check per sample vs location by array')
qplot(WELL_COL, WELL_ROW, fill=factor(STATUS), facets = ~PLATE, data=sex_check, geom='tile') + ggtitle('Sex check per sample vs location by plate')

#' Verify probes intensity on X and Y chromosomes
xy_b_allele_intensity <- read.table(paste(genotype_dir, 'easd_nh_xy_cr95_gc6.txt', sep='/'), header = T, skip=9, sep='\t')
xy_b_allele_intensity <- subset(xy_b_allele_intensity, GC.Score >= 0.15)
mean_xy_intensity <- ddply(xy_b_allele_intensity, .(Sample.ID, Chr), summarize, mean=mean(Y, na.rm = T))
mean_xy_intensity <- merge(mean_xy_intensity, sex_check, by.x='Sample.ID', by.y='IID')
tmp <- dcast(mean_xy_intensity, Sample.ID + STATUS + SNPSEX + COUNTRY ~ Chr, value.var = 'mean')
qplot(X, Y, color=STATUS, data=tmp, geom='point')

#' Sex check on replicates
replicates <- subset(phenotype, grepl('.*_.*', phenotype$OMICRON_ID), select=c('FAMILY_ID', 'OMICRON_ID'))
genotypes <- plink('--keep', to_file(replicates), infile=updated_genotypes)
genotypes <- plink('--merge-x', infile=genotypes)
genotypes <- plink('--split-x', 'hg19', infile=genotypes)
genotypes <- plink('--check-sex', infile=genotypes)
replicates_sex_check <- read.table(paste(genotypes, 'sexcheck', sep='.'), header = T)
replicates_sex_check <- merge(replicates_sex_check, phenotype, by.x='IID', by.y='OMICRON_ID', all.x=T)
replicates_sex_check$SAMPLE_ID <- factor(replicates_sex_check$SAMPLE_ID)
table(replicates_sex_check$STATUS)

sex_check_reproducibility <- ddply(replicates_sex_check, .(SAMPLE_ID), function (df){ all(df$STATUS == df$STATUS[1]) })
sum(sex_check_reproducibility$V1)/length(unique(replicates_sex_check$SAMPLE_ID))
sum(!sex_check_reproducibility$V1)
replicates_sex_check[replicates_sex_check$SAMPLE_ID %in% subset(sex_check_reproducibility, V1==F)$SAMPLE_ID,1:6]

#' # Filtering
#' MAF fitering
genotypes <- plink('--maf', '0.01', infile=prefiltered_genotypes)

#' Filter by missingnes per marker
genotypes <- plink('--geno', '0.1', infile=genotypes)

#' HWE filtering
# TODO record - do not remove
genotypes <- plink('--hwe', '1e-5', infile=genotypes)

#' Filter by missingnes per subject
genotypes <- plink('--mind', '0.02', infile=genotypes)

# TODO Remove tmp dir when done?!?