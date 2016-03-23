#' Configuration
#+ echo=FALSE
library(printr)
data_dir <- '/data'
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
  stdout <- system(paste(plink_path, '--noweb', '--bfile', infile, '--make-bed', '--out', outfile, ...), intern = T)
  print(grep('done\\.$', stdout, invert = T, value = T))
  return(outfile)
}

to_file <- function(x) {
  file_name <- tempfile(tmpdir = temp_dir)
  write.table(x, file=file_name, col.names = F, row.names = F, quote = F)
  return(file_name)
}

king_path <- '/usr/bin/king'
king <- function (infile) {
  outfile = tempfile(tmpdir = temp_dir)
  print(system(paste(king_path, '-b', paste0(infile, '.bed'), '--kinship', '--prefix', outfile)))
  return(list('all'=read.table(paste0(outfile, '.kin0'), sep='\t', header = T), 'family'=read.table(paste0(outfile, '.kin'), sep='\t', header = T)))
}

#' # Update GenomeStudio result with data from phentoype.csv
phenotype <- read.csv(paste(data_dir, 'phenotype.csv', sep='/'), sep='\t')

#' Convert exported text file to binary file
if (!file.exists(paste(raw_genotypes, 'bed', sep='.'))) {
  system(paste(plink_path, '--noweb', '--file', paste(genotype_dir, 'easd_nh_2015-11-13-101820', sep='/'), '--make-bed', '--out', raw_genotypes))
}

#' Update FIDs
genotypes <- plink('--update-ids',
                   to_file(subset(phenotype, select=c(OMICRON_ID, OMICRON_ID, FAMILY_ID, OMICRON_ID))),
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
  # Exclude markers on Y chromosome
  missing_per_marker <- subset(missing_per_marker, CHR != 24)
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
table(sex_check$PLATE, sex_check$STATUS)
qplot(seq_along(F), F, color=factor(STATUS), data=sex_check) + geom_hline(yintercept=0.2)

qplot(ARRAY_COL, ARRAY_ROW, fill=factor(STATUS), facets = ~ARRAY_ID, data=sex_check, geom='tile') + ggtitle('Sex check per sample vs location by array')
qplot(WELL_COL, WELL_ROW, fill=factor(STATUS), facets = ~PLATE, data=sex_check, geom='tile') + ggtitle('Sex check per sample vs location by plate')

# TODO load this data from beeline final report!
# #' Plot probes intensities on X and Y chromosomes
# xy_snps <- read.csv(pipe("tail -n +9 /resources/arrays/HumanCore-12-v1/HumanCore-12-v1-0-B.csv | awk -F , '$10 == \"X\" || $10 == \"Y\"' | cut -d , -f 2,10,11"), header=F)
# names(xy_snps) <- c('name', 'chr', 'position')
#
# TODO fixme
# xy_b_allele_freq <- read.csv(pipe(paste("zcat", paste(genotype_dir, 'easd_nh_final_report.txt_Final.csv.gz', sep='/'), "| tail -n +9 | cut -d , -f 1,2,4,10,11")), header=T)
# mean_xy_intensity <- ddply(xy_b_allele, .(Sample.ID, Chr), summarize, mean=mean(Y, na.rm = T))
# mean_xy_intensity <- merge(mean_xy_intensity, sex_check, by.x='Sample.ID', by.y='IID')
# library(reshape2)
# qplot(X, Y, color=STATUS, data=dcast(mean_xy_intensity, Sample.ID + STATUS + SNPSEX + COUNTRY ~ Chr, value.var = 'mean'), geom='point')

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
sex_check_reproducibility_failed <- replicates_sex_check[replicates_sex_check$SAMPLE_ID %in% subset(sex_check_reproducibility, V1==F)$SAMPLE_ID,1:6]
sex_check_reproducibility_failed

# TODO fixme
# #' Plot B allele frequency for samples flagged by replicate sex check
# xy_b_allele_freq <- subset(xy_b_allele, Sample.ID %in% sex_check_reproducibility_failed$IID)
# xy_b_allele_freq <- merge(xy_b_allele_freq, replicates_sex_check[,c('IID', 'PEDSEX', 'SAMPLE_ID')], by.x='Sample.ID', by.y='IID', all.x=T)
# #+ fig.height=10
# qplot(Position, B.Allele.Freq, data=xy_b_allele_freq, geom='point') + facet_grid(Sample.ID + SAMPLE_ID  ~ Chr, scales = 'free', space='free')

#' Drop all samples flagged by sex check
genotypes <- plink('--remove', to_file(subset(sex_check, STATUS=='PROBLEM', select=c(FAMILY_ID, IID))), infile=prefiltered_genotypes)

#' # Filtering
#' MAF fitering
genotypes <- plink('--maf', '0.01', infile=genotypes)

#' Filter by missingnes per marker
genotypes <- plink('--geno', '0.1', infile=genotypes)

#' Filter by missingnes per subject
genotypes <- plink('--mind', '0.01', infile=genotypes)

#' Remove SNPs with het. haploid genotypes
het_hap_to_remove <- read.table(paste(genotypes, 'hh', sep='.'), header = F)
names(het_hap_to_remove) <- c('FID', 'IID', 'marker')
filtered_genotypes <- plink('--exclude', to_file(as.character(unique(het_hap_to_remove$marker))), infile=genotypes)

#' # Cleanup
#' Analyze heterozygosity
genotypes <- plink('--het', infile=filtered_genotypes)
heterozygosity <- read.table(paste(genotypes, 'het', sep='.'), header = T)
het_mean <- mean(heterozygosity$F)
het_sd <- sd(heterozygosity$F)
i <- 2
length(which(heterozygosity$F < het_mean-i*het_sd | heterozygosity$F > het_mean+i*het_sd))
qplot(F, data=heterozygosity) + geom_vline(xintercept=c(het_mean+i*het_sd, het_mean-i*het_sd), color='red')

#' Plot relatedness - plink ibd/King
#' Kinship coefficient! See: http://cphg.virginia.edu/quinlan/?p=300
kinship <- king(genotypes)
qplot(Kinship, data=kinship[['all']])
qplot(seq_along(Kinship), Kinship, data=kinship[['all']]) +
  geom_hline(yintercept=c(.5, .375, .25, .125, .0938, .0625, .0313, .0078, .002), color='blue') +
  scale_y_log10()
# Cryptic duplicates
subset(kinship[['all']], Kinship == 0.5)

qplot(Kinship, data=kinship[['family']])
qplot(seq_along(Kinship), Kinship, data=kinship[['family']]) +
  geom_hline(yintercept=c(.5, .375, .25, .125, .0938, .0625, .0313, .0078, .002), color='blue') +
  scale_y_log10()
# Cryptic duplicates
subset(kinship[['family']], Kinship == 0.5)

#' Population structure

#' HWE filtering
# TODO record and plot - do not remove
# genotypes <- plink('--hwe', '1e-5', infile=genotypes)
