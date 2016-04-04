#' Configuration
library(RcppEigen)
library(printr)
library(plyr)
library(ggplot2)
library(plink2R) # install_github("gabraham/plink2R", subdir = 'plink2R')

data_dir <- '/data'
phenotype_dir <- paste(data_dir, 'phenotype', sep='/')
genotype_dir <- paste(data_dir, 'genotype', sep='/')
raw_genotypes <- paste(data_dir, 'genotype', 'easd_nh', sep='/')
resources_dir <- '/resources/arrays'
temp_dir = tempdir()

#' PLINK wrapper
plink_path <- '/usr/bin/plink'
plink <- function(..., infile, outfile = NULL){
  if (is.null(outfile)) {
    outfile = tempfile(tmpdir = temp_dir)
  }
  stdout <- system(paste(plink_path, '--noweb', '--bfile', infile, '--make-bed', '--out', outfile, ...), intern = T)
  print(grep('done\\.$', stdout, invert = T, value = T))
  return(outfile)
}

#' data.frame to TSV file
to_file <- function(x) {
  file_name <- tempfile(tmpdir = temp_dir)
  write.table(x, file=file_name, col.names = F, row.names = F, quote = F)
  return(file_name)
}

#' Compare two steps of the analysis and print sample statistics
sample_stats <- function(genotype_a, genotype_b) {
  stats <- list()
  a <- read_plink(genotype_a)
  b <- read_plink(genotype_b)
  stats$snps <- c(dim(a$bim)[1], dim(b$bim)[1])
  stats$individuals <- c(dim(a$fam)[1], dim(b$fam)[1])
  if (stats$individuals[1] != stats$individuals[2]) {
    stats$individuals_removed <- setdiff(apply(a$fam[,1:2], 1, function (row) paste(row[1], row[2])), apply(b$fam[,1:2], 1, function (row) paste(row[1], row[2])))
  }
  if (stats$snps[1] != stats$snps[2]) {
    stats$snps_removed <- stats$snps[1] - stats$snps[2]
  }
  stats
}

#' # Update raw genotypes with data from phentoype.csv
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
genotypes_with_phenotype <- plink('--keep',
                   to_file(subset(phenotype, select=c(FAMILY_ID, OMICRON_ID))),
                   infile=genotypes)
sample_stats(genotypes, genotypes_with_phenotype)

#' Remove SNPs with unknown locaton
genotypes_no_chr0 <- plink('--not-chr', 0, infile = genotypes_with_phenotype)
sample_stats(genotypes_with_phenotype, genotypes_no_chr0)

#' Plot missinges
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

missing_per_sample <- missingnes(plink('--missing', infile=genotypes_no_chr0))

#' Drop failed plate #16
genotypes_no_plate_16 <- plink('--remove', to_file(subset(phenotype, PLATE == 16, select=c('FAMILY_ID', 'OMICRON_ID'))), infile=genotypes_no_chr0)
sample_stats(genotypes_no_chr0, genotypes_no_plate_16)

#' Drop failed arrays
# TODO remove - does not remove any arrays
mean_missingnes <- ddply(missing_per_sample, .(ARRAY_ID), summarize, mean=mean(F_MISS))
genotypes_no_failures <- plink('--remove', to_file(subset(phenotype, ARRAY_ID %in% mean_missingnes[mean_missingnes$mean > 0.5, 1], select=c('FAMILY_ID', 'OMICRON_ID'))), infile=genotypes_no_plate_16)
sample_stats(genotypes_no_plate_16, genotypes_no_failures)

#' Plot missinges after removing replicates
invisible(missingnes(plink('--missing', infile=genotypes_no_failures)))

#' Filter by missingnes per subject
genotypes_mind_1pct <- plink('--mind', '0.01', infile=genotypes_no_failures)
sample_stats(genotypes_no_failures, genotypes_mind_1pct)

#' Sex check
genotypes <- plink('--merge-x', infile=genotypes_mind_1pct)
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
genotypes <- plink('--keep', to_file(replicates), infile=genotypes_mind_1pct)
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

#' Drop all samples flagged by sex check
genotypes_no_sex_errors <- plink('--remove', to_file(subset(sex_check, STATUS=='PROBLEM', select=c(FAMILY_ID, IID))), infile=genotypes_mind_1pct)
sample_stats(genotypes_mind_1pct, genotypes_no_sex_errors)

#' Drop chr >= 23
genotypes_autosomes <- plink('--chr', '1-22', infile = genotypes_no_sex_errors)
sample_stats(genotypes_no_sex_errors, genotypes_autosomes)

#' Analyze heterozygosity
genotypes <- plink('--het', infile=genotypes_autosomes)
heterozygosity <- read.table(paste(genotypes, 'het', sep='.'), header = T)
heterozygosity$H <- heterozygosity$O.HOM/heterozygosity$N.NM*100
qplot(H, data=heterozygosity) +
  geom_vline(xintercept=mean(heterozygosity$H)+2*sd(heterozygosity$H), color='red') +
  geom_vline(xintercept=mean(heterozygosity$H)-2*sd(heterozygosity$H), color='red')
# TODO remove samples with extreme heterozygosity

#' Remove replicates by missingnes
genotypes <- plink('--missing', infile=genotypes_autosomes)
missing_per_sample <- read.table(paste(genotypes, 'imiss', sep='.'), header = T)
missing_per_sample <- merge(missing_per_sample, phenotype, by.x='IID', by.y='OMICRON_ID', all.x=T)
replicates_to_remove <- ddply(subset(missing_per_sample, grepl('.*_.*', missing_per_sample$IID)), .(SAMPLE_ID), function (df){ df[order(df$F_MISS)[2:nrow(df)], c('FAMILY_ID', 'IID')] })[,2:3]
replicates_to_remove <- replicates_to_remove[complete.cases(replicates_to_remove),]
genotypes_unique <- plink('--remove', to_file(replicates_to_remove), infile=genotypes_autosomes)
sample_stats(genotypes_autosomes, genotypes_unique)

#' KING wrapper
king_path <- '/usr/bin/king'
king <- function (infile) {
  outfile = tempfile(tmpdir = temp_dir)
  print(system(paste(king_path, '-b', paste0(infile, '.bed'), '--kinship', '--prefix', outfile)))
  return(list('all'=read.table(paste0(outfile, '.kin0'), sep='\t', header = T), 'family'=read.table(paste0(outfile, '.kin'), sep='\t', header = T)))
}

#' Kinship coefficient! See: http://cphg.virginia.edu/quinlan/?p=300
kinship <- king(genotypes_unique)
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

# TODO remove cryptic duplicates

#' Filter by missingnes per marker
genotypes_geno_05 <- plink('--geno', '0.05', infile=genotypes_unique)
sample_stats(genotypes_unique, genotypes_geno_05)

# TODO filter by geno & maf!
#' HWE filtering
# genotypes <- plink('--hwe', '1e-5', infile=genotypes)

#' MAF fitering
genotypes_maf_01 <- plink('--maf', '0.01', infile=genotypes_geno_05)
sample_stats(genotypes_geno_05, genotypes_maf_01)

#' # Population structure
convertf_path <- '/usr/bin/convertf'
smartpca_path <- '/usr/bin/smartpca'

#' Extract HapMap SNPs
hapmap_snps <- plink('--extract', paste(resources_dir, 'hapmap3', 'hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt', sep='/'), infile=genotypes_maf_01)

#' LD pruning
genotypes_pruned <- plink('--exclude', paste(resources_dir, 'high_ld_hg19.txt', sep='/'), '--range', '--indep-pairwise 50 5 0.2', infile=hapmap_snps)

genotypes_with_hapmap <- plink('--bmerge', paste(resources_dir, 'hapmap3',
                    'hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps', sep='/'),
                   '--extract', paste(genotypes_pruned, 'prune.in', sep='.'), infile=genotypes_pruned)

#' Flip failed SNPs
genotypes_flipped <- plink('--extract', paste(resources_dir, 'hapmap3', 'hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt', sep='/'),
                   '--flip', paste0(genotypes_with_hapmap, '-merge.missnp'), infile=hapmap_snps)

genotypes_with_hapmap_2 <- plink('--bmerge', paste(resources_dir, 'hapmap3',
                                                 'hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps', sep='/'),
                               '--extract', paste(genotypes_pruned, 'prune.in', sep='.'), infile=genotypes_flipped)

genotypes <- plink('--exclude', paste0(genotypes_with_hapmap_2, '-merge.missnp'), '--extract', paste(genotypes_pruned, 'prune.in', sep='.'), infile = genotypes_flipped)

genotypes_with_hapmap_3 <- plink('--bmerge', paste(resources_dir, 'hapmap3',
                                                   'hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps', sep='/'),
                                 '--extract', paste(genotypes_pruned, 'prune.in', sep='.'), infile=genotypes)

# Convert to PED/MAP
genotypes_merged <- tempfile(tmpdir = temp_dir)
system(paste(plink_path, '--noweb', '--bfile', genotypes_with_hapmap_3, '--recode', '--output-missing-phenotype', '1', '--out', genotypes_merged))

# Setup convertf
convertf_params <- tempfile(tmpdir = temp_dir)
eigenstrat_input <- tempfile(tmpdir = temp_dir)
cat(paste0('genotypename:    ', genotypes_merged, '.ped\n'), file=convertf_params)
cat(paste0('snpname:         ', genotypes_merged, '.map\n'), file=convertf_params, append = T)
cat(paste0('indivname:       ', genotypes_merged, '.ped\n'), file=convertf_params, append = T)
cat('outputformat:    EIGENSTRAT\n', file=convertf_params, append = T)
cat(paste0('genotypeoutname: ', eigenstrat_input, '.eigenstratgeno\n'), file=convertf_params, append = T)
cat(paste0('snpoutname:      ', eigenstrat_input, '.snp\n'), file=convertf_params, append = T)
cat(paste0('indivoutname:    ', eigenstrat_input, '.ind\n'), file=convertf_params, append = T)
cat('familynames:     NO\n', file=convertf_params, append = T)

# Convert PED/MAP to EIGENSTRAT
system(paste(convertf_path, '-p', convertf_params))

# Setup smartcpa
smartcpa_params <- tempfile(tmpdir = temp_dir)
eigenstrat_output <- tempfile(tmpdir = temp_dir)

cat(paste0('genotypename: ', eigenstrat_input, '.eigenstratgeno\n'), file=smartcpa_params)
cat(paste0('snpname:      ', eigenstrat_input, '.snp\n'), file=smartcpa_params, append = T)
cat(paste0('indivname:    ', eigenstrat_input, '.ind\n'), file=smartcpa_params, append = T)
cat(paste0('evecoutname:  ', eigenstrat_output, '.evec\n'), file=smartcpa_params, append = T)
cat(paste0('evaloutname:  ', eigenstrat_output, '.eval\n'), file=smartcpa_params, append = T)
cat('altnormstyle: NO\n', file=smartcpa_params, append = T)
cat('numoutevec:     10\n', file=smartcpa_params, append = T)
cat('numoutlieriter: 0\n', file=smartcpa_params, append = T)
cat('numoutlierevec: 2\n', file=smartcpa_params, append = T)
cat('outliersigmathresh: 8.0\n', file=smartcpa_params, append = T)
cat('qtmode: 0\n', file=smartcpa_params, append = T)
cat('nsnpldregress: 2\n', file=smartcpa_params, append = T)
# cat(paste0('evaloutname:  ', eigenstrat_output, '.outlier\n'), file=smartcpa_params, append = T)

# Run smartpca
system(paste(smartpca_path, '-p', smartcpa_params))

population_structure <- read.table(paste0(eigenstrat_output, '.evec'), header = F, skip=1)
names(population_structure) <- c('IID', paste0('PC', 1:10), 'group')

#' Add group labels from HapMap3
relationships_w_pops <- read.table(paste(resources_dir, 'hapmap3', 'relationships_w_pops_121708.txt', sep='/'), header = T)
relationships_w_pops <- relationships_w_pops[,c(2,7)]
population_structure <- merge(population_structure, relationships_w_pops, by='IID', all.x=T)
population_structure$population <- as.character(population_structure$population)
population_structure[is.na(population_structure$population),'population'] <- 'EASD'
population_structure$population <- factor(population_structure$population)
population_structure <- subset(population_structure, select=-group)

#' Population structure PC1 vs PC2
qplot(PC1, PC2, color=population, data=population_structure)
