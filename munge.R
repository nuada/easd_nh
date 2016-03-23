library(xlsx)
library(lubridate)

data_dir <- '/data'
phenotype_dir <- paste(data_dir, 'phenotype', sep='/')

# CZ - Czech Republic
columns_cz <- list(
  c('FAMILY_ID',                   'character'),
  c('SAMPLE_ID',                   'character'),
  c('SEX',                         'character'),
  c('DOB',                         'Date'),
  c('AGE_AT_DIAGNOSIS',            'numeric'),
  c('AGE_AT_CLINICAL_EXAM',        'numeric'),
  c('HEIGHT',                      'numeric'),
  c('WEIGHT',                      'numeric'),
  c('BMI_AT_CLINICAL_EXAM',        'numeric'),
  c('MUTATION_AA_CHANGE',          'character'),
  c('EXON',                        'character'),
  c('MUTATION_CODING',             'character'),
  c('MUTATION_INHERITANCE_FATHER', 'character'),
  c('MUTATION_INHERITANCE_MOTHER', 'character'),
  c('INSULIN_START_AGE',           'numeric'),
  c('THERAPY_AT_DIAGNOSIS',        'character'),
  c('THERAPY_AT_REFERRAL',         'character')
)
phenotype_cz <- read.xlsx(paste(phenotype_dir, 'ClinicalDatabase-Czech2015.xls', sep='/'),
                          sheetName = 'CR',
                          stringsAsFactors = F,
                          colIndex = 2:18,
                          colClasses = sapply(columns_cz, function(c)c[2]))
names(phenotype_cz) <- sapply(columns_cz, function(c)c[1])
phenotype_cz$COUNTRY <- 'CZ'

phenotype_cz$EXON_NO <- as.numeric(sub('^.* ([0-9]+)', '\\1', phenotype_cz$EXON, perl=T))
phenotype_cz$EXON <- grepl('^exon', trimws(phenotype_cz$EXON), ignore.case = T)
phenotype_cz$EXON[is.na(phenotype_cz$EXON_NO)] <- NA

phenotype_cz$MUTATION_INHERITANCE <- NA
phenotype_cz$MUTATION_INHERITANCE[phenotype_cz$MUTATION_INHERITANCE_FATHER == '1'] <- 'F'
phenotype_cz$MUTATION_INHERITANCE[phenotype_cz$MUTATION_INHERITANCE_MOTHER == '1'] <- 'M'
phenotype_cz <- subset(phenotype_cz, select=-c(MUTATION_INHERITANCE_FATHER, MUTATION_INHERITANCE_MOTHER))

# FR - France
columns_fr <- list(
  c('FAMILY_ID',            'character'),
  c('SAMPLE_ID',            'character'),
  c('SEX',                  'character'),
  c('PROBAND',              'logical'), # NA <- False
  c('RELATIVES',            'character'),
  c('DOB',                  'Date'),
  c('AGE_AT_DIAGNOSIS',     'numeric'),
  c('BMI_AT_DIAGNOSIS',     'numeric'),
  c('EXON_NO',              'character'), # Split into EXON and EXON_NO, convert to numeric
  c('MUTATION_CODING',      'character'),
  c('MUTATION_AA_CHANGE',   'character'),
  c('MUTATION_INHERITANCE', 'character'),
  c('OHA',                  'character'),
  c('INSULIN',              'logical'),
  c('INSULIN_DELAY',        'numeric')    # TODO unit!?!
)
phenotype_fr <- read.xlsx(paste(phenotype_dir, 'DNAsending_KLUPA_140414_CBellanne.xlsx', sep='/'),
                          sheetIndex = 1,
                          stringsAsFactors = F,
                          colIndex=3:17,
                          colClasses = sapply(columns_fr, function(c)c[2]))
names(phenotype_fr) <- sapply(columns_fr, function(c)c[1])
phenotype_fr$COUNTRY <- 'FR'

phenotype_fr <- subset(phenotype_fr, select=-c(OHA, RELATIVES))

phenotype_fr$PROBAND[is.na(phenotype_fr$PROBAND)] <- F

phenotype_fr$EXON <- substr(toupper(trimws(phenotype_fr$EXON_NO)), 1, 1) == 'E'
phenotype_fr$EXON_NO <- as.numeric(sub('^.* ([0-9]+)', '\\1', phenotype_fr$EXON_NO, perl=T))

# PL - Poland
columns_pl <- list(
  c('SAMPLE_ID',                   'character'),
  c('SEX',                         'character'),
  c('DOB_YEAR',                    'numeric'),
  c('AGE_AT_DIAGNOSIS',            'numeric'),
  c('AGE_AT_CLINICAL_EXAM',        'numeric'),
  c('HEIGHT',                      'numeric'),
  c('WEIGHT',                      'numeric'),
  c('BMI_AT_DIAGNOSIS',            'numeric'),
  c('MUTATION_AA_CHANGE',          'character'),
  c('MUTATION_INHERITANCE_FATHER', 'character'),
  c('MUTATION_INHERITANCE_MOTHER', 'character')
)
phenotype_pl <- read.xlsx(paste(phenotype_dir, 'MODY - microarrays_uzup_MS_2015-11-10.xls', sep='/'),
                          sheetIndex = 1,
                          stringsAsFactors = F,
                          colIndex=1:11,
                          colClasses = sapply(columns_pl, function(c)c[2]))
names(phenotype_pl) <- sapply(columns_pl, function(c)c[1])
phenotype_pl$COUNTRY <- 'PL'

phenotype_pl$FAMILY_ID <- sapply(phenotype_pl$SAMPLE_ID, function (x) unlist(strsplit(x, '-', fixed = T))[1])

phenotype_pl$MUTATION_INHERITANCE <- NA
phenotype_pl$MUTATION_INHERITANCE[phenotype_pl$MUTATION_INHERITANCE_FATHER == 'Y'] <- 'F'
phenotype_pl$MUTATION_INHERITANCE[phenotype_pl$MUTATION_INHERITANCE_MOTHER == 'Y'] <- 'M'

# SK - Slovakia
columns_sk <- list(
  c('FAMILY_ID',                      'character'),
  c('SAMPLE_ID',                      'character'),
  c('SEX',                            'character'),
  c('DOB',                            'Date'),
  c('AGE_AT_DIAGNOSIS',               'numeric'),
  c('AGE_AT_CLINICAL_EXAM',           'numeric'),
  c('HEIGHT',                         'numeric'),
  c('WEIGHT',                         'numeric'),
  c('BMI_AT_CLINICAL_EXAM',           'numeric'),
  c('MUTATION_AA_CHANGE',             'character'),
  c('EXON_NO',                        'numeric'),
  c('MUTATION_CODING',                'character'),
  c('MUTATION_INHERITANCE',           'character'),
  c('MUTATION_INHERITANCE_CONFIRMED', 'character'), # 'Not ...' <- False
  c('INSULIN_START_AGE',              'numeric'),
  c('THERAPY_AT_DIAGNOSIS',           'character'),
  c('THERAPY_AT_REFERRAL',            'character')
)
phenotype_sk <- read.xlsx(paste(phenotype_dir, 'ClinicalDatabase-Slovakia_MS_27022015.xls', sep='/'),
                          sheetName = 'Slovakia',
                          stringsAsFactors = F,
                          colIndex=2:17,
                          colClasses = sapply(columns_sk, function(c)c[2]),
                          na.strings='N/A')
names(phenotype_sk) <- sapply(columns_sk, function(c)c[1])
phenotype_sk$COUNTRY <- 'SK'

phenotype_sk$SAMPLE_ID <- sub(' ', '', phenotype_sk$SAMPLE_ID, fixed=T)

phenotype_sk$MUTATION_INHERITANCE_CONFIRMED[is.na(phenotype_sk$MUTATION_INHERITANCE_CONFIRMED)] <- 'FALSE'
phenotype_sk$MUTATION_INHERITANCE_CONFIRMED <- as.logical(phenotype_sk$MUTATION_INHERITANCE_CONFIRMED)

phenotype_sk$EXON <- T

# Obvious typo
year(phenotype_sk$DOB[phenotype_sk$SAMPLE_ID == 'P350']) <- 1981

# UK - United Kingdom
columns_uk <- list(
  c('SAMPLE_ID',               'character'), # extract FAMILY_ID
  c('SEX',                     'character'),
  c('DOB',                     'Date'),
  c('AGE_AT_DIAGNOSIS',        'numeric'),
  c('AGE_AT_CLINICAL_EXAM',    'numeric'),
  c('STATUS',                  'character'),
  c('RELATIONSHIP_TO_PROBAND', 'character'),
  c('THERAPY_AT_DIAGNOSIS',    'character'), # Common dictionary
  c('THERAPY_AT_REFERRAL',     'character'), # Common dictionary
  c('ETHNIC_ORIGIN',           'character'),
  c('EXON',                    'character'), # Convert to logical
  c('EXON_NO',                 'character'), # Remove non-numeric stuff, convert to numeric
  c('MUTATION_AA_CHANGE',      'character'),
  c('MUTATION_CODING',         'character'),
  c('PROTEIN_EFFECT',          'character'),
  c('P291_FS_INS_C_RESULT',    'character'),
  c('HNF1A_RESULT',            'character'),
  c('GENE_NAME',               'character'),
  c('HEIGHT',                  'numeric'),   # Convert to cm
  c('WEIGHT',                  'numeric'), 
  c('BMI_AT_CLINICAL_EXAM',    'numeric'),
  c('MOTHER_MUTATION',         'character'),
  c('FATHER_MUTATION',         'character'),
  c('MOTHER_DM',               'character'), # Convert to logical
  c('MOTHER_AGE_AT_DIAGNOSIS', 'numeric'),
  c('FATHER_DM',               'character'), # Convert to logical
  c('INSULIN_DELAY',           'character')  # Convert to years
)
phenotype_uk <- read.xlsx(paste(phenotype_dir, 'HNF1A_June2014_Poland_Final-Hattersley.xlsx', sep='/'),
                          sheetName = 'clinical info',
                          stringsAsFactors = F,
                          colIndex=c(1, 5:30),
                          colClasses = sapply(columns_uk, function(c)c[2]))
names(phenotype_uk) <- sapply(columns_uk, function(c)c[1])
phenotype_uk$COUNTRY <- 'UK'

phenotype_uk <- phenotype_uk[grep('IGNORE', phenotype_uk$SAMPLE_ID, invert = T),]

phenotype_uk$FAMILY_ID <- sub('^([0-9]+[A-Z]+)([0-9]+) *.*', '\\1', phenotype_uk$SAMPLE_ID, perl=T)

phenotype_uk$PROBAND <- toupper(trimws(phenotype_uk$STATUS)) == 'PROBAND'

phenotype_uk$EXON <- trimws(toupper(phenotype_uk$EXON)) == 'EXON'
phenotype_uk$EXON_NO <- as.numeric(sub('^.* ([0-9]+)', '\\1', phenotype_uk$EXON_NO, perl=T))

phenotype_uk$HEIGHT <- phenotype_uk$HEIGHT*100

phenotype_uk$MOTHER_DM <- trimws(toupper(phenotype_uk$MOTHER_DM)) == 'YES'
phenotype_uk$FATHER_DM <- trimws(toupper(phenotype_uk$FATHER_DM)) == 'YES'

phenotype_uk$MUTATION_INHERITANCE <- ''
phenotype_uk$MUTATION_INHERITANCE <- NA

# TODO should we round to whole years?!?
phenotype_uk$INSULIN_DELAY <- as.numeric(phenotype_uk$INSULIN_DELAY)/12

# US - USA
columns_us <- list(
  c('FAMILY_ID',            'character'),
  c('SAMPLE_ID',            'character'),
  c('SEX',                  'character'),
  c('AGE_AT_DIAGNOSIS',     'numeric'),
  c('AGE_AT_CLINICAL_EXAM', 'numeric'),
  c('HEIGHT',               'numeric'),
  c('WEIGHT',               'numeric'), 
  c('BMI_AT_CLINICAL_EXAM', 'numeric'),
  c('MUTATION_INHERITANCE', 'character'),
  c('MUTATION_AA_CHANGE',   'character'),
  c('THERAPY',              'character'), # Common dictionary
  c('INSULIN_START_AGE',    'numeric'),
  c('DOB_YEAR',             'numeric')
)
phenotype_us <- read.xlsx(paste(phenotype_dir, 'MODY Aliquots_Phenotype file - Alessandro.xlsx', sep='/'),
                          sheetIndex = 1,
                          stringsAsFactors = F,
                          colIndex=c(1:12, 18),
                          colClasses = sapply(columns_us, function(c)c[2]))
names(phenotype_us) <- sapply(columns_us, function(c)c[1])
phenotype_us$COUNTRY <- 'US'

# Initialize merged phenotype with correct data types
phenotype <- data.frame(
  AGE = numeric(),
  AGE_AT_CLINICAL_EXAM = numeric(),
  AGE_AT_DIAGNOSIS = numeric(),
  AGE_AT_CLINICAL_EXAM = numeric(),
  BMI_AT_CLINICAL_EXAM = numeric(),
  BMI_AT_DIAGNOSIS = numeric(),
  COUNTRY = character(),
  DOB = as.Date(character()),
  DOB_YEAR = numeric(),
  ETHNIC_ORIGIN = character(),
  EXON = logical(),
  EXON_NO = numeric(),
  FAMILY_ID = character(),
  FATHER_DM = logical(),
  FATHER_MUTATION = character(),
  GENE_NAME = character(),
  HEIGHT = numeric(),
  HNF1A_RESULT = character(),
  INSULIN = logical(),
  INSULIN_DELAY = numeric(),
  INSULIN_START_AGE = numeric(),
  MOTHER_AGE_AT_DIAGNOSIS = numeric(),
  MOTHER_DM = logical(),
  MOTHER_MUTATION = character(),
  MUTATION_AA_CHANGE = character(),
  MUTATION_CODING = character(),
  MUTATION_INHERITANCE = character(),
  MUTATION_INHERITANCE_CONFIRMED = logical(),
  P291_FS_INS_C_RESULT = character(),
  PROBAND = logical(),
  PROTEIN_EFFECT = character(),
  RELATIONSHIP_TO_PROBAND = character(),
  SAMPLE_ID = character(),
  SEX = character(),
  STATUS = character(),
  THERAPY = character(),
  THERAPY_AT_DIAGNOSIS = character(),
  THERAPY_AT_REFERRAL = character(),
  WEIGHT = numeric(),
  stringsAsFactors = F
)

# Merge
library(plyr)
phenotype <- ldply(list(phenotype, phenotype_uk, phenotype_cz, phenotype_fr, phenotype_pl, phenotype_sk, phenotype_us), data.frame)

# Final cleanup
phenotype$COUNTRY <- factor(phenotype$COUNTRY) # TODO add labels with full country names
phenotype$FAMILY_ID <- factor(phenotype$FAMILY_ID)
phenotype$SAMPLE_ID <- factor(phenotype$SAMPLE_ID)

phenotype$SEX[phenotype$SEX == 'F'] <- 'Female'
phenotype$SEX[phenotype$SEX == 'M'] <- 'Male'
phenotype$SEX[phenotype$SEX == '0'] <- 'Female'
phenotype$SEX[phenotype$SEX == '1'] <- 'Male'
phenotype$SEX <- factor(phenotype$SEX, levels=c('Male', 'Female'))
phenotype$SEX_PLINK <- as.numeric(phenotype$SEX)

phenotype$DOB_YEAR[is.na(phenotype$DOB_YEAR)] <- year(phenotype$DOB[is.na(phenotype$DOB_YEAR)])
phenotype <- subset(phenotype, select=-DOB)

# TODO clean up MUTATION_AA_CHANGE, use proper HGVS!!!

phenotype$MUTATION_INHERITANCE[phenotype$MUTATION_INHERITANCE == 'N/A'] <- NA
phenotype$MUTATION_INHERITANCE <- toupper(substr(trimws(phenotype$MUTATION_INHERITANCE), 1, 1))
phenotype$MUTATION_INHERITANCE[phenotype$MUTATION_INHERITANCE == 'P'] <- 'F'
phenotype$MUTATION_INHERITANCE <- factor(phenotype$MUTATION_INHERITANCE, labels = c('Both parents', 'De novo', 'Father', 'Mother'))

library(reshape2)
library(ggplot2)
phenotype$BMI_CALCULATED <- phenotype$WEIGHT/((phenotype$HEIGHT/100)^2)
bmi_xcheck <- subset(phenotype, BMI_AT_CLINICAL_EXAM > BMI_CALCULATED+0.2 | BMI_AT_CLINICAL_EXAM < BMI_CALCULATED-0.2)[,c('SAMPLE_ID', 'BMI_AT_CLINICAL_EXAM', 'BMI_CALCULATED')]
bmi_xcheck$diff <- bmi_xcheck$BMI_AT_CLINICAL_EXAM-bmi_xcheck$BMI_CALCULATED
bmi_xcheck[order(bmi_xcheck$diff, decreasing = T),]
qplot(SAMPLE_ID, value, data=melt(bmi_xcheck, id.vars=c(1,4)), color=variable) + coord_flip()

# Expected table columns
columns_final <- c(
  'COUNTRY',              # two letter country code
  'FAMILY_ID',
  'SAMPLE_ID',
  'SEX',                  # Female, Male
  'SEX_PLINK',            # Plink coded (1=male, 2=female) SEX field
  'DOB_YEAR',             # Date of birth, year
  'AGE_AT_DIAGNOSIS',     # [years]
  'AGE_AT_CLINICAL_EXAM', # [years]
  'HEIGHT',               # [cm]
  'WEIGHT',               # [kg]
  'BMI_AT_CLINICAL_EXAM',
  'MUTATION_AA_CHANGE',   # HGVS protein notation with transcript ID
  'MUTATION_INHERITANCE', # Both parents, De novo, Father, Mother
  'INSULIN_START_AGE'     # [years]
)

phenotype_final <- subset(phenotype, select=columns_final)

# Split by country
for (part in split(phenotype_final, phenotype_final$COUNTRY))
  write.xlsx(part, paste0(paste(data_dir, 'phenotype_', sep='/'), part$COUNTRY[1], '.xlsx'), col.names=T, row.names = F, showNA = F)

sample_sheet <- read.csv(paste(data_dir, 'sample_sheet.csv', sep='/'), skip=8)
sample_sheet <- subset(sample_sheet, select=1:6)
names(sample_sheet) <- c('OMICRON_ID', 'SAMPLE_ID', 'PLATE', 'WELL', 'ARRAY_ID', 'ARRAY_POSITION')
sample_sheet$ARRAY_ROW <- as.numeric(substr(as.character(sample_sheet$ARRAY_POSITION), 2, 3))
sample_sheet$ARRAY_COL <- as.numeric(substr(as.character(sample_sheet$ARRAY_POSITION), 5, 6))
sample_sheet$WELL_ROW <- as.factor(substr(as.character(sample_sheet$WELL), 1, 1))
sample_sheet$WELL_COL <- as.numeric(substr(as.character(sample_sheet$WELL), 2, 3))

phenotype_and_sample_sheet <- unique(merge(sample_sheet, phenotype_final, by='SAMPLE_ID', all.x=T))

# Drop phenotypes without sex - they are useless
phenotype_and_sample_sheet <- subset(phenotype_and_sample_sheet, !is.na(SEX_PLINK))
write.table(phenotype_and_sample_sheet, paste(data_dir, 'phenotype.csv', sep='/'), row.names = F, na = '', quote = F, sep='\t')

# Samples that were not genotyped
setdiff(phenotype_final$SAMPLE_ID, sample_sheet$SAMPLE_ID)

# Samples without phenotype information
setdiff(sample_sheet$SAMPLE_ID, phenotype_final$SAMPLE_ID)
