library(reshape2)
library(ggplot2)
library(plyr)
library(printr)

# Load final phenotype
phenotype <- read.csv('/data/phenotype.csv', sep='\t')

#' # Metadata statistical analysis
#' ## Number of samples
nrow(phenotype)

#' ## Count per country
table(phenotype$COUNTRY)

table(phenotype$COUNTRY, phenotype$SEX)

#' ## Plots
p <- subset(phenotype, select=c(1, 11, 13, 15:ncol(phenotype)))
d_ply(melt(p), .(variable), function (phenotype){
  print(qplot(variable, value, data=phenotype, color=COUNTRY, facets = ~SEX, geom='boxplot'))
})

#' ## Basic stats per sample type and group
ddply(phenotype, .(COUNTRY), function (md) {
  ldply(names(md)[c(15:20, 23)], function (col) { 
    data.frame(name=col,
               min=min(md[,col], na.rm=T),
               max=max(md[,col], na.rm=T),
               median=median(md[,col], na.rm=T),
               mean=mean(md[,col], na.rm=T),
               sd=sd(md[,col], na.rm=T))
  })
})

#' ## Normality tests
# http://stackoverflow.com/a/4357932
qqplot.data <- function (vec)
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  d <- data.frame(y = vec)
  ggplot(d, aes(sample = y)) + stat_qq() + geom_abline(slope = slope, intercept = int)
}

diagnose_normality <- function (metadata) {
  md <- metadata[,c(15:20, 23)]
  cols <- names(md)[sapply(md, function (col){!all(is.na(col))})]
  for (name in cols) {
    if (!all(is.na(md[,name]))) {
      print(qqplot.data(md[,name]) + labs(title=paste('QQ plot for', name, 'in', metadata$COUNTRY[1])))
    }
  }
  try({
    print(as.character(metadata$COUNTRY[1]))
    tests <- ldply(cols, function (col) {
      t <- shapiro.test(md[,col])
      data.frame(name=col, statistic=t$statistic, p.value=t$p.value)
    })
    print(tests[order(tests$p.value),])
  })
}

#' ### Normality
d_ply(phenotype, .(COUNTRY), function (md) { diagnose_normality(md) })
