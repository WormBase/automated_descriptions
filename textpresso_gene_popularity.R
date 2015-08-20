#!/usr/bin/env R

p <- function(..., sep='') {
    paste(..., sep=sep, collapse=sep)
}

pdf("textpresso_gene_popularity.pdf")

args<-commandArgs(TRUE)

release <- args[1]

popularity_url <- p("http://textpresso-dev.caltech.edu/concise_descriptions/release/", release, "/c_elegans/gene_lists/textpresso_gene_popularity.text")

popularity <- read.table(url(popularity_url), header = TRUE)

log_data_hist <- hist(log10(popularity$papers[!popularity$papers==0]), plot=FALSE)

plot(log_data_hist,main="Textpresso Popularity of C. elegans Genes",xlab="log(Papers)",col="green",border="blue",ylim=c(0,20000), xlim=c(0,4))
