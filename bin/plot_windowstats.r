#!/usr/bin/env Rscript

library("optparse")

args <- commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option(c("-b", "--bed"),
              type = "character",
              help = "BED file with sample-wise window stats"),
  make_option(c("-o", "--outfile"),
              type = "character",
              help = "Output PDF file")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Plots sample-wise window stats based on BED input")
opt <- parse_args(opt_parser, args=args)


data <- read.table(opt$bed, header=TRUE, sep="\t", comment.char="")
colnames(data) <- sub("^X\\.?", "", colnames(data))
samples <- colnames(data)[9:11]
samples <- sub("_.*$", "", samples)

pdf(file=opt$outfile, width=12, height=9)
par(mfrow=c(3,1), mar=c(2.5, 4.4, 3, 2.1), oma=c(2.5, 0, 3, 0), cex=0.7)

cols <- c("black", "blue", "red")
chroms <- levels(factor(data$chrom))
#chroms <- grep("chr", chroms, value=TRUE)

for (chrom in chroms){

  chrom_data <- data[data$chrom==chrom,]

  print(paste0("Working on chromosome ", chrom, " ..."))
  
  plot(chrom_data$start, chrom_data$gc_content, type="l", col="black", xlim=c(0, max(chrom_data$start)), ylim=c(0.2, 0.8), xlab=NA, ylab=NA)
  mtext(side=2, text="GC content", cex=0.8, line=3)
  mtext(side=1, text="Position [bp]", cex=0.8, line=1.2, outer=TRUE)
  title(chrom, cex.main=1.5, line=1.5, outer=FALSE)

  plot(NULL, type='n', xlim=c(0, max(chrom_data$start)), ylim=c(0, 2.5), xlab=NA, ylab=NA)
  abline(h=0.5, lwd=1, lty=5)
  abline(h=1.0, lwd=1, lty=5)
  for (i in 1:length(samples)) {
    lines(chrom_data$start, chrom_data[, paste0(samples[i], "_corr")], type='l', col=cols[i])
  }
  mtext(side=2, text="Depth", cex=0.8, line=3)
  mtext(side=1, text="Position [bp]", cex=0.8, line=1.2, outer=TRUE)
  legend('topright', samples, cex=1.2, lty=rep.int(1, length(samples)), lwd=rep.int(2, length(samples)), col=cols, ncol=1)

  plot(NULL, type='n', xlim=c(0, max(chrom_data$start)), ylim=c(0, 0.01), xlab=NA, ylab=NA)
  for (i in 1:length(samples)) {
    lines(chrom_data$start, chrom_data[, paste0(samples[i], "_het")] / chrom_data[, paste0(samples[i], "_val")], type='l', col=cols[i])
  }
  mtext(side=2, text="Heterozygosity", cex=0.8, line=3)
  mtext(side=1, text="Position [bp]", cex=0.8, line=1.2, outer=TRUE)
  legend('topright', samples, cex=1.2, lty=rep.int(1, length(samples)), lwd=rep.int(2, length(samples)), col=cols, ncol=1)

}

dev.off()
