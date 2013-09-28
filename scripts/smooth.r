# the program is used to smooth the read coverage of plus and minus strand of each contig, it is called in by chipdenovo.py
library(seqinr)
args <- commandArgs(trailingOnly = TRUE)
fp <- read.fasta(args[1], as.string=TRUE)
fm <- read.fasta(args[2], as.string=TRUE)
for (i in 1:length(fp)){
	covp <- as.numeric(strsplit(fp[[i]][1], " ")[[1]])
	covm <- as.numeric(strsplit(fm[[i]][1], " ")[[1]])
	x <- 1:length(covp)
	x.points <- x
	covp.sm <- ksmooth(x,covp,'normal',bandwidth=100)$y
	covm.sm <- ksmooth(x,covm,'normal',bandwidth=100)$y
	fp[[i]][1] <- paste(covp.sm, collapse=' ')
	fm[[i]][1] <- paste(covm.sm, collapse=' ')
}
write.fasta(sequences=fp,names=names(fp),file.out=paste(args[1],'.smoothed.fa',sep=''))
write.fasta(sequences=fm,names=names(fm),file.out=paste(args[2],'.smoothed.fa',sep=''))
