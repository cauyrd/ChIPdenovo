# the program is used to smooth the read coverage of plus and minus strand of each contig, it is called in by chipdenovo.py
library(Biostrings)
args <- commandArgs(trailingOnly = TRUE)
fp <- readFASTA(args[1])
fm <- readFASTA(args[2])
for (i in 1:length(fp)){
	covp <- as.numeric(strsplit(fp[[i]]$seq, " ")[[1]])
	covm <- as.numeric(strsplit(fm[[i]]$seq, " ")[[1]])
	x <- 1:length(covp)
	x.points <- x
	covp.sm <- ksmooth(x,covp,'normal',bandwidth=100)$y
	covm.sm <- ksmooth(x,covm,'normal',bandwidth=100)$y
	fp[[i]]$seq <- paste(covp.sm, collapse=' ')
	fm[[i]]$seq <- paste(covm.sm, collapse=' ')
}
writeFASTA(fp,file=paste(args[1],'.smoothed.fa',sep=''))
writeFASTA(fm,file=paste(args[2],'.smoothed.fa',sep=''))
