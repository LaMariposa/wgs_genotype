setwd("tlens")

files=list.files(pattern='tlen.txt')

pdf("tlen_dist.pdf",height=20, width=20)
par(mfrow=c(5,5))

for (infile in files)
  {
    print(infile)
    tlen=read.table(infile, header=F)
    prefix=substr(infile, 1,(nchar(infile)-9))
    tlensub=subset(tlen, V1>0 & V1<1000)
    hist(tlensub$V1, main=prefix)
  }

dev.off()

