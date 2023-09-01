library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, cyto.exclude<-c("chrX", "chrY"),tracks.inside=10, tracks.outside=0 )
out.file<-"RCircosDemoHumanGenome.pdf"#output file name
pdf(file=out.file,height=8,width=8,compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data) 
name.col <- 4 
side <- "in" 
track.num <- 1 
data = read.csv('mac-m6a.csv')
head(data)
RCircos.Gene.Connector.Plot(data,
                              + track.num, side)
track.num <- 2
RCircos.Gene.Name.Plot(data,
                         + name.col,track.num, side)#gene name
dev.off()