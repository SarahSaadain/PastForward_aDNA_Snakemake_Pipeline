# Rscript visualize-plotable.R input.plotable output.png
suppressPackageStartupMessages({
  library(tidyverse)  
})

debug=FALSE
if(!debug)
{
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) <2) {
     cat("Please provide an input file and an output file; Usage Rscript visualize-plotable.R input.plotable output.png")
      quit("no", 1)
    }
    # get parameters
    file<-args[1]
    outfile<-args[2]
}else{
  # debug 
  rm(list = ls())
  file<-"/Users/robertkofler/gh/teplotter/test2/mdg1#LTR_Gypsy_te"
  outdir <- tempdir()
}

#
# some parameters; feel free to modify
#
mindeletion=10 # minimum length of the internal deletions
width=8       # plot width 
height=5      # plot height
dpi=300       # plot dpi
#
# end parameters
#

data <- read_tsv(file,col_names = FALSE,cols(.default = col_character()))

# split of coverage
cov <- data |> filter(X3== "cov")
cov <- cov |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,covy=X5) 
cov <- cov |> mutate(pos = as.double(pos),covy= as.double(covy))

# split of ambcoverge
ambcov <- data |> filter(X3== "ambcov")
ambcov <- ambcov |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,ambcovy=X5) 
ambcov <- ambcov |> mutate(pos = as.double(pos),ambcovy= as.double(ambcovy))

# split of snps
snp <- data |> filter(X3=="snp")
snp <- snp |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,refc=X5,base=X6,count=X7) 
snp <- snp |>  mutate(pos = as.double(pos),count= as.double(count))

# split of deletion
deletion <- data |> filter(X3== "del")
deletion <- deletion |> rename(seqid=X1,sampleid=X2,feature=X3,start=X4,end=X5,startcov=X6,endcov=X7,count=X8) 
deletion <- deletion |> mutate(start = as.double(start),end= as.double(end),startcov = as.double(startcov),endcov= as.double(endcov),count= as.double(count))

# split of insertion
insertion <- data |> filter(X3=="ins")
insertion <- insertion |> rename(seqid=X1,sampleid=X2,feature=X3,pos=X4,length=X5,count=X6) 
insertion <- insertion |> mutate(pos = as.double(pos), length= as.double(length), count= as.double(count))

# prepare insertions
# filter min size of insertion
deletion<- deletion |> filter(end-start>mindeletion)
# size of scaling
deletion$scale=log(deletion$count)

theme_set(theme_bw())
plo<-ggplot()+
  geom_polygon(data = cov, mapping = aes(x = pos, y = covy), fill = 'darkgrey', color = 'darkgrey') +
  geom_polygon(data = ambcov, aes(x = pos, y = ambcovy), fill = 'lightgrey', color = 'lightgrey')+
  geom_curve(data = deletion, mapping = aes(x = start, y = startcov, xend = end, yend = endcov, linewidth = scale),  curvature = -0.15, ncp=5,show.legend = FALSE)+
  scale_linewidth(range = c(0.3, 2))+xlab("position") + ylab("coverage")+
  geom_bar(data=snp,aes(x=pos,y=count,fill=base),stat="identity",width=2)+
  geom_bar(data=insertion,aes(x=pos,y=count),stat="identity",color="grey50",width=4)

# legend position and style; 
plo<-plo +
  theme(
    legend.position   = "top",              # or "bottom" if preferred
    legend.direction  = "horizontal",       # default anyway, but explicit
    legend.justification = "center",        # center the legend horizontally
    legend.box        = "horizontal",       # keep keys + labels in one row
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.background = element_rect(fill = "white", colour = "grey80"),  # optional: nicer box
    legend.margin     = margin(2, 0, 2, 0)  # reduce space around legend
  )

# faceting
nseq<-n_distinct(cov$seqid)
nsample<-n_distinct(cov$sampleid)
if (nseq > 1 & nsample>1) {
  plo<-plo+facet_grid(sampleid~seqid,scales = "free_x", space = "free_x")
} else if (nseq>1){
  plo<-plo+facet_grid(~seqid,scales = "free_x", space = "free_x")
}else if (nsample>1){
  plo<-plo+facet_grid(sampleid~.)
}

if(!debug){
  ggsave(outfile, plot = plo, width = width, height = height, dpi = dpi)
  
}else{
  plot(plo)
}


