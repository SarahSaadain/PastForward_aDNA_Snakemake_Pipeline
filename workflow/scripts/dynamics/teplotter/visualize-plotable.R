# Rscript visualize-plotable.R input.plotable output.png [--log]
suppressPackageStartupMessages({
  library(tidyverse)  
})
debug=FALSE
log_scale=FALSE
log_auto_threshold=NULL

if(!debug)
{
    args <- commandArgs(trailingOnly = TRUE)

    # detect --log flag or --log=N threshold
    log_arg <- args[grepl("^--log", args)]
    args    <- args[!grepl("^--log", args)]   # remove from positional args

    if (length(log_arg) > 0) {
      if (grepl("^--log=", log_arg[1])) {
        # --log=N form: auto-switch if max coverage exceeds N
        log_auto_threshold <- as.numeric(sub("^--log=", "", log_arg[1]))
      } else {
        # plain --log: always use log scale
        log_scale <- TRUE
      }
    }

    if (length(args) < 2) {
      cat("Please provide an input file and an output file;\nUsage: Rscript visualize-plotable.R input.plotable output.png [--log] [--log=N]")
      quit("no", 1)
    }

    file    <- args[1]
    outfile <- args[2]
}else{
  rm(list = ls())
  file      <- "/Users/robertkofler/gh/teplotter/test2/mdg1#LTR_Gypsy_te"
  outdir    <- tempdir()
  log_scale <- FALSE
  log_auto_threshold <- NULL
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

# prepare insertions
# filter min size of insertion
deletion <- deletion |> filter(end - start > mindeletion)
# size of scaling
deletion$scale = log(deletion$count)

# after cov is loaded, resolve auto log threshold
if (!is.null(log_auto_threshold)) {
  if (max(cov$covy, na.rm = TRUE) > log_auto_threshold) {
    log_scale <- TRUE
    message("Auto-switching to log scale: max coverage ", 
            round(max(cov$covy, na.rm=TRUE)), 
            " exceeds threshold ", log_auto_threshold)
  }
}


# log transform BEFORE building the plot
if (log_scale) {
  cov      <- cov      |> mutate(covy     = log10(covy     + 1))
  ambcov   <- ambcov   |> mutate(ambcovy  = log10(ambcovy  + 1))
  deletion <- deletion |> mutate(startcov = log10(startcov + 1),
                                  endcov   = log10(endcov   + 1))
}

theme_set(theme_bw())
plo <- ggplot() +
  geom_polygon(data = cov,    mapping = aes(x = pos, y = covy),    fill = 'darkgrey',  color = 'darkgrey') +
  geom_polygon(data = ambcov, aes(x = pos, y = ambcovy),           fill = 'lightgrey', color = 'lightgrey') +
  geom_curve(data = deletion, mapping = aes(x = start, y = startcov, xend = end, yend = endcov, linewidth = scale), curvature = -0.15, ncp = 5, show.legend = FALSE) +
  scale_linewidth(range = c(0.3, 2)) + xlab("position") + ylab("coverage")

if (log_scale) {
  # stack snps in linear space, then log-transform
  snp_stacked <- snp |>
    group_by(seqid, sampleid, pos) |>
    arrange(base, .by_group = TRUE) |>
    mutate(
      cum_top    = cumsum(count),
      cum_bottom = cumsum(count) - count
    ) |>
    ungroup() |>
    mutate(
      yend   = log10(cum_top    + 1),
      ystart = log10(cum_bottom + 1)
    )

  plo <- plo +
    geom_segment(data = snp_stacked,
                 aes(x = pos, xend = pos, y = ystart, yend = yend, color = base),
                 linewidth = 0.1) +
    geom_segment(data = insertion,
                 aes(x = pos, xend = pos, y = 0, yend = log10(count + 1)),
                 color = "grey50", linewidth = 0.5)

  max_y  <- ceiling(max(cov$covy, na.rm = TRUE))
  fixed_vals   <- c(0, 1, 2, 5, 10, 50, 100, 1000, 5000, 10000, 100000, 1000000)
  fixed_breaks <- log10(fixed_vals + 1)
  
  # only keep breaks up to max_y
  keep   <- fixed_breaks <= max_y
  breaks <- fixed_breaks[keep]
  labels <- sapply(fixed_vals[keep], \(v) {
    if      (v == 0)    "0"
    else if (v >= 1e6)  paste0(round(v/1e6, 1), "M")
    else if (v >= 1e3)  paste0(round(v/1e3, 1), "K")
    else                as.character(v)
  })

  plo <- plo +
    scale_y_continuous(breaks = breaks, labels = labels, 
                       limits = c(0, max_y), 
                       expand = expansion(mult = c(0, 0.01))) +
    ylab("coverage (log10)") 

} else {
  plo <- plo +
    geom_bar(data = snp,       aes(x = pos, y = count, fill = base), stat = "identity", width = 2) +
    geom_bar(data = insertion, aes(x = pos, y = count),              stat = "identity", color = "grey50", width = 4)
}

# legend position and style; 
plo<-plo +
  theme(
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.justification = "center",
    legend.box        = "horizontal",
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 10),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    legend.margin     = margin(2, 0, 2, 0)
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

nplots <- nseq + nsample

# change plot height dynamically
if (nplots>1) {
  height<-2*(nplots)
}

if(!debug){
  ggsave(outfile, plot = plo, width = width, height = height, dpi = dpi)
}else{
  plot(plo)
}