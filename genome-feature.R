#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse", quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library("tools", quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library("data.table", quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library("GenomicRanges", quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library("ggplot2", quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE))

option_list <- list(
	make_option(c("-d","--density-filter"), help="density < density filter will be kept", default = 1.0),
	make_option(c("-s","--score-filter"), help="socre < score filter will be kept", default = 1.0),
	make_option(c("-l","--range-length"), help="downstream & upstream range length will be cut to bins", default = 5000),
	make_option(c("-b","--bin-size"), help="bin size", default = 100)
)

arguments <- parse_args(OptionParser(usage = "%prog [options] genefile regionfile", option_list = option_list), positional_arguments = 2)
opt <- arguments$options

gDensityFilter <- opt$`density-filter`
gScoreFilter <- opt$`score-filter`
gRangeLength <- opt$`range-length`
gBinSize <- opt$`bin-size`

gGeneFile <- arguments$args[1]
gRegionFile <- arguments$args[2]

if(!file.exists(gGeneFile)){
	stop("gene file \"", gGeneFile ,"\" does not exist.")
}

if(!file.exists(gRegionFile)){
	stop("region file \"", gRegionFile ,"\" does not exist.")
}

filebase <- unlist(strsplit(basename(gRegionFile), ".", fixed = TRUE))
tissue <- filebase[1]
rtype <- tail(filebase, 2)[1]

# main 

genes <- fread(gGeneFile, header = TRUE, sep = "\t")
genes.seg <- data.table(chr = genes$chrom, name = genes$name, name2 = genes$name2, strand = genes$strand, start = genes$txStart, end = genes$txEnd)

seg.range <- gRangeLength
bin.size <- gBinSize
bin.count <- as.integer(seg.range / bin.size)
genes.seg$sStart <- genes.seg$start - seg.range
genes.seg$sEnd <- genes.seg$start
genes.seg$eStart <- genes.seg$end
genes.seg$eEnd <- genes.seg$end + seg.range

genes.seg.bins <- genes.seg[, .(name2, strand, end, bStart = c(head(seq(sStart, sEnd, bin.size), -1), start, head(seq(eStart, eEnd, bin.size), -1)), bEnd = c(tail(seq(sStart, sEnd, bin.size), -1), end, tail(seq(eStart, eEnd, bin.size), -1)), bId = c(-bin.count:-1, 0, 1:bin.count)), by = c("chr", "name", "start")]
genes.seg.bins$bId[genes.seg.bins$strand == '-'] <- -genes.seg.bins$bId[genes.seg.bins$strand == '-']
genes.rg <- GRanges(seqnames = genes.seg.bins$chr, ranges = IRanges(start = genes.seg.bins$bStart, end = genes.seg.bins$bEnd), strand = "+")

regions <- fread(gRegionFile, header = TRUE, sep = "\t")
regions <- regions[meanScore < gScoreFilter & meanDensity < gDensityFilter, ]
regions$center <- as.integer((regions$start + regions$end) / 2)
regions.rg <- GRanges(seqnames = regions$chrom, ranges = IRanges(start = regions$center, end = regions$center), strand = "+")

hits <- findOverlaps(regions.rg, genes.rg, select = "first", type = "within")
genes.seg.bins$hits <- 0
stats <- as.data.frame(table(hits))
genes.seg.bins$hits[stats[, 1]] <- stats[, 2]

total.hits <- sum(genes.seg.bins$hits)
dist <- genes.seg.bins[, .(p = sum(hits) / total.hits), by = "bId"]
dist$tissue <- tussue
out.filename <- paste("gene.pos", tissue, rtype, "csv", sep = ".")
write.csv(dist, out.filename, row.names = FALSE, quote = FALSE)


