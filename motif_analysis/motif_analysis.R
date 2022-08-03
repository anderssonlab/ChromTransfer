library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)
library(BiocParallel)
library(ggplot2)
library(viridis)

theme_set(theme_bw(base_size=10))

## Parameters
motif.thres <- 500
motif.cnt.thres <- 100
DHS.width <- 600

classA.file <- "grad-x-input_data/HEPG2_pos_pred_finetuned.grad-x-input"
classB.file <- "grad-x-input_data/HEPG2_pos_pred_pretrained.grad-x-input"

classA.name <- "HEPG2_finetuned"
classB.name <- "HEPG2_pretrained"

project.name <- "HEPG2_finetuned_vs_HEPG2_pretrained"
JASPAR.file <- "JASPAR2022_hg38.bb"
workers <- 30

register(MulticoreParam(workers=workers))

plot.dir <- paste(project.name,"_",motif.thres,"_plots/",sep="")
if (!file.exists(plot.dir))
    dir.create(plot.dir)

olap.dir <- paste(project.name,"_",motif.thres,"_motif_olaps/",sep="")
if (!file.exists(olap.dir))
    dir.create(olap.dir)

data.dir <- paste(project.name,"_",motif.thres,"_data/",sep="")
if (!file.exists(data.dir))
    dir.create(data.dir)


## Read in data

message("reading feature attribution data")

classA <- scan(classA.file,sep="\n",what=character(0),quiet=TRUE)
n <- gsub(">","",classA[seq(1,length(classA),by=3)])
classA <- t(sapply(classA[seq(3,length(classA),by=3)],function(x) as.numeric(strsplit(x," ")[[1]])))
rownames(classA) <- n
colnames(classA) <- NULL

classB <- scan(classB.file,sep="\n",what=character(0),quiet=TRUE)
n <- gsub(">","",classB[seq(1,length(classB),by=3)])
classB <- t(sapply(classB[seq(3,length(classB),by=3)],function(x) as.numeric(strsplit(x," ")[[1]])))
rownames(classB) <- n
colnames(classB) <- NULL

## Plot mean gradient*input

df.means <- data.frame(position=1:DHS.width,GI=c(colMeans(classA),colMeans(classB)),model=rep(c(classA.name,classB.name),each=DHS.width))
ggplot(df.means,aes(position,GI,color=model)) + geom_line(size=1) + ylab("average gradient x input") + scale_colour_viridis(discrete=TRUE)
ggsave(paste(plot.dir,"mean_gradient_x_input_vs_model.pdf",sep=""),width=7,height=4)


## Build GRanges object from DHS positions

classA.DHS.gr <- GRanges(seqnames=as.character(sapply(rownames(classA), function(n) strsplit(n,":")[[1]][1])),
                         ranges=IRanges(as.numeric(sapply(rownames(classA), function(n) strsplit(strsplit(n,":")[[1]][2],"-")[[1]][1])),
                                        as.numeric(sapply(rownames(classA), function(n) strsplit(strsplit(n,":")[[1]][2],"-")[[1]][2]))-1,
                                        names=rownames(classA)),
                         strand=strand("*"))

classB.DHS.gr <- GRanges(seqnames=as.character(sapply(rownames(classB), function(n) strsplit(n,":")[[1]][1])),
                         ranges=IRanges(as.numeric(sapply(rownames(classB), function(n) strsplit(strsplit(n,":")[[1]][2],"-")[[1]][1])),
                                        as.numeric(sapply(rownames(classB), function(n) strsplit(strsplit(n,":")[[1]][2],"-")[[1]][2]))-1,
                                        names=rownames(classB)),
                         strand=strand("*"))

## Associate predicted TF binding sites (JASPAR 2022) with model GI scores

message("associating predicted TFBSs with model feature attribution scores")

bbf <- BigBedFile(JASPAR.file)

if (file.exists("motif_GI_score.RData")) {
    load(paste(data.dir,"motif_GI_score.RData",sep=""))
} else {

    classA.motif.GI.score <- mclapply(rownames(classA), function(n) {
        gr <- classA.DHS.gr[n]
        m <- import.bb(bbf, selection=BigBedSelection(gr, colnames = c("name","score")))
        strand(m) <- "*"
        m <- subset(m,score >= motif.thres)
        mn <- paste(as.character(ranges(m)),m$name,sep="_")
        umn <- unique(mn)
        m <- m[match(umn,mn),]
        m.gr <- mapToTranscripts(m,gr,ignore.strand=TRUE)
        s <- start(ranges(m.gr))
        e <- end(ranges(m.gr))

        if (length(s)==0)
            return(GRanges())

        xHits <- m.gr$xHits
        mcols(m.gr) <- NULL
        m.gr$motif <- m$name[xHits]
        m.gr$score <- m$score[xHits]

        f <- as.numeric(classA[n,1:ncol(classA)])
        fd <- sapply(1:length(s), function(j) c(max(f[s[j]:e[j]]),mean(f[s[j]:e[j]])))

        m.gr$max.GI <- fd[1,]
        m.gr$mean.GI <- fd[2,]

        m.gr
    },mc.cores=workers)
    names(classA.motif.GI.score) <- rownames(classA)

    classB.motif.GI.score <- mclapply(rownames(classB), function(n) {

        gr <- classB.DHS.gr[n]
        m <- import.bb(bbf, selection=BigBedSelection(gr, colnames = c("name","score")))
        strand(m) <- "*"
        m <- subset(m,score >= motif.thres)
        mn <- paste(as.character(ranges(m)),m$name,sep="_")
        umn <- unique(mn)
        m <- m[match(umn,mn),]
        m.gr <- mapToTranscripts(m,gr,ignore.strand=TRUE)
        s <- start(ranges(m.gr))
        e <- end(ranges(m.gr))

        if (length(s)==0)
            return(GRanges())

        xHits <- m.gr$xHits
        mcols(m.gr) <- NULL
        m.gr$motif <- m$name[xHits]
        m.gr$score <- m$score[xHits]

        p <- as.numeric(classB[n,1:ncol(classB)])
        pd <- sapply(1:length(s), function(j) c(max(p[s[j]:e[j]]),mean(p[s[j]:e[j]])))

        m.gr$max.GI <- pd[1,]
        m.gr$mean.GI <- pd[2,]

        m.gr
    },mc.cores=workers)
    names(classB.motif.GI.score) <- rownames(classB)

    save(classA.motif.GI.score, classB.motif.GI.score,file=paste(data.dir,"motif_GI_score.RData",sep=""))
}

## Combine motif gradient * input scores into one large data frame

if (file.exists("motif_GI_score_df.RData")) {
    load(paste(data.dir,"motif_GI_score_df.RData",sep=""))
} else {
    classA.motif.GI.score.df <- do.call("rbind",
                                           mclapply(classA.motif.GI.score, function(gr) {
                                               if (length(gr)==0)
                                                   return(data.frame())
                                               d <- as.data.frame(mcols(gr)[,-1])
                                               m <- mcols(gr)[,1]
                                               chr <- as.character(seqnames(classA.DHS.gr[as.character(seqnames(gr))]))
                                               start <- start(classA.DHS.gr[as.character(seqnames(gr))])
                                               d$chr <- chr
                                               d$start <- start+start(ranges(gr))-1
                                               d$end <- start+end(ranges(gr))-1
                                               d$pos <- paste(chr,":",d$start,"-",d$end,sep="")
                                               d$region <- as.character(seqnames(gr))
                                               d$motif <- m
                                               d
                                           },mc.cores=workers))

    classB.motif.GI.score.df <- do.call("rbind",
                                            mclapply(classB.motif.GI.score, function(gr) {
                                                if (length(gr)==0)
                                                    return(data.frame())
                                                d <- as.data.frame(mcols(gr)[,-1])
                                                m <- mcols(gr)[,1]
                                                chr <- as.character(seqnames(classB.DHS.gr[as.character(seqnames(gr))]))
                                                start <- start(classB.DHS.gr[as.character(seqnames(gr))])
                                                d$chr <- chr
                                                d$start <- start+start(ranges(gr))-1
                                                d$end <- start+end(ranges(gr))-1
                                                d$pos <- paste(chr,":",d$start,"-",d$end,sep="")
                                                d$region <- as.character(seqnames(gr))
                                                d$motif <- m
                                                d
                                            },mc.cores=workers))

    motif.GI.score.df <- merge(classA.motif.GI.score.df, classB.motif.GI.score.df, by=c("region","motif","pos","chr","start","end"), all=TRUE, suffixes = c(".classA",".classB"))

    rm(classA.motif.GI.score.df, classB.motif.GI.score.df)

    ## Calculate differences in GI scores between models

    motif.GI.score.df$max.GI.diff <- motif.GI.score.df$max.GI.classA - motif.GI.score.df$max.GI.classB
    motif.GI.score.df$mean.GI.diff <- motif.GI.score.df$mean.GI.classA - motif.GI.score.df$mean.GI.classB

    motif.GI.score.df$ID <- seq(1:nrow(motif.GI.score.df))
    motif.GI.score.df.classA <- subset(motif.GI.score.df, !is.na(max.GI.classA))
    motif.GI.score.df.classB <- subset(motif.GI.score.df, !is.na(max.GI.classB))
    motif.GI.score.df.diff <- subset(motif.GI.score.df, !is.na(max.GI.diff))

    save(motif.GI.score.df,motif.GI.score.df.classA,motif.GI.score.df.classB,motif.GI.score.df.diff,file=paste(data.dir,"motif_GI_score_df.RData",sep=""))
}

classA.motif.cnt <- table(subset(motif.GI.score.df,!is.na(max.GI.classA))$motif)
classB.motif.cnt <- table(subset(motif.GI.score.df,!is.na(max.GI.classB))$motif)
motifs <- union(names(classA.motif.cnt[classA.motif.cnt>=motif.cnt.thres]),
                names(classB.motif.cnt[classB.motif.cnt>=motif.cnt.thres]))

## Plot motif score distributions
df <- subset(data.frame(GI=c(motif.GI.score.df$max.GI.classA,motif.GI.score.df$max.GI.classB),model=rep(c(classA.name,classB.name),each=nrow(motif.GI.score.df))), !is.na(GI))
ggplot(df,aes(GI,color=model)) + geom_density() + xlab("gradient x input") + scale_colour_viridis(discrete=TRUE)
ggsave(paste(plot.dir,"gradient_x_input_vs_model_density.pdf",sep=""),width=5,height=3)

ggplot(subset(motif.GI.score.df, !is.na(max.GI.diff)),aes(max.GI.diff)) + geom_density() + xlab(paste("gradient x input diff.\n(",classA.name," - ",classB.name,")",sep=""))
ggsave(paste(plot.dir,"gradient_x_input_model_difference_density.pdf",sep=""),width=5,height=3)

## Calculate overlaps between motif sets

message("finding overlaps between predicted TFBSs")

motif.gr <- GRanges(seqnames=motif.GI.score.df$chr,
                    ranges=IRanges(motif.GI.score.df$start,
                                   motif.GI.score.df$end,
                                   names=motif.GI.score.df$pos),
                    strand=strand("*"))
motif.gr$motif <- motif.GI.score.df$motif

vec <- 1:length(motif.gr)

safe_motif_name <- function(n) {
    gsub("^[.]+$","_",gsub("[/\\?<>\\:*|\":]","_",n))
}

olaps <- mclapply(motifs, function(n) {
    file <- paste(olap.dir,safe_motif_name(n),sep="")
    if (!file.exists(file)) {
        m <- which(motif.gr$motif==n)
        olap <- vec[-m][unique(queryHits(findOverlaps(motif.gr[-m],motif.gr[m,])))]
        cat(olap,sep="\n",file=file)
        length(olap)
    }
},mc.cores=workers)

## Calcualte rank-based enrichments of motifs

message("calculating rank-based enrichments of motif feature attribution scores")

if (file.exists("motif_ks.RData")) {
    load(paste(data.dir,"motif_ks.RData",sep=""))
} else {

    rf <- frank(motif.GI.score.df.classA$max.GI.classA)
    rp <- frank(motif.GI.score.df.classB$max.GI.classB)
    rd <- frank(motif.GI.score.df.diff$max.GI.diff)

    motif.ks <- mclapply(motifs, function(n) {

        m.f <- which(motif.GI.score.df.classA$motif==n)
        m.p <- which(motif.GI.score.df.classB$motif==n)
        m.d <- which(motif.GI.score.df.diff$motif==n)

        o <- scan(paste(olap.dir,safe_motif_name(n),sep=""),what=numeric(0),sep="\n",quiet=TRUE)
        o.f <- match(intersect(o,motif.GI.score.df.classA$ID), motif.GI.score.df.classA$ID)
        o.p <- match(intersect(o,motif.GI.score.df.classB$ID), motif.GI.score.df.classB$ID)
        o.d <- match(intersect(o,motif.GI.score.df.diff$ID), motif.GI.score.df.diff$ID)

        fm <- motif.GI.score.df.classA$max.GI.classA[m.f]
        fnm <- motif.GI.score.df.classA$max.GI.classA[-c(m.f,o.f)]
        pm <- motif.GI.score.df.classB$max.GI.classB[m.p]
        pnm <- motif.GI.score.df.classB$max.GI.classB[-c(m.p,o.p)]
        dm <- motif.GI.score.df.diff$max.GI.diff[m.d]
        dnm <- motif.GI.score.df.diff$max.GI.diff[-c(m.d,o.d)]

        res1 <- ks.test(fm,fnm)
        res2 <- ks.test(pm,pnm)
        res3 <- try(ks.test(dm,dnm),silent=TRUE)
        if (class(res3) == "try-error")
            res3 <- list(statistic=NA, p.value=NA)

        sf <- ifelse(mean(rf[m.f])<mean(rf[-c(m.f,o.f)]),-1,1)
        sp <- ifelse(mean(rp[m.p])<mean(rp[-c(m.p,o.p)]),-1,1)
        sd <- ifelse(mean(rd[m.d])<mean(rd[-c(m.d,o.d)]),-1,1)

        c(sf*res1$statistic, res1$p.value, sp*res2$statistic, res2$p.value, sd*res3$statistic, res3$p.value)
    },mc.cores=workers)

    names(motif.ks) <- motifs
    save(motif.ks, file=paste(data.dir,"motif_ks.RData",sep=""))
}

if (file.exists("motif_ks_df.RData")) {
    load(paste(data.dir,"motif_ks_df.RData",sep=""))
} else {

    ks.df <- data.frame(motif=motifs,
                        classA.D=sapply(motif.ks,function(x) x[1]),
                        classA.pval=sapply(motif.ks,function(x) x[2]),
                        classB.D=sapply(motif.ks,function(x) x[3]),
                        classB.pval=sapply(motif.ks,function(x) x[4]),
                        diff.D=sapply(motif.ks,function(x) x[5]),
                        diff.pval=sapply(motif.ks,function(x) x[6]))

    ks.df$classA.fdr <- p.adjust(ks.df$classA.pval,method="BH")
    ks.df$classB.fdr <- p.adjust(ks.df$classB.pval,method="BH")
    ks.df$diff.fdr <- p.adjust(ks.df$diff.pval,method="BH")
    ks.df$sig <- ks.df$classB.fdr < 0.001 | ks.df$classA.fdr < 0.001

    s <- which(ks.df$sig)
    ks.df$classA.rank <- NA
    ks.df$classB.rank <- NA
    ks.df$diff.rank <- NA
    ks.df[s,]$classA.rank <- rank(ks.df[s,]$classA.D)
    ks.df[s,]$classB.rank <- rank(ks.df[s,]$classB.D)
    ks.df[s,]$diff.rank <- rank(ks.df[s,]$diff.D)

    save(ks.df, file=paste(data.dir,"motif_ks_df.RData",sep=""))
}

## Plot KS results

message("plotting results")

lim <- c(-1,1)*max(abs(ks.df$diff.D))

focus <- unique(c(head(subset(ks.df,sig)[order(subset(ks.df,sig)$classA.D),"motif"],5),
                  head(subset(ks.df,sig)[order(subset(ks.df,sig)$classB.D),"motif"],5),
                  tail(subset(ks.df,sig)[order(subset(ks.df,sig)$classA.D),"motif"],10),
                  tail(subset(ks.df,sig)[order(subset(ks.df,sig)$classB.D),"motif"],10),
                  tail(subset(ks.df,sig)[order(subset(ks.df,sig)$diff.D),"motif"],10),
                  head(subset(ks.df,sig)[order(subset(ks.df,sig)$diff.D),"motif"],10),
                  head(subset(ks.df,sig & classB.D>0)[order(subset(ks.df,sig & classB.D>0)$classA.D),"motif"],5),
                  tail(subset(ks.df,sig & classB.D<0)[order(subset(ks.df,sig & classB.D<0)$classA.D),"motif"],5),
                  head(subset(ks.df,sig & classA.D>0)[order(subset(ks.df,sig & classA.D>0)$classB.D),"motif"],5),
                  tail(subset(ks.df,sig & classA.D<0)[order(subset(ks.df,sig & classA.D<0)$classB.D),"motif"],5)))

ggplot(subset(ks.df,sig), aes(classA.D, classB.D)) +
    geom_vline(xintercept=0,color="lightgrey",size=0.5) +
    geom_hline(yintercept=0,color="lightgrey",size=0.5) +
    geom_point(aes(color=diff.D),size=1) +
    geom_text_repel(data=subset(ks.df,motif %in% focus),aes(label=motif), size=1.5, max.overlaps=25,
                    min.segment.length=unit(0, 'lines'), segment.size=0.2, segment.color="grey50") +
    coord_fixed() +
    scale_colour_distiller(type="div",palette="RdBu",limits=lim) +
    labs(colour = paste("K-S D (",classA.name," -\n",classB.name,")",sep="")) +
    xlab(paste("K-S D (",classA.name," model)",sep="")) +
    ylab(paste("K-S D (",classB.name," model)",sep=""))
ggsave(paste(plot.dir,"motif_gradient_x_input_KS_D_",classB.name,"_vs_",classA.name,"_scatter.pdf",sep=""),width=8,height=4)

tmp <- subset(ks.df,sig)

ggplot(tmp, aes(classA.rank,classA.D)) + geom_point(size=1) +
    geom_text_repel(data=subset(tmp,classA.rank <= 10),
                    aes(label=motif), size=2, max.overlaps=20,
                    nudge_x=-100-subset(tmp,classA.rank <= 10)$classA.rank,
                    min.segment.length=unit(0, 'lines'), segment.size=0.2, segment.color="grey50",
                    hjust=1, direction="y") +
    geom_text_repel(data=subset(tmp,classA.rank >= max(tmp$classA.rank)-9),
                    aes(label=motif), size=2, max.overlaps=20,
                    nudge_x=max(tmp$classA.rank)+100-subset(tmp,classA.rank >= max(tmp$classA.rank)-9)$classA.rank,
                    min.segment.length=unit(0, 'lines'), segment.size=0.2, segment.color="grey50",
                    hjust=0, direction="y") +
    ylab(paste("K-S D (",classA.name," model)",sep="")) +
    xlab(paste("rank (",classA.name," K-S D)",sep="")) +
    scale_x_continuous(breaks=seq(0,200*ceiling(max(tmp$classA.rank)/200),by=200),
                       limits=c(-300,max(tmp$classA.rank)+300))
ggsave(paste(plot.dir,"motif_gradient_x_input_KS_D_",classA.name,".pdf",sep=""),width=4,height=4)

ggplot(tmp, aes(classB.rank,classB.D)) + geom_point(size=1) +
    geom_text_repel(data=subset(tmp,classB.rank <= 10),
                    aes(label=motif), size=2, max.overlaps=20,
                    nudge_x=-100-subset(tmp,classB.rank <= 10)$classB.rank,
                    min.segment.length=unit(0, 'lines'), segment.size=0.2, segment.color="grey50",
                    hjust=1, direction="y") +
    geom_text_repel(data=subset(tmp,classB.rank >= max(tmp$classB.rank)-9),
                    aes(label=motif), size=2, max.overlaps=20,
                    nudge_x=max(tmp$classB.rank)+100-subset(tmp,classB.rank >= max(tmp$classB.rank)-9)$classB.rank,
                    min.segment.length=unit(0, 'lines'), segment.size=0.2, segment.color="grey50",
                    hjust=0, direction="y") +
    ylab(paste("K-S D (",classB.name," model)",sep="")) +
    xlab(paste("rank (",classB.name," K-S D)",sep="")) +
    scale_x_continuous(breaks=seq(0,200*ceiling(max(tmp$classB.rank)/200),by=200),
                       limits=c(-300,max(tmp$classB.rank)+300))
ggsave(paste(plot.dir,"motif_gradient_x_input_KS_D_",classB.name,".pdf",sep=""),width=4,height=4)


diff.lim <- c(min(motif.GI.score.df$max.GI.diff),max(motif.GI.score.df$max.GI.diff))
GI.lim <- c(min(c(motif.GI.score.df$max.GI.classA,motif.GI.score.df$max.GI.classB)),
            max(c(motif.GI.score.df$max.GI.classA,motif.GI.score.df$max.GI.classB)))

ecdf.data.2 <- function(n) {
    m.f <- which(motif.GI.score.df.classA$motif==n)
    m.p <- which(motif.GI.score.df.classB$motif==n)
    m.d <- which(motif.GI.score.df.diff$motif==n)

    o <- scan(paste(olap.dir,safe_motif_name(n),sep=""),what=numeric(0),sep="\n",quiet=TRUE)
    o.f <- match(intersect(o,motif.GI.score.df.classA$ID), motif.GI.score.df.classA$ID)
    o.p <- match(intersect(o,motif.GI.score.df.classB$ID), motif.GI.score.df.classB$ID)
    o.d <- match(intersect(o,motif.GI.score.df.diff$ID), motif.GI.score.df.diff$ID)

    r.f <- c(m.f,o.f)
    r.p <- c(m.p,o.p)
    r.d <- c(m.d,o.d)

    fm <- motif.GI.score.df.classA$max.GI.classA[m.f]
    efm <- ecdf(fm)
    fm <- unique(fm)

    fnm <- motif.GI.score.df.classA$max.GI.classA[-r.f]
    efnm <- ecdf(fnm)
    fnm <- unique(fnm)

    pm <- motif.GI.score.df.classB$max.GI.classB[m.p]
    epm <- ecdf(pm)
    pm <- unique(pm)

    pnm <- motif.GI.score.df.classB$max.GI.classB[-r.p]
    epnm <- ecdf(pnm)
    pnm <- unique(pnm)

    dm <- motif.GI.score.df.diff$max.GI.diff[m.d]
    edm <- function(x) rep(NA,length(x))
    if (length(dm)==0)
        dm <- NA
    else {
        edm <- ecdf(dm)
        dm <- unique(dm)
    }

    dnm <- motif.GI.score.df.diff$max.GI.diff[-r.d]
    ednm <- function(x) rep(NA,length(x))
    if (length(dnm)==0)
        dnm <- NA
    else {
        ednm <- ecdf(dnm)
        dnm <- unique(dnm)
    }

    list(data.frame(max.GI.diff=c(dm,dnm),
                    ECDF=c(edm(dm),ednm(dnm)),
                    motif=c(rep(n,length(dm)),rep(paste("not",n),length(dnm)))),
         data.frame(GI=c(fm,fnm,pm,pnm),
                    ECDF=c(efm(fm),efnm(fnm),epm(pm),epnm(pnm)),
                    motif=factor(c(rep(n,length(fm)),rep(paste("not",n),length(fnm)),rep(n,length(pm)),rep(paste("not",n),length(pnm))),levels=c(n, paste("not",n))),
                    model=factor(c(rep(classA.name,length(c(fm,fnm))),rep(classB.name,length(c(pm,pnm)))),levels=c(classB.name,classA.name))))
}

tmp <- subset(ks.df,sig)
focus.motifs <- unique(c(subset(tmp,classA.rank <= 5)$motif,
                         subset(tmp,classA.rank >= max(tmp$classA.rank)-4)$motif,
                         subset(tmp,classB.rank <= 5)$motif,
                         subset(tmp,classB.rank >= max(tmp$classA.rank)-4)$motif,
                         subset(tmp,diff.rank <= 10)$motif,
                         subset(tmp,diff.rank >= max(tmp$classA.rank)-9)$motif))

for (m in focus.motifs) {
    tmp <- ecdf.data.2(m)
    ggplot(tmp[[2]],aes(GI,ECDF,colour=paste(motif,model))) + geom_line(size=0.5) + xlab("gradient_x_input") + ylab("ECDF") + labs(colour="") + xlim(GI.lim[1],GI.lim[2]) + scale_colour_viridis(discrete=TRUE)
    ggsave(paste(plot.dir,safe_motif_name(m),"_ECDF_gradient_x_input.pdf",sep=""),width=4,height=3)

    if (!is.na(tmp[[1]]$max.GI.diff[1])) {
        ggplot(tmp[[1]],aes(max.GI.diff,ECDF,color=motif)) + geom_line(size=0.5) +
            xlab(paste("gradient x input diff.\n(",classA.name," - ",classB.name," model)",sep="")) + ylab("ECDF") + xlim(diff.lim[1],diff.lim[2]) + scale_colour_viridis(discrete=TRUE)
        ggsave(paste(plot.dir,safe_motif_name(m),"_ECDF_max_gradient_x_input_diff.pdf",sep=""),width=4,height=3)
    }
 }


## Plot example regions

genome <- BSgenome.Hsapiens.UCSC.hg38

example.regions <-
    lapply(focus.motifs, function(n) {
        tmp <- subset(motif.GI.score.df,motif==n & !is.na(max.GI.diff))
        list(tmp[which.max(tmp$max.GI.diff),c("region","pos")],
             tmp[which.min(tmp$max.GI.diff),c("region","pos")])
    })

regions <- unlist(example.regions,recursive=TRUE)[seq(1,length(example.regions)*4,by=2)]
ids <- match(unique(regions),regions)
regions <- regions[ids]
mregions <- unlist(example.regions,recursive=TRUE)[seq(2,length(example.regions)*4,by=2)][ids]

for (i in 1:length(regions)) {

    r <- regions[i]
    print(r)
    m <- mregions[i]
    gr <- classA.motif.GI.score[[r]]

    chr <- as.character(seqnames(classA.DHS.gr[as.character(seqnames(gr))[1]]))
    start <- start(classA.DHS.gr[as.character(seqnames(gr))[1]])
    end <- end(classA.DHS.gr[as.character(seqnames(gr))[1]])

    mstart <- as.numeric(strsplit(strsplit(m,":")[[1]][2],"-")[[1]][1])
    mend <- as.numeric(strsplit(strsplit(m,":")[[1]][2],"-")[[1]][2])

    seqlevels(gr) <- chr
    seqlengths(gr) <-  seqlengths(genome)[chr]
    genome(gr) <- "hg38"
    ranges(gr) <- IRanges(start=start(gr)+start-1,end=end(gr)+start-1)

    dgr <- GRanges(seqnames=chr,ranges = IRanges(start=start:(start+DHS.width-1),end=start:(start+DHS.width-1)))
    dgr$classA <- as.numeric(classA[r,1:ncol(classA)])
    dgr$classB <- as.numeric(classB[r,1:ncol(classB)])
    genome(dgr) <- "hg38"
    seqlengths(dgr) <-  seqlengths(genome)[chr]

    mtrack <- AnnotationTrack(gr,chromosome=chr,start=start(gr),genome="hg38",name="JASPAR 2022 motifs",)
    feature(mtrack) <- gr$motif
    dtrack <- DataTrack(dgr,chromosome=chr,start=start(dgr),end=end(dgr),genome="hg38",name="GI",type="h",groups=c(classA.name,classB.name))
    strack <- SequenceTrack(Hsapiens,chromosome=chr,genome="hg38",fontsize=8,fontface=1)
    atrack <- GenomeAxisTrack()
    tracks <- list(atrack,dtrack,strack,mtrack)

    pdf(paste(plot.dir,sub(":","_",r),".pdf",sep=""),width=5,height=2)
    plotTracks(list(atrack,dtrack),from=start,to=end,chromosome=chr)
    dev.off()

    pdf(paste(plot.dir,sub(":","_",r),"_motif.pdf",sep=""),width=5,height=5)
    try(plotTracks(tracks,from=mstart-10,to=mend+10,chromosome=chr,featureAnnotation = "feature",fontcolor.feature = 1,cex.feature = 0.3,type="b"))
    dev.off()
}
