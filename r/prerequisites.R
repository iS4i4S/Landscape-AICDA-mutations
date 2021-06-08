# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load required packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


required.packages <- c('data.table','ggplot2','cowplot','RColorBrewer',
                       'ggsignif','binom','scales','forestplot','TCGAutils','ggpubr',"sqldf","seqinr",'annotate','TxDb.Hsapiens.UCSC.hg19.knownGene','stats','corrplot','reshape2',
                       'ggrepel','Hmisc','Rcpp','pheatmap','ComplexHeatmap','lawstat','e1071','survminer','matrixStats','IHW','DOSE','enrichplot','ggbeeswarm',
                       'rmarkdown','GSVA','survival','clusterProfiler','circlize',"tidyr")

lapply(required.packages, require, character.only = TRUE)
missing.packages <- required.packages[!required.packages %in% (.packages())]
if(length(missing.packages)>0) stop(paste('Could not load required packages:',paste(missing.packages,collapse=', ')))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define commonly used objects and helper functions (This functions were adapted from Gorelick et al [Composite mutations article]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

important_classes <- c('Missense_Mutation','Nonsense_Mutation','Splice_Site','In_Frame_Del','In_Frame_Ins',
                       'Frame_Shift_Del','Frame_Shift_Ins','TERT promoter','Translation_Start_Site')
truncating_classes <- c('Nonsense_Mutation','Splice_Site','Frame_Shift_Del','Frame_Shift_Ins')



## ggplot theme
theme_std <- function(base_size = 11, base_line_size = base_size/22, base_rect_size = base_size/22) {
    require(ggplot2)
    #theme_classic(base_size = base_size, base_family = 'ArialMT', base_line_size = base_line_size, base_rect_size = base_rect_size)  %+replace%
    theme_classic(base_size = base_size, base_family = 'ArialMT')  %+replace%
    theme(
          line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "round"),
          text = element_text(family = 'ArialMT', face = "plain",
                              colour = "black", size = base_size, lineheight = 0.9,
                              hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug=F),
          axis.text = element_text(colour = "black", family='ArialMT', size=rel(0.8)),
          axis.ticks = element_line(colour = "black", size=rel(1)),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = rel(1)),
          legend.key = element_blank(),
          strip.background = element_blank())
}


## function to force axis break in ggplot
break_axis <- function(y, maxlower, minupper=NA, lowerticksize, upperticksize, ratio_lower_to_upper) {
    if(is.na(minupper)) {
        breakpos <- maxlower
        lowerticklabels <- seq(0,breakpos,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(breakpos+upperticksize,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- breakpos + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > breakpos
        newy[ind] <- breakpos + uppertickspacing*((newy[ind]-breakpos) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    } else {
        lowerticklabels <- seq(0,maxlower,by=lowerticksize); lowerticklabels
        upperticklabels <- seq(minupper,max(y)+upperticksize,by=upperticksize); upperticklabels
        ticklabels <- c(lowerticklabels, upperticklabels); ticklabels
        lowertickpos <- lowerticklabels
        uppertickspacing <- ratio_lower_to_upper * lowerticksize
        uppertickpos <- maxlower + 0.5*lowerticksize + ((1:length(upperticklabels))*uppertickspacing)
        tickpos <- c(lowertickpos, uppertickpos)
        newy <- as.numeric(y)
        ind <- newy > maxlower
        newy[ind] <- maxlower + 0.5*lowerticksize + 1*uppertickspacing + uppertickspacing*((newy[ind]-minupper) / upperticksize)
        list(newy=newy, breaks=tickpos, labels=ticklabels, limits=range(tickpos))
    }
}


## function to convert ggplot object into separate objects for plot vs legend
extract_gglegend <- function(p){
    require(ggplot2)
    require(cowplot)

    ## extract the legend from a ggplot object
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) > 0) leg <- tmp$grobs[[leg]]
    else leg <- NULL
    leg

    ## return the legend as a ggplot object
    legend <- cowplot::ggdraw() + cowplot::draw_grob(grid::grobTree(leg))
    plot <- p + theme(legend.position='none')
    list(plot=plot,legend=legend)
}


## function to split HGVSp_Short field in MAF into Amino_Acid_Position, Reference_Amino_Acid, Variant_Amino_Acid, Amino_acids
HGVSp_Short_parse <- function(x,cpus=1) {
    require(parallel)

    x1 <- gsub('p.','',x)
    m <- gregexpr('[0-9]+',x1)
    nums <- regmatches(x1,m)
    get.aa <- function(num) num[1]
    aas <- sapply(nums,get.aa,USE.NAMES=F)
    tmp <- data.table(HGVSp_Short=x1,aa=aas)
    tmp$i <- 1:nrow(tmp)

    split <- function(d) {
        s <- strsplit(d$HGVSp_Short,d$aa)[[1]]
        s[2] <- gsub('_sice','splice',s[2])
        list(Reference_Amino_Acid=s[1],Variant_Amino_Acid=s[2])
    }
    info <- tmp[,split(.SD),by=i]
    info$HGVSp_Short <- x
    info$Amino_Acid_Position=as.integer(aas)
    info <- info[,c(4,2,5,3),with=F]
    info$Amino_acids <- paste(info$Reference_Amino_Acid,info$Variant_Amino_Acid,sep='/')
    info$Amino_acids[is.na(info$Amino_Acid_Position)] <- NA
    info
}


## creates a data.table with histogram of values in a vector 
table.freq <- function(value) {
    if(is.null(value) | length(value)==0) {
        tbl <- data.table(value=NA,N=NA)
    } else {
        tbl <- adt(table(value))
        tbl <- tbl[order(tbl$N,decreasing=T),]
    }
    tbl
}


## shortcut for sort(unique(...))
sortunique <- function(x,...) {
        sort(unique(x),na.last=T,...)
}


## shortcut to write tab-delimited data with consistent format 
write.tsv <- function(d, file, sep = "\t", quote = F, row.names = F, ...) {
    write.table(d, file = file, sep = sep, quote = quote, row.names = row.names, ...)
}


## shortcut for as.data.table
adt <- function(d) as.data.table(d)


## prompt the user to enter number of CPUs to use for parallel processing
## (for some scripts in r/ directory parallelization is advised)
cpu_prompt <- function() {
    cpus_available <- detectCores(all.tests = FALSE, logical = TRUE)
    prompt_text <- paste0('Wait! Parallization is recommended. You have ',cpus_available,' CPUs available. Enter number to use [1-',cpus_available,']: ')

    if (interactive() ) {
        input <- readline(prompt_text)
    } else {
        cat(prompt_text);
        input <- readLines("stdin",n=1);
    }

    cpus <- as.integer(input)
    if(!is.na(cpus) & cpus>0 & cpus<=cpus_available) {
        message('Proceeding with ',cpus,' CPUs ...')
    } else {
        stop('Invalid entry!')
    }   
    return(cpus)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define commonly used objects and helper functions defined here
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


signature_exposure_plot2<-function (mutSign_nums=NULL, mutSign_props, signature_colours = NULL, order=NULL,regression_line=NULL) 
{
    scale <- 1
    .theme_ss <- theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90, 
                                                                             vjust = 0.5, size = 8 * scale, family = "mono"), axis.text.y = element_text(hjust = 0.5, 
                                                                                                                                                         size = 12 * scale, family = "mono"), axis.text = element_text(size = 12 * 
                                                                                                                                                                                                                           scale, family = "mono"))
    if (is.null(order)=="TRUE") {
        ordering <- order(colSums(t(mutSign_nums)), decreasing = T)
    }else{
        ordering<-order
    }
    
    #mutSign_nums <- t(mutSign_nums)
    mutSign_props <- t(mutSign_props)
    #mutSign_nums <- mutSign_nums[, ordering]
    mutSign_props <- mutSign_props[, ordering]
    sample.ordering <- colnames(mutSign_props)
    #x1 <- melt(mutSign_nums)
    x2 <- melt(mutSign_props)
    #colnames(x1) <- c("Signature", "Sample", "Activity")
    colnames(x2) <- c("Signature", "Sample", "Activity")
    x2[, "class0"] <- c("Proportions")
    df2 <-  x2
    df2$class0 <- factor(df2$class0, c("Proportions"))
    df2$Sample <- factor(df2$Sample, sample.ordering)
    p = ggplot(df2, aes(x = factor(Sample), y = Activity, fill = Signature))
    p = p + geom_bar(stat = "identity", position = "stack")
    p = p + scale_fill_manual(values = signature_colours)
    p = p + ggtitle("Mutational Signature Exposures")
    if (length(intersect(df2$Signature, regression_line))>0) {
        p = p +  geom_smooth(data=df2[df2$Signature == regression_line,],aes(x = factor(Sample), y = Activity, group=1),color = "firebrick1", size = 0.5,show.legend = F,fill="#69b3a2") #Draw a regression line for the assigned signature
    }
    
    p = p + theme(plot.title = element_text(lineheight = 1, face = "bold", 
                                            size = 15 * scale))
    p = p + xlab("Samples") + ylab("Mutational Signature Content")
    p = p + theme(axis.title.x = element_text(face = "bold", 
                                              colour = "black", size = 15 * scale))
    p = p + theme(axis.title.y = element_text(face = "bold", 
                                              colour = "black", size = 15 * scale))
    p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                             size = 12 * scale, face = "bold", colour = "black"))
    p = p + theme(axis.text.y = element_text(size = 10 * scale, 
                                             face = "bold", colour = "black"))
    p = p + theme(legend.title = element_blank())
    p = p + .theme_ss
    p = p + theme(legend.position = "top")+theme( axis.ticks.x=element_blank())+scale_y_continuous(expand = c(0,0))
    return(p)
}


cor.test.p_2 <- function(data, extract=NULL, ...){
  p_vals<-list(); Tumor_type<-unique(data$type)
  cor.test.p <- function(x){
    FUN <- function(x, y) cor.test(x, y)[["p.value"]]
    z <- outer(
        colnames(x), 
        colnames(x), 
        Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}
  for(i in 1:length(Tumor_type)){
    x<-subset(data[,!(colnames(data) %in% "type")], data$type %in% Tumor_type[i])
     p_vals[[i]]<-cor.test.p(x)
  }
  if (extract=="NULL") {
     return(p_vals)
  }else{
  test_d<-matrix(data = 0, nrow =length(p_vals), ncol = ncol(p_vals[[1]]),dimnames = list(Tumor_type,colnames(p_vals[[1]])) ) #Create variable

for (i in 1:length(p_vals)) {
        test_d[i,]<- p_vals[[i]][,extract]
}
  #
  return(test_d)
  }
 
} #Function

cor.test.cor <- function(data, extract=NULL, ...){
  p_vals<-list(); Tumor_type<-unique(data$type)
  cor.test.p <- function(x){
    FUN <- function(x, y) cor.test(x, y)[["estimate"]][["cor"]]
    z <- outer(
        colnames(x), 
        colnames(x), 
        Vectorize(function(i,j) FUN(x[,i], x[,j]))
    )
    dimnames(z) <- list(colnames(x), colnames(x))
    z
}
  for(i in 1:length(Tumor_type)){
    x<-subset(data[,!(colnames(data) %in% "type")], data$type %in% Tumor_type[i])
     p_vals[[i]]<-cor.test.p(x)
  }
  if (extract=="NULL") {
     return(p_vals)
  }else{
  test_d<-matrix(data = 0, nrow =length(p_vals), ncol = ncol(p_vals[[1]]),dimnames = list(Tumor_type,colnames(p_vals[[1]])) ) #Create variable

for (i in 1:length(p_vals)) {
        test_d[i,]<- p_vals[[i]][,extract]
}
  #
  return(test_d)
  }
 
} #Function














