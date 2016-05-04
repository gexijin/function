library(gplots)
hclust2 <- function(x, method="average", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))
  #####################################
# This function reads a gene list file in the GMT format. See GSEA web site for info.			 
readGMT <- function (file1) {
	# Read in the first file 
	x <- scan(file1, what="", sep="\n")
	x <- gsub("\t\t.","",x)     # GMT files saved by Excel has a lot of empty cells "\t\t\t\t"   "\t." means one or more tab
	x <- gsub(" ","",x)  # remove white space
	x <- toupper(x)    # convert to upper case

	#----Process the first file
	# Separate elements by one or more whitespace
	y <- strsplit(x, "\t")
	# Extract the first vector element and set it as the list element name
	names(y) <- sapply(y, `[[`, 1)
	#names(y) <- sapply(y, function(x) x[[1]]) # same as above
	# Remove the first vector element from each list element
	y <- lapply(y, `[`, -c(1,2))
	#y <- lapply(y, function(x) x[-1]) # same as above
	# remove duplicated elements
	for ( i in 1:length(y) )  y[[i]] <- unique(y[[i]])
	# check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
	if( max( sapply(y,length) ) <5) cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
	return(y)
}


myheatmap <- function (x,n=-1) {
if(n == -1) n=dim(x)[1]
geneSD = apply(x,1,sd)
x = x[order(-geneSD),]
# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
labRow = rownames(x)
if(dim(x)[1] > 300) labRow = F 
hy <-  heatmap.2(x, distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
#,Colv=FALSE,
,labRow = labRow
,key=F
,margins = c(6, 8)
)}
myheatmap3 <- function (x,bar = bar) {
if(length(bar) == 0 | length(bar) != dim(x)[1]) {
cat( "length does not match!"); stop() }
# this will cutoff very large values, which could skew the color 
x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff

nclass = length(unique(bar))
set.seed(2)
ncolors = sample( rainbow(nclass) )

heatmap.2(x,#distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
,dendrogram ="column"
,key=F, labRow = F
,Rowv = FALSE
,RowSideColors = ncolors[bar]
,margins = c(6, 8)
)
}
myheatmap2 <- function (x,bar = bar) {
if(length(bar) == 0 | length(bar) != dim(x)[1]) {
cat( "length does not match!"); stop() }
# this will cutoff very large values, which could skew the color 
x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
x = x[length(bar):1,]; bar = bar[length(bar):1] # reverse order
nclass = length(unique(bar))
set.seed(2)
ncolors = sample( rainbow(nclass) )
# use heatmap as it is faster
hy <-  heatmap(x,distfun = dist2,hclustfun=hclust2
, col=greenred(75),labRow = F,Rowv = NA
,RowSideColors = ncolors[bar],margins = c(6, 8))
}





###########################################################
# compute overlap 
################
# This program computes overlaps between two gene sets files in the GMT format. 
overlap <- function (file1, file2, total_elements = 35000, minFDR=0.05, minP=0.01 ) {
Pval_cutoff <- .01
#total_elements <- 30750 #24389 #31269   # This needs to be changed
#total_elements <- dim(genes)[1]
#total_elements <- length(universe)
#total_elements <- 35000
Min_overlap <- 1
Depeletion =FALSE    # Change to TRUE if want to identify depletion
minSetSize = 3; 
#file1 = "Clusters.gmt";
#file2 = "KEGG_Glycine_max_soybean_EnsemblID.gmt";
#file2 ="Soy_GO.txt"

# Read in the first file 
x <- scan(file1, what="", sep="\n")
#x <- gsub("\t\t.","",x)     # GMT files saved by Excel has a lot of empty cells "\t\t\t\t"   "\t." means one or more tab
#x <- gsub(" ","",x)  # remove white space
#x <- toupper(x)    # convert to upper case

#----Process the first file
# Separate elements by one or more whitespace
y <- strsplit(x, "\t")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -c(1,2))
#y <- lapply(y, function(x) x[-1]) # same as above
# remove duplicated elements
for ( i in 1:length(y) )  y[[i]] <- unique(y[[i]])
# check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
if( max( sapply(y,length) ) <5) cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
y = y[which(sapply(y,length) > minSetSize)]  # gene sets smaller than 10 is ignored!!!

#---Process second file
x2 <- scan(file2, what="", sep="\n")
#x2 <- gsub("\t\t.","",x2)  
#x2 <- gsub(" ","",x2)
#x2 <- toupper(x2)
# Separate elements by one or more whitepace
y2 <- strsplit(x2, "\t")
# Extract the first vector element and set it as the list element name
names(y2) <- sapply(y2, `[[`, 1)
y2 <- lapply(y2, `[`, -c(1,2))
# remove duplicated elements
for ( i in 1:length(y2) )  y2[[i]] <- unique(y2[[i]])
if( max( sapply(y2,length)) <5) cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
y2 = y2[which(sapply(y2,length) > minSetSize)]  # gene sets smaller than 10 is ignored!!!

#initialize a matrix to hold results
results <- matrix(ncol=8)
results_count <- matrix(ncol=length(y2), nrow=length(y))
Pval_table <- matrix(ncol=length(y2), nrow=length(y))
rownames(results_count) = names(y);   colnames(results_count) = names(y2)
rownames(Pval_table) = names(y);   colnames(Pval_table) = names(y2)
colnames(results) <- c("Cluster/Set1","n1","Set2","n2","#overlaped","Genes","Pval","FDR")
Pvalues = c()
Ntest=0;
#compute overlaps and P values
for ( i in 1:length(y) ) {
   cat(paste(i,"/",length(y),"\n"))
   if( length(y[[i]]) ==0    ) next; 
   for(j in 1:length(y2) ) {
       if( length(y2[[j]]) ==0    ) next; 
       ooo <- intersect(y[[i]], y2[[j]])
	   if (length(ooo) == 0) next; 
       results_count[i,j] = length(ooo)
       Ntest = Ntest+1;
       # if (length(ooo) <= Min_overlap) next  this may cause problem for detecting depletions
       xx <- length(ooo)
       mm <- length(y[[i]])
       nn <- total_elements - mm
       kk <- length(y2[[j]])
       if(nn<0 || kk> total_elements) { cat(paste(" Check total elements")); stop(); }
	   if (xx> 100) ooo <- ooo[1:100] # only keep the first 100 genes, if there is more
        # Note that it needs to be i-1 http://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c
       Pval_deplete=phyper(xx,mm,nn,kk, lower.tail=TRUE );
       Pval_enrich=phyper(xx-1,mm,nn,kk, lower.tail=FALSE );
       #Pval_table[i,j] = (-log10(Pval_enrich));   # if enrichment a positive number otherwise negative
       #if(Pval_deplete < Pval_enrich) Pval_table[i,j] = -(-log10(Pval_deplete))^(1/2);
	 if(Depeletion)   Pval=Pval_deplete  else Pval=Pval_enrich;
       if( Pval < Pval_cutoff ) {
           newO =c(names(y)[i],mm,names(y2)[j],kk,xx,paste(ooo,collapse=";"),Pval,"NA")
           results <- rbind(results,newO) 
   	     Pvalues = c(Pvalues,Pval)
        }
   }
}

results <- results[-1,]   # remove the first row as it is empty
if(dim(results)[1] <1) { cat("\nNo significant overlap found! Please "); stop() }
results <- as.data.frame(results)  #convert to data frame
results$FDR = (p.adjust(Pvalues,method="fdr",n=Ntest))
results <- results[ order( as.numeric(results[,1]),results[,8])  ,]  # sort according to FDR

results <- results[ which(results[,8]<minFDR)  ,]  # filter by FDR

return(results)
#results = results[,-6] # remove genes as there are too many


}




########################################################

myPGSEA  <- function (exprs, cl, range = c(25, 500), ref = NULL, center = TRUE, 
    p.value = 0.005, weighted = TRUE, nPermutation=100, enforceRange = TRUE, ...) 
{
    if (is(exprs, "ExpressionSet")) 
        exprs <- exprs(exprs)
    if (!is.list(cl)) 
        stop("cl need to be a list")
    if (!is.null(ref)) {
        if (!is.numeric(ref)) 
            stop("column index's required")
    }
    if (!is.null(ref)) {
        if (options()$verbose) 
            cat("Creating ratios...", "\n")
        ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
        exprs <- sweep(exprs, 1, ref_mean, "-")
    }
    if (center) 
        exprs <- scale(exprs, scale = FALSE)         # column centering is done
    results <- matrix(NA, length(cl), ncol(exprs))
    rownames(results) <- names(cl)
    colnames(results) <- colnames(exprs)
    mode(results) <- "numeric"
	Setsize = c(rep(0,length(cl)))     # gene set size vector
	mean2 = c(rep(0,length(cl)))     # mean of the range of means 
	meanSD = c(rep(0,length(cl)))     # SD of the range of means	
    if (is.logical(p.value)) 
        { p.results <- results; mean.results <- results;}
    for (i in 1:length(cl)) {              # for each gene list
		cat("\nProcessing gene set",i);
        if (class(cl[[i]]) == "smc") {
            clids <- cl[[i]]@ids
        }
        else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
            clids <- cl[[i]]@geneIds
        }
        else {
            clids <- cl[[i]]
        }
        if (options()$verbose) 
            cat("Testing region ", i, "\n")
        ix <- match(clids, rownames(exprs))
        ix <- unique(ix[!is.na(ix)])
        present <- sum(!is.na(ix))
		Setsize[i] <- present 
        if (present < range[1]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too small-", 
                  present, ",\n")
            next
        }
        if (present > range[2]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too large-", 
                  present, "\n")
            next
        }
        texprs <- exprs[ix, ]           # expression matrix for genes in gene set
        if (any(is.na(texprs))) 
            cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
        if (!is.matrix(texprs)) 
            texprs <- as.matrix(texprs)
                            
        stat <- try(apply(texprs, 2, t.test, ...))
		means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
		ps <- unlist(lapply(stat, function(x) x$p.value))
        stat <- unlist(lapply(stat, function(x) x$statistic))
        p.results[i, ] <- ps
		mean.results[i,] <- means
        results[i, ] <- as.numeric(stat)
		
		# permutation of gene sets of the same size
		if(nPermutation > 2 )  { # no permutation if <=2
			meansRanges = c(0, rep(nPermutation))
			for( k in 1:nPermutation ) {
				ix <- sample.int( dim(exprs)[1], length(ix) )
				texprs <- exprs[ix, ] 
				means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
				meansRanges[k] = dynamicRange(means)
			}
			mean2[i] = mean(meansRanges)
			meanSD[i]= sd(meansRanges,na.rm=TRUE)   # NA are removed before calculating standard deviation
		}
    }
    return(list(results = results, p.results = p.results, means = mean.results, size=Setsize, mean2=mean2, meanSD=meanSD))
    
}

dynamicRange <- function( x ) {
y = sort(x)
   if(length(x)>=4)  k =2 else k =1;
   return( y[length(x)-k+1] - y[k]) 
}



entropy <- function(x)
{  x = x/sum(x)
 x1 = log2(x+1e-50)
 return( sum(-x*x1) )
}
