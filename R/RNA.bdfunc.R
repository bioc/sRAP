ttest.pvalue <- function(arr)
	{
		pos.text <- arr[1]
		neg.text <- arr[2]
		pos.values <- as.numeric(unlist(strsplit(pos.text, ",")))
		neg.values <- as.numeric(unlist(strsplit(neg.text, ",")))
		
		if((length(pos.values) >= 3) & (length(neg.values) >= 3))
			{
				result <- t.test(pos.values, neg.values)
				return(result$p.value)
			}
		else
			{
				return(1)
			}

	}#end def ttest2 calc.pvalue

ttest.est <- function(arr)
	{
		pos.text <- arr[1]
		neg.text <- arr[2]
		pos.values <- as.numeric(unlist(strsplit(pos.text, ",")))
		neg.values <- as.numeric(unlist(strsplit(neg.text, ",")))

		if((length(pos.values) >= 3) & (length(neg.values) >= 3))
			{
				result <- t.test(pos.values, neg.values)
				return(result$statistic)
			}
		else
			{
				return(NA)
			}
	}#end def ttest2 calc.pvaluer
	
ks.est <- function(arr)
	{
		pos.text <- arr[1]
		neg.text <- arr[2]
		pos.values <- as.numeric(unlist(strsplit(pos.text, ",")))
		neg.values <- as.numeric(unlist(strsplit(neg.text, ",")))
		
		if((length(pos.values) >= 3) & (length(neg.values) >= 3))
			{
				result <- ks.test(pos.values, neg.values)
				return(result$statistic)
			}
		else
			{
				return(NA)
			}
	}#end def ttest2 calc.pvalue

mwu.est <- function(arr)
	{
		pos.text <- arr[1]
		neg.text <- arr[2]
		pos.values <- as.numeric(unlist(strsplit(pos.text, ",")))
		neg.values <- as.numeric(unlist(strsplit(neg.text, ",")))

		if((length(pos.values) >= 3) & (length(neg.values) >= 3))
			{
				result <- wilcox.test(pos.values, neg.values)
				return(result$statistic)
			}
		else
			{
				return(NA)
			}
	}#end def ttest2 calc.pvalue

mwu.pvalue <- function(arr)
	{
		pos.text <- arr[1]
		neg.text <- arr[2]
		pos.values <- as.numeric(unlist(strsplit(pos.text, ",")))
		neg.values <- as.numeric(unlist(strsplit(neg.text, ",")))

		if((length(pos.values) >= 3) & (length(neg.values) >= 3))
			{
				result <- wilcox.test(pos.values, neg.values)
				return(result$p.value)
			}
		else
			{
				return(1)
			}
	}#end def ttest2 calc.pvalue
	
ks.pvalue <- function(arr)
	{
		pos.text <- arr[1]
		neg.text <- arr[2]
		pos.values <- as.numeric(unlist(strsplit(pos.text, ",")))
		neg.values <- as.numeric(unlist(strsplit(neg.text, ",")))

		if((length(pos.values) >= 3) & (length(neg.values) >= 3))
			{
				result <- ks.test(pos.values, neg.values)
				return(result$p.value)
			}
		else
			{
				return(1)
			}
	}#end def ttest2 calc.pvalue

map.genes <- function(arr, gene.list, all.genes)
	{
		#print(length(arr))
		#print(all.genes)
		gene.names <- as.character(unlist(strsplit(as.character(gene.list), ",")))
		#print(gene.names)
		gene.signal <- arr[match(gene.names, all.genes, nomatch=0)]
		#print(gene.signal)
		gene.signal <- gene.signal[gene.signal != 0]
		gene.signal <- paste(gene.signal, collapse=",")
		return(gene.signal)
	}#end def ttest2 calc.pvalue
	
`RNA.bdfunc.fc` <-function (stat.table, project.name, project.folder, species=NULL, enrichment.file=NULL, p.method="t-test", p.adjust.method="fdr", plot.flag=TRUE)
{
	#stat table should not have duplicates
	bdfunc.folder<-file.path(project.folder,"BD-Func")
	dir.create(bdfunc.folder, showWarnings=FALSE)
	
	genes <- stat.table[[1]]
	log2.ratio <- NULL
	log2.ratio.name <- colnames(stat.table)[[2]] 
	if(length(grep(".log2.ratio",log2.ratio.name)) == 1)
		{
			log2.ratio <- stat.table[[2]]
		}
	else
		{
			stop("BD-Func requires a log2 ratio value for the primary variable.  This table doesn't look like a statistic table from RNA.deg()")
		}
		
	if(is.null(enrichment.file))
		{
			if(species == "human"){
				data(bdfunc.enrichment.human)
				enrichment.table <- bdfunc.enrichment.human
			} else if(species == "mouse"){
				data(bdfunc.enrichment.mouse)
				enrichment.table <- bdfunc.enrichment.mouse
			} else if(is.null(species)){
				stop("Need to specify a species!")
			} else {
				stop(paste(species," does not have an existing enrichment file.  Please specify one using the enrichment.file parameter",sep=""))
			}		
		}#end if(is.null(enrichment.file))
	else
		{
			enrichment.table <- read.table(enrichment.file, header=T, sep="\t")
		}
		
	pos.values <- array(dim=nrow(enrichment.table))
	neg.values <- array(dim=nrow(enrichment.table))
	
	if(plot.flag)
		{
			fig.folder<-file.path(bdfunc.folder,paste(project.name,"fc",species,sep="_"))
			dir.create(fig.folder, showWarnings=FALSE)
		}
	
	for (i in 1:nrow(enrichment.table))
		{
			pos.genes <- as.character(unlist(strsplit(as.character(enrichment.table[i,2]), ",")))
			neg.genes <- as.character(unlist(strsplit(as.character(enrichment.table[i,3]), ",")))
			
			pos.signal <- log2.ratio[match(pos.genes, genes,nomatch=0)]
			pos.signal <- pos.signal[pos.signal != 0]
			neg.signal <- log2.ratio[match(neg.genes, genes,nomatch=0)]
			neg.signal <- neg.signal[neg.signal != 0]
			
			if((plot.flag) & (length(pos.signal) > 3) & (length(neg.signal) > 3))
				{
					total.min <- min(pos.signal)
					total.max <- max(pos.signal)
					if(min(neg.signal) < total.min)
						{
							total.min <- min(neg.signal)
						}

					if(max(neg.signal) > total.max)
						{
							total.max <- max(neg.signal)
						}
					
					category.name <- enrichment.table[i,1]
					category.name <- gsub('\\s+','_',category.name,perl=TRUE)
					category.name <- gsub('\\\\','_',category.name)
					category.name <- gsub('/','_',category.name)
					density.file <- file.path(fig.folder, paste(category.name,"pdf",sep="."))
					pdf(file = density.file)
					den <- density(pos.signal, from =total.min, to=total.max)
					Signal <- den$x
					freq <- den$y
					plot(Signal, freq, type="l", lwd = 4, col="red", ylab = "Density", ylim=c(0,3))
					den <- density(neg.signal, from =total.min, to=total.max)
					Signal <- den$x
					freq <- den$y
					lines(Signal, freq, lwd = 4, col="green", type = "l")

					legend("topright", legend=c("Activated Genes","Inhibited Genes"), lwd = 4, col=c("red","green"))
					dev.off()
				}#end if(plot.flag)
			
			if((length(pos.signal) >= 3) & (length(neg.signal) >= 3))
				{
						pos.values[i] <- paste(pos.signal, collapse=",")
						neg.values[i] <- paste(neg.signal, collapse=",")
				}
			else
				{
						pos.values[i] <- NA
						neg.values[i] <- NA					
				}
		}#end for (i in 1:nrow(enrichment.table))
		
	data.table <- data.frame(pos=pos.values, neg=neg.values)
	filtered.categories <- enrichment.table[[1]]
	filtered.categories <- filtered.categories[!is.na(pos.values)]
	data.table <- data.table[!is.na(pos.values),]
	
	output.table <- data.frame(Category = filtered.categories)
	
	p.values <- array(dim=nrow(data.table))
	if(p.method == "t-test") {
		p.values <- apply(data.table, 1, ttest.pvalue)
	} else if(p.method == "mann-whitney") {
		p.values <- apply(data.table, 1, mwu.pvalue)
		mwu.statistic <- apply(data.table, 1, mwu.est)
	} else {
		p.values <- apply(data.table, 1, ks.pvalue)
		ks.statistic <- apply(data.table, 1, ks.est)
	}
			
	if(p.adjust.method == "q-value") {
		qobj <- qvalue(p.values)
		output.table <- data.frame(output.table, p.value = p.values, fdr = qobj$qvalues)
	} else if (p.adjust.method == "fdr") {
		fdr.values <- p.adjust(p.values, method="fdr")
		output.table <- data.frame(output.table,p.value = p.values, fdr = fdr.values)
	} else {
		output.table<- data.frame(output.table,p.value = p.values)
	}
	
	output.file <- file.path(bdfunc.folder, paste(project.name,"_fc_",species,".txt",sep=""))
	write.table(output.table, output.file, quote=F, row.names=F, sep="\t")
	
	print(warnings())
}#end def RNA.bdfunc.fc

`RNA.bdfunc.signal` <-function (expression.table, sample.file, project.name, project.folder, species=NULL, enrichment.file=NULL, p.method="t-test", p.adjust.method="fdr", plot.flag=TRUE, color.palette = c("red","blue","green","orange","purple","cyan","pink","maroon","yellow","grey","black",colors()))
{
	#signal.file should not have duplicates
	bdfunc.folder<-file.path(project.folder,"BD-Func")
	dir.create(bdfunc.folder, showWarnings=FALSE)
	
	sample.table <- read.table(sample.file, header=T, sep = "\t")
	samples <- as.character(sample.table[[1]])
	group.id <- sample.table[[2]]

	gene.names <- names(expression.table)[-1]
	sample.names <- expression.table[[1]]
	data.table <- expression.table[,2:ncol(expression.table)]
	
	if(length(samples) != length(sample.names[match(samples, sample.names, nomatch=0)]))
	{
		stop("Some samples in sample description file are not present in the expression file!")
	}

	if(length(samples)>2)
		{
			data.table <- data.table[match(samples, sample.names, nomatch=0),]
		}
	else
		{
			stop("You need at least 2 samples to conduct differential expression analysis!")
		}
	
	if(is.null(enrichment.file))
		{
			if(species == "human"){
				data(bdfunc.enrichment.human)
				enrichment.table <- bdfunc.enrichment.human
			} else if(species == "mouse"){
				data(bdfunc.enrichment.mouse)
				enrichment.table <- bdfunc.enrichment.mouse
			} else if(is.null(species)){
				stop("Need to specify a species!")
			} else {
				stop(paste(species," does not have an existing enrichment file.  Please specify one using the enrichment.file parameter",sep=""))
			}		
		}#end if(is.null(enrichment.file))
	else
		{
			enrichment.table <- read.table(enrichment.file, header=T, sep="\t")
		}
	
	output.table <- data.frame(Category = enrichment.table[[1]])
	p.values <- array(dim=nrow(enrichment.table))
	
	if(plot.flag)
		{
			fig.folder<-file.path(bdfunc.folder,paste(project.name,"signal",species,sep="_"))
			dir.create(fig.folder, showWarnings=FALSE)
		}
	
	for (i in 1:nrow(enrichment.table))
		{
			#print(dim(data.table))
			pos.mat <- apply(data.table, 1, map.genes, gene.list=enrichment.table[i,2], all.genes=gene.names)
			neg.mat <- apply(data.table, 1, map.genes, gene.list=enrichment.table[i,3], all.genes=gene.names)
			
			temp.data.frame <- data.frame(pos=pos.mat, neg=neg.mat)
			
			if(p.method == "t-test") {
				test.statistic <- apply(temp.data.frame, 1, ttest.est)
			} else if(p.method == "mann-whitney") {
				test.statistic <- apply(temp.data.frame, 1, mwu.est)
			} else {
				test.statistic <- apply(temp.data.frame, 1, ks.est)
			}
			
			box.data <- data.frame(group=group.id, score=as.numeric(test.statistic))
			box.factors <- as.factor(as.character(box.data$group[!is.na(box.data$score)]))
			
			if(length(levels(box.factors)) >= 2)
				{
					test <- aov(box.data$score~box.data$group) 
					result <- summary(test)[[1]][['Pr(>F)']][1]
					p.values[i] <- result
				}
			else
				{
					p.values[i] <- 1
				}
			
			if((plot.flag) & (length(levels(box.factors)) >= 2))
				{
					category.name <- enrichment.table[i,1]
					category.name <- gsub('\\s+','_',category.name,perl=TRUE)
					category.name <- gsub('\\\\','_',category.name)
					category.name <- gsub('/','_',category.name)
					
					if(typeof(group.id) == "double")
						{
							labelColors <- NULL
						}
					else
						{
							groups <- levels(group.id)
							if(p.method == "t-test") {
								labelColors <-array(dim=length(groups))
		
								for(i in 1:length(groups))
									{
										group.median <- median(box.data$score[box.data$group == groups[i]])
										if(group.median >= 2){
											labelColors[i] <- "red"
										} else if(group.median <= -2){
											labelColors[i] <- "green"
										}else{
											labelColors[i] <- "grey"
										}
									}#end for(i in 1:length(groups))
							} else {
								print(paste("Group: ",groups,sep=""))
								print(paste("Color: ",color.palette[1:length(groups)], sep=""))
								
								clusMember <- group.id
								labelColors <- as.character(clusMember)
								for (i in 1:length(groups))
									{
										#print(heatmap.label.colors[i])
										labelColors[clusMember == groups[i]] = color.palette[i]
									}
							}
						}
				
				box.file <- file.path(fig.folder, paste(category.name,"pdf",sep="."))
				pdf(file = box.file)
				test <- aov(box.data$score~box.data$group) 
				result <- summary(test)[[1]][['Pr(>F)']][1]
					
				boxplot(score~group, data=box.data, ylab="Test Statistic",col=labelColors, main = paste("ANOVA p-value = ",round(result, digits=2),sep=""))
				#stripchart(score~group, data=box.data, method='jitter', vertical=T, add=T, col=dot.colors, pch=16)
				dev.off()
				
				if((length(groups) == 2) & (tolower(groups)[1] == "negative") & (tolower(groups)[2] == "positive"))
					{
						rocr.labels <- as.character(group.id)
						rocr.labels[rocr.labels == "negative"] = 0
						rocr.labels[rocr.labels == "positive"] = 1
						
						ROCPred <- prediction(box.data$score, rocr.labels)
						ROCPerf<-performance(ROCPred,measure="tpr",x.measure="fpr")
						auc <- performance(ROCPred, "auc")
						auc <- unlist(slot(auc, "y.values"))
						roc.file <- file.path(fig.folder, paste(category.name,"_ROC.pdf",sep=""))
						pdf(file = roc.file)
						plot(ROCPerf,col="blue",xlim=c(0,1),ylim=c(0,1), main = paste("AUC = ",round(auc, digits=2),sep=""))
						dev.off()
					}#end if(roc.levels == groups)
				}#end if(plot.flag)
		}#for (i in 1:nrow(enrichment.table))

	if(p.adjust.method == "q-value") {
		qobj <- qvalue(p.values)
		output.table <- data.frame(output.table, p.value = p.values, fdr = qobj$qvalues)
	} else if (p.adjust.method == "fdr") {
		fdr.values <- p.adjust(p.values, method="fdr")
		output.table <- data.frame(output.table,p.value = p.values, fdr = fdr.values)
	} else {
		output.table<- data.frame(output.table,p.value = p.values)
	}
	
	output.file <- file.path(bdfunc.folder, paste(project.name,"_signal_",species,".txt",sep=""))
	write.table(output.table, output.file, quote=F, row.names=F, sep="\t")
		
	print(warnings())
}#end def RNA.bdfunc.signal