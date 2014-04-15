annova.pvalue <- function(y, x.mat)
{
	#print(arr)
	#print(grp.levels)
	
	y.no.na <- (y[!is.na(y)])
	if(length(y.no.na) >= 3)
		{
			test <- aov(y ~ ., data=data.frame(x.mat))
			#print(summary(test))
			result <- summary(test)[[1]][['Pr(>F)']]
			#print(length(result))
			result <- result[-length(result)]
			#print(result)
			if(is.null(result))
				{
					return(1)
				}
			else
				{
					return(result)
				}
		}
	else
		{
			return(1)
		}
}#end def annova.pvalue

lm.annova.pvalue <- function(y, x.mat)
{
	#print(arr)
	#print(grp.levels)
	
	y.no.na <- (y[!is.na(y)])
	if(length(y.no.na) >= 3)
		{
			reg <- lm(y ~ ., data=data.frame(x.mat)) 
			test <- anova(reg)
			result <- as.numeric(unlist(test['Pr(>F)']))
			#print(length(result))
			result <- result[-length(result)]
			#print(result)
			if(is.null(result))
				{
					return(1)
				}
			else
				{
					return(result)
				}
		}
	else
		{
			return(1)
		}
}#end def annova.pvalue

`RNA.deg` <-function (sample.file, expression.table, project.name, project.folder, log2.fc.cutoff=0.58, pvalue.cutoff=0.05, fdr.cutoff=0.05, box.plot = TRUE, ref.group=FALSE, ref="none",method = "lm", color.palette = c("red","blue","green","orange","purple","cyan","pink","maroon","yellow","grey","black",colors()))
{
	deg.folder<-file.path(project.folder,"DEG")
	dir.create(deg.folder, showWarnings=FALSE)
	
	data.folder<-file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	sample.table <- read.table(sample.file, header=T, sep = "\t")
	samples <- as.character(sample.table[[1]])
	sample.group <- sample.table[[2]]
	x.matrix <- sample.table[,2:ncol(sample.table)]
	x.names <- names(sample.table)[-1]

	gene.names <- names(expression.table)[-1]
	sample.names <- expression.table[[1]]
	expression.values <- t(expression.table[,2:ncol(expression.table)])
	expr.max <- ceiling(max(expression.values))
	expr.min <- ceiling(min(expression.values))
	
	if(length(samples) != length(sample.names[match(samples, sample.names, nomatch=0)]))
	{
		stop("Some samples in sample description file are not present in the expression file!")
	}

	if(length(samples)>2)
		{
			expression.values <- expression.values[,match(samples, sample.names, nomatch=0)]
			colnames(expression.values)=samples
		}
	else
		{
			stop("You need at least 2 samples to conduct differential expression analysis!")
		}
	
	stat.table <- data.frame(Gene=as.character(gene.names))
	
	if(expr.max > 20)
		{
			warning("Maximum expression value seems high.  Please make sure you are analysing data on a log2 scale.")
		}
	
	if(method == "lm"){
		stat.pvalues <- apply(expression.values, 1, lm.annova.pvalue, x.mat=x.matrix)
	} else if(method == "aov"){
		stat.pvalues <- apply(expression.values, 1, annova.pvalue, x.mat=x.matrix)
	} else{
	stop(paste(method," is not a valid method type!",sep=""))
	}

	#print(dim(stat.pvalues))
	#print(stat.pvalues[1,1])
	
	diff.index <- NULL
	fc.index <- NULL
	if(length(x.names) == 1)
		{
			if(ref.group)
				{
					group.ids <- x.matrix
					#print(group.ids)
					ref.mat <- expression.values[,group.ids == ref]
					#print(ref.mat[1:3,])
					trt.mat <- expression.values[,group.ids != ref]
					#print(trt.mat[1:3,])
					
					ref.avg.expr <- apply(ref.mat, 1, mean, na.rm = T)
					#print(ref.avg.expr[1:3])
					trt.avg.expr <- apply(trt.mat, 1, mean, na.rm = T)
					
					log2.ratio <- trt.avg.expr - ref.avg.expr
					temp.matrix <- data.frame(fc=log2.ratio)
					colnames(temp.matrix) <- c(paste(x.names[1],".log2.ratio",sep=""))
					rownames(temp.matrix) <- NULL
					stat.table <- data.frame(stat.table,temp.matrix)
				}#end if(ref.group)
			stat.fdr <- p.adjust(stat.pvalues, method="fdr")
			#print(dim(stat.fdr))
			#print(length(stat.fdr))
			temp.matrix <- data.frame(pvalue=stat.pvalues, fdr=stat.fdr)
			diff.index <- ((stat.pvalues < pvalue.cutoff) & (stat.fdr < fdr.cutoff))
			#print(diff.index)
			colnames(temp.matrix) <- c(paste(x.names[1],".pvalue",sep=""),paste(x.names[1],".fdr",sep=""))
			rownames(temp.matrix) <- NULL
			#print(dim(temp.matrix))
			stat.table <- data.frame(stat.table,temp.matrix)
		}#end if(length(x.names) == 1)
	else
		{
			if(ref.group)
				{
					group.ids <- x.matrix[1]
					#print(group.ids)
					ref.mat <- expression.values[,group.ids == ref]
					#print(ref.mat[1:3,])
					trt.mat <- expression.values[,group.ids != ref]
					#print(trt.mat[1:3,])
					
					ref.avg.expr <- apply(ref.mat, 1, mean, na.rm = T)
					#print(ref.avg.expr[1:3])
					trt.avg.expr <- apply(trt.mat, 1, mean, na.rm = T)
					
					log2.ratio <- trt.avg.expr - ref.avg.expr
					temp.matrix <- data.frame(fc=log2.ratio)
					colnames(temp.matrix) <- c(paste(x.names[1],".log2.ratio",sep=""))
					rownames(temp.matrix) <- NULL
					stat.table <- data.frame(stat.table,temp.matrix)
				}#end if(ref.group)
				
			for(i in 1:nrow(stat.pvalues))
				{
					temp.stat.pvalue <- as.array(stat.pvalues[i,])
					stat.fdr <- p.adjust(temp.stat.pvalue, method="fdr")
					#print(dim(stat.fdr))
					#print(length(stat.fdr))
					temp.matrix <- data.frame(pvalue=temp.stat.pvalue, fdr=stat.fdr)
					diff.index <- ((temp.stat.pvalue < pvalue.cutoff) & (stat.fdr < fdr.cutoff))
					#print(diff.index)
					colnames(temp.matrix) <- c(paste(x.names[i],".pvalue",sep=""),paste(x.names[i],".fdr",sep=""))
					rownames(temp.matrix) <- NULL
					print(dim(temp.matrix))
					#print(dim(temp.matrix))
					stat.table <- data.frame(stat.table,temp.matrix)
				}#end for(i in 1:ncol(stat.pvalues))
		}#end else
	
	#print(dim(stat.table))
	
	xlsfile <- file.path(data.folder, paste(project.name,"_all_data.xlsx",sep=""))
	WriteXLS("stat.table", ExcelFileName = xlsfile)
	
	diff.genes <- stat.table[diff.index,]

	if(ref.group)
		{
			fc.values <- diff.genes[[2]]

			up.genes <- na.omit(diff.genes[fc.values > log2.fc.cutoff,])
			up.genes <- up.genes[order(up.genes[[2]], decreasing = TRUE),]
			xlsfile <- file.path(deg.folder, paste(project.name,"_up.xlsx",sep=""))
			WriteXLS("up.genes", ExcelFileName = xlsfile)
			
			down.genes <- na.omit(diff.genes[fc.values < -log2.fc.cutoff,])
			down.genes <- down.genes[order(down.genes[[2]]),]
			xlsfile <- file.path(deg.folder, paste(project.name,"_down.xlsx",sep=""))
			WriteXLS("down.genes", ExcelFileName = xlsfile)
			
			any.change.genes <- na.omit(diff.genes[abs(fc.values) > log2.fc.cutoff,])
			print(dim(any.change.genes))
			any.change.genes <- any.change.genes[order(any.change.genes[[3]]),]
			print(dim(any.change.genes))
			xlsfile <- file.path(deg.folder, paste(project.name,"_DEG.xlsx",sep=""))
			WriteXLS("any.change.genes", ExcelFileName = xlsfile)
			
			diff.genes  <- any.change.genes
		}
	else
		{
			diff.genes <- na.omit(diff.genes)
			diff.genes<- diff.genes[order(diff.genes[[2]]),]
			xlsfile <- file.path(deg.folder, paste(project.name,"_DEG.xlsx",sep=""))
			WriteXLS("diff.genes", ExcelFileName = xlsfile)
		}
	
	diff.gene.names <- t(diff.genes[[1]])
	if(length(diff.gene.names) == 0)
		{
			print("There are no differentially expressed genes!")
			print(warnings())
			return(NULL)
		}
	else
		{
			diff.expr <- expression.values[match(diff.gene.names,gene.names,nomatch=0), ]
			print(dim(diff.expr))
			if(box.plot)
				{
					box.folder<-file.path(deg.folder,"Box_Plots")
					dir.create(box.folder, showWarnings=FALSE)
					
					for(i in 1:length(diff.gene.names))
						{
							gene.name <- diff.gene.names[i]
							gene.expr <- diff.expr[i,]
							
							box.file <- file.path(box.folder, paste(project.name,gene.name,"box_plot.pdf",sep="_"))
							pdf(file = box.file)
							boxplot(gene.expr ~ sample.group, col=color.palette)
							dev.off()
						}#end for(i in 1:length(diff.gene.names))
				}#end if(box.plot)
			
			diff.expr <- stdize(t(diff.expr))
			rownames(diff.expr) <- rep("",nrow(diff.expr))
			if(ncol(diff.expr) > 20)
				{
					colnames(diff.expr) <- rep("",ncol(diff.expr))
				}
			
			if(typeof(sample.group) == "double")
				{
					labelColors <- NULL
				}
			else
				{
					groups <- levels(sample.group)
					print(paste("Group: ",groups,sep=""))
					print(paste("Color: ",color.palette[1:length(groups)], sep=""))
					
					clusMember <- sample.group
					labelColors <- as.character(clusMember)
					for (i in 1:length(groups))
						{
							#print(heatmap.label.colors[i])
							labelColors[clusMember == groups[i]] = color.palette[i]
						}
				}
			#print(dim(diff.expr))
			heatmap.file <- file.path(deg.folder, paste(project.name,"_heatmap.pdf",sep=""))
			pdf(file = heatmap.file)
			heatmap.2(diff.expr, col=redgreen(33),density.info="none", trace="none", key=F, RowSideColors=labelColors)
			dev.off()
			
			print(warnings())
			return(stat.table)		
		}#end else
}#end def RNA.deg
