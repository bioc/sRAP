colLab <- function(n, labelColors, clusMember) { 
   if(is.leaf(n)) { 
       a <- attributes(n) 
	   #print(a)
       # clusMember - vector of sample names (ordered to match label color.palette)
       # labelColors - a vector of color.palette for the above grouping 
       labCol <- labelColors[clusMember == a$label]
	   #print(labCol)
       attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol) 
   } 
   n 
}

`RNA.qc` <-function (sample.file, expression.table, project.name, project.folder, plot.legend=TRUE, color.palette = c("red","blue","green","orange","purple","cyan","pink","maroon","yellow","grey","black",colors()))
{
qc.folder<-file.path(project.folder,"QC")
dir.create(qc.folder, showWarnings=FALSE)

sample.table <- read.table(sample.file, header=T, sep = "\t")
samples <- as.character(sample.table[[1]])
sample.group <- sample.table[[2]]
#print(expression.file)

sample.names <- expression.table[[1]]
expression.values <- t(expression.table[,2:ncol(expression.table)])
expr.max <- ceiling(max(expression.values))
expr.min <- ceiling(min(expression.values))

#print(samples)
#print(sample.names)

if(length(samples) != length(sample.names[match(samples, sample.names, nomatch=0)]))
	{
		print("Some samples in sample description file are not present in the expression file!")
		print(paste(length(samples),"items in sample description file",sep=" "))
		print(paste(length(sample.names),"items in gene expression file",sep=" "))
		print(paste(length(sample.names[match(samples, sample.names, nomatch=0)]),"matching items in gene expression file",sep=" "))
		#print(sample.names[match(samples, sample.names, nomatch=0)])
		stop()
	}

if(length(samples)>1)
	{
		expression.values <- expression.values[,match(samples, sample.names, nomatch=0)]
		colnames(expression.values)=samples
	}
#print(dim(expression.values))
print(samples)
#print(sample.names)
#print(colnames(expression.values))
#print(expression.values[1,])
rm(expression.table)

print(typeof(sample.group))


if(typeof(sample.group) == "double")
	{
		color.palette <- c("black")
	}
else
	{
		groups <- levels(sample.group)
		print(paste("Group: ",groups,sep=""))
		color.palette <- color.palette[1:length(groups)]
		print(paste("Color: ",color.palette[1:length(groups)], sep=""))
	}

#setwd(output.folder)
#print("Working directory has been set")
#sample histogram
q0 <- array(dim=length(samples))
q25 <- array(dim=length(samples))
q50 <- array(dim=length(samples))
q75 <- array(dim=length(samples))
q100 <- array(dim=length(samples))

#hist.file <- paste(project.name,"_hist.pdf")
hist.file <- file.path(qc.folder,paste(project.name,"_hist.pdf",sep=""))
pdf(file = hist.file)
#print(samples)
#print(length(samples))
for (i in 1:length(samples))
	{
		#print(i)
		#print(samples[i])
		#print(expression.values[,i])
		data <- -1
		if(length(samples)>1)
			{
				data <- as.numeric(t(expression.values[,i]))
			}
		else
			{
				data <- as.numeric(expression.values)
			}
		#print(length(data))
		quant <- quantile(data, na.rm=T)
		q0[i] <- quant[1]
		q25[i] <- quant[2]
		q50[i] <- quant[3]
		q75[i] <- quant[4]
		q100[i] <- quant[5]
		
		col <- "black"
		if(typeof(sample.group) != "double")
			{
				expression.group <- sample.group[i]
				#print(expression.group)
				for (j in 1:length(groups))
					{
						if(expression.group == groups[j])
							{
							col = color.palette[j]
							}
					}
			}#end if(typeof(sample.group) != "double")
		
		if(i == 1)
			{
				den <- density(data, na.rm=T,from=expr.min, to=expr.max)
				expr <- den$x
				freq <- den$y
				plot(expr, freq, type="l", xlab = "Log2(RPKM)", ylab = "Density", xlim=c(expr.min,expr.max), col=col)
				if(plot.legend)
					{
						legend("topright",legend=groups,col=color.palette,  pch=19)
					}
			}#end if(i == 1)
		else
			{
				den <- density(data, na.rm=T,from=expr.min, to=expr.max)
				expr <- den$x
				freq <- den$y
				lines(expr, freq, type = "l", col=col)
			}#end else
	}#end for (i in 1:length(bed.indices))
dev.off()
#print(samples)
#print(q50)
hist.table <- data.frame(sample = samples, min=q0, bottom25=q25, median=q50, top25=q75, max=q100)
hist.text.file <- file.path(qc.folder,paste(project.name,"_descriptive_statistics.txt",sep=""))
write.table(hist.table, hist.text.file, quote=F, row.names=F, sep="\t")
rm(hist.table)

if(length(samples) > 1)
	{
		#cluster tree
		#print(as.matrix(t(expression.values))[,1])
		dist1 <- dist(as.matrix(t(expression.values)))
		clusMember <- sample.group
		labelColors <- as.character(clusMember)
		for (i in 1:length(groups))
			{
			labelColors[clusMember == groups[i]] = color.palette[i]
			}
		hc <- hclust(dist1)
		rm(dist1)
		dend1 <- as.dendrogram(hc)
		rm(hc)
		cluster.file <- file.path(qc.folder,paste(project.name,"_cluster.pdf",sep=""))
		pdf(file = cluster.file)
		#print(labelColors)
		#print(clusMember)
		dend1 <- dendrapply(dend1, colLab, labelColors=labelColors, clusMember=samples) 
		a <- attributes(dend1) 
		attr(dend1, "nodePar") <- c(a$nodePar, lab.col = labelColors) 
		op <- par(mar = par("mar") + c(0,0,0,10)) 
		plot(dend1, horiz=T)
		par(op) 
		dev.off()
		rm(dend1)
		
		#box plot
		box.file <- file.path(qc.folder,paste(project.name,"_box_plot.pdf",sep=""))
		pdf(file=box.file)
		boxplot(expression.values, col=labelColors,xaxt='n')
		if(plot.legend)
			{
				legend("topright",legend=groups,col=color.palette,  pch=19)
			}
		dev.off()
		
		#PCA
		#print(dim(expression.values))
		#print(dim(na.omit(data.matrix(expression.values))))
		#print(data.matrix(expression.values)[1,])
		#print(na.omit(data.matrix(expression.values)[1,]))
		pca.values <- prcomp(na.omit(data.matrix(expression.values)))
		#print(attributes(pca.values))
		#print(pca.values$rotation)
		pc.values <- data.frame(pca.values$rotation)
		variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
		pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
		pca.text.file <- file.path(qc.folder,paste(project.name,"_pca.txt",sep=""))
		write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")
		pca.file <- file.path(qc.folder,paste(project.name,"_pca.pdf",sep=""))
		pdf(file=pca.file)
		plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19)
		if(plot.legend)
			{
				legend("topright",legend=groups,col=color.palette,  pch=19)
			}
		dev.off()
	}#end 
print(warnings())
}#end def RNA.qc