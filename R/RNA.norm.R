`RNA.norm` <-function (input.file, project.name, project.folder, RPKM.cutoff = 0.1)
{
	data.folder<-file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	input.table <- read.table(input.file, header=T, sep="\t")
	sample.ids <- input.table[[1]]
	data.mat <- input.table[,2:ncol(input.table)]
	data.mat <- log(data.mat + RPKM.cutoff, base=2)
	
	norm.mat <- data.frame(ID=sample.ids, data.mat)
	
	if(nrow(norm.mat) < 1000000)
		{
			norm.expr <- t(data.mat)
			colnames(norm.expr) <- sample.ids
			gene.names <- rownames(norm.expr)
			norm.expr <- data.frame(Gene = gene.names, norm.expr)
			xlsfile <- file.path(data.folder, paste(project.name,"_normalized_expression.xlsx",sep=""))
			WriteXLS("norm.expr", ExcelFileName = xlsfile)
		}#end if(nrow(norm.mat) < 1000000)
	print(warnings())
	return(norm.mat)
}#end def RNA.norm