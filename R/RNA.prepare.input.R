`RNA.prepare.input` <-function (sample.list, output.file, gene.index=1, rpkm.index=10)
{	
	Perl.Path <- file.path(path.package("sRAP"), "Perl")
	perl.script <- file.path(Perl.Path , "tabulate_RPKM.pl")
	cmd <- paste("perl \"",perl.script, "\" \"", sample.list,"\" \"", output.file,"\" ", gene.index," ", rpkm.index, sep="")
	system(cmd, intern=TRUE, wait=TRUE)
	
	warning(warnings())
}#end def RNA.deg