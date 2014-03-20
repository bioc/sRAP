\name{RNA.prepare.input}
\alias{RNA.prepare.input}

\title{
Prepare sRAP InputFile
}
\description{
Reads a table of samples containing RPKM (Read Per Kilobase per Million reads [1])
gene expression values and tabules results into a single table that can be read by 
RNA.norm().

Please note:
1) The default settings are designed for cufflinks [2] .fpkm_tracking files
2) The first file determines the set of geneIDs to be defined in the final table.
If this first file is missing RPKM/FPKM values for any genes, those genes will be
ignored in all subsequent samples.
}
\usage{
RNA.prepare.input(sample.list, output.file, gene.index=1, rpkm.index=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.list}{
	Table of samples RPKM expression values.
	Sample ID (used in sample description file) should be in the first column.
	Complete path to file should be in second column.
	Table should have column headers
}
  \item{output.file}{
	Table of RPKM/FPKM values for all samples.
	Genes are defined in columns, samples are defined in rows.
}
  \item{gene.index}{
	1-based index for gene ID in files described in sample.list.
	
	Default setting assumes use of cufflinks to define FPKM values.
}
  \item{rpkm.index}{
	1-based index for RPKM/FPKM values in files described in sample.list.
	
	Default setting assumes use of cufflinks to define FPKM values.
}
}
\value{
	Tab-delimited text file to be used for subsequent RNA.norm() step.
}
\references{

[1] Mortazavi A, Williams BA, McCue K, Schaeffer L, and Wold B. (2008).
Mapping and quantifying mammalian transcriptomes by RNA-Seq.
\emph{Nat Meth},
5:621-628.

[2] Trapnell C, Williams BA, Pertea G, Mortazavi A, Kwan G,
van Baren MJ, Salzberg SL, Wold BJ, Pachter L. (2010).
Transcript assembly and quantification by RNA-Seq reveals unannotated
transcripts and isoform switching during cell differentiation.
\emph{Nat Biotechnol},
28(5):511-5.
}
\author{
Charles Warden <cwarden45@gmail.com>
}

\examples{
	
library("sRAP")

dir <- system.file("extdata", package="sRAP")
cufflinks.folder <- file.path(dir,"cufflinks")
sample.table <- file.path(dir,"MiSeq_Sample_Description.txt")

samples <- c("SRR493372", "SRR493373","SRR493374","SRR493375","SRR493376","SRR493377")
cufflinks.files <- paste(samples,"_truncated.fpkm_tracking",sep="")
cufflinks.files <- file.path(cufflinks.folder, cufflinks.files)

project.folder <- getwd()
sample.mat <- data.frame(sample=samples, file=cufflinks.files)
sample.list <- file.path(project.folder, "cufflinks_files.txt")
write.table(sample.mat, file = sample.list, sep="\t", quote=FALSE, row.names=FALSE)
#You can view the "sample.list" file to see what it looks like
#For example, this sort of file can be created using Excel

rpkm.file <- file.path(project.folder, "sRAP_input.txt")
RNA.prepare.input(sample.list, rpkm.file)
}

\seealso{
Please post questions on the sRAP discussion group: http://sourceforge.net/p/bdfunc/discussion/srap/
}