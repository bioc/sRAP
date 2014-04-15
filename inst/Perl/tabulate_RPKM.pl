#Code written by Charles Warden (cwarden@coh.org, x60233)

use warnings;
use strict;

use Cwd 'abs_path'; 

$| =1;

my $os = $^O;
my $os_name;
if (($os eq "MacOS")||($os eq "darwin")||($os eq "linux"))
	{
		#Mac
		$os_name = "MAC";
	}#end if ($os eq "MacOS")
elsif ($os eq "MSWin32")
	{
		#PC
		$os_name = "PC";
	}#end if ($os eq "MacOS")
else
	{
		print "Need to specify folder structure for $os!\n";
		exit;
	}#end else
	
my $inputfile = $ARGV[0];
my $gene_file = $ARGV[1];
my $gene_index = $ARGV[2] - 1;
my $RPKM_index = $ARGV[3] - 1;

my %sample_hash=gene_fpkm_hash( $inputfile, 0, 1);

my %gene_hash;

		
		open(GENE, ">$gene_file") || die("Could not open $gene_file!");
		
		foreach my $sample (keys %sample_hash)
			{
				my $gene_file = $sample_hash{$sample};
				
				if(-f $gene_file)
					{
						print "Working on $sample....\n";
						my %temp_gene_hash=gene_fpkm_hash( $gene_file, $gene_index, $RPKM_index);
						
						if(scalar(keys %gene_hash) == 0)
							{
								%gene_hash = %temp_gene_hash;
								print GENE "Sample\t",join("\t",keys(%gene_hash)),"\n";;
							}#end if(scalar(keys %gene_hash) == 0)
							
						print GENE "$sample";
						
						foreach my $gene (keys %gene_hash)
							{
								my $FPKM = "NA";
								
								if(defined($temp_gene_hash{$gene}))
									{
										$FPKM=$temp_gene_hash{$gene};
									}
								
								print GENE "\t$FPKM";
							}
						print GENE "\n";
					}#end if(-f $gene_file)
			}#end foreach my $sample (keys %sample_hash)
		close(GENE);

		
exit;

sub gene_fpkm_hash
	{
		my ($inputfile, $gene_index, $RPKM_index)=@_;
		
		my %hash;
		
		my $line_count=0;
		open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
		while (<INPUTFILE>)
			{
				 $line_count++;
				 my $line = $_;
				 chomp $line;
				 if ($line_count > 1)
				 	{
				 		my @line_info = split("\t",$line);
						my $gene=$line_info[$gene_index];
						my $FPKM = $line_info[$RPKM_index];
						
						$hash{$gene}=$FPKM;
				 	}#end if ($line_count > 1)
			}#end while (<INPUTFILE>)
		close(INPUTFILE);
		
		return(%hash);
	}#end gene_fpkm_hash