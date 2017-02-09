#!/usr/bin/perl -w
use strict;

print "This script summarize expression cluster info from current WS release into automated description.\n";

if ($#ARGV !=0) {
    die "usage: $0 series_file ace\n";
}
my $specode = $ARGV[0];
my %speName = ("cbg" => "Caenorhabditis briggsae",
	       "cbn" => "Caenorhabditis brenneri",
	       "cja" => "Caenorhabditis japonica",
	       "cre" => "Caenorhabditis remanei",
	       "ce" => "Caenorhabditis elegans",
	       "ppa" => "Pristionchus pacificus",
	       "bma" => "Brugia malayi", 
	       "ovo" => "Onchocerca volvulus", 
	       "sra" => "Strongyloides ratti");

if ($speName{$specode}) {
    print "***** create genomic expression summary table for $speName{$specode} *****\n";    
} else {
    die "Species code $specode is not recognized.\n";
}

my @tmp;
my ($line, $tmp_length, $ec, $method, $AceFile);

#----------- build gene name hash ------------------
my ($g, $gName);
my %geneName;
open (ID, "/home/wen/fragmine/WBGeneName.csv") || die "cannot open !\n";
while ($line = <ID>) {
    chomp($line);
    @tmp = split /\t/, $line;
    $g = $tmp[0];
    $gName = $tmp[1];
    $geneName{$g} = $gName;
}
close (ID);

#------------------build  EC method hash  ---------------------

my %ecMethod;
my %methodName = ("Microarray_experiment" => "microarray", 
		  "Mass_spectrometry" => "proteomic study",
		  "RNASeq" => "RNA-seq",
		  "Tiling_array" => "tiling array",
		  "qPCR" => "quantitative PCR",
		  "Expr_pattern" => "Chronogram analysis");

$AceFile = join "", "ace_files/", $specode, "ECmethod.ace";
open (METHOD, "$AceFile") || die "cannot open $AceFile!\n";
while ($line = <METHOD>) {
    chomp($line);
    @tmp=(); 
    if  ($line =~ /^Expression_cluster/) {
	@tmp = split '"', $line;
	$ec = $tmp[1];
    } elsif ($line ne "") {
	@tmp = split /\s+/, $line;
	$method = $tmp[0];
	if ($methodName{$method}) {
	    $ecMethod{$ec} = $methodName{$method};
	} else {
	    $ecMethod{$ec} = "N.A.";
	    print "$ec -- no name for \"$method\"\n"; #for script testing
	}
	#if ($ec eq "WBPaper00046509:lin-35(n745)_downregulated") {
	#    print "$ec -- \"$method\" -- $methodName{$method}\n"; #for script testing
	#}
    }
}
close (METHOD);
$ec = "WBPaper00028802:intestine_unique";
$ecMethod{$ec} = "SAGE analysis";

#----------------- build EC GeneReg hash ----
my %ecGeneReg;
my @ecList;
my ($i, $ecGeneRegDes);

$AceFile = join "", "ace_files/", $specode, "ECreg.ace";
open (EC, "$AceFile") || die "cannot open $AceFile!\n";

$line = <EC>;
while ($line = <EC>) {
    chomp($line);
    @tmp= (); 
    if  ($line =~ /^Expression_cluster/) {
	@tmp = split '"', $line;
	$ec = $tmp[1];
	$i = 0;
	@ecList = ();
    } elsif ($line =~ /^Regulated_by_gene/) {
	#print "$line\n"; # for script testing
	@tmp = split '"', $line;
	$g = $tmp[1];
	if ($geneName{$g}) {
	    $ecList[$i] = $geneName{$g};
	    $i++;
	} else {
	    print "ERROR: no gene name for $g!\n";
	}
    } elsif ($line eq "") {
       
	$ecGeneRegDes = join ",", @ecList;
	if ($ecGeneRegDes ne "") {
	    $ecGeneReg{$ec} = $ecGeneRegDes;
	    #print "$ec -- $ecGeneRegDes\n"; #for script testing
        }
    }
}
close (EC);

#----------- print genomic expression summary -----------

my $OutputFile1 = join "", $specode, "ECsummary_geneReg.csv";
open (OUT, ">$OutputFile1") || die "cannot open $OutputFile1!\n";
print OUT "GeneID\tPublic Name\tGene Regulation Summary\n";

$AceFile = join "", "ace_files/", $specode, "GeneWithEC.ace";
open (GENE, "$AceFile") || die "cannot open $AceFile!\n";
#open (GENE, "ace_files/test.ace") || die "cannot open $AceFile!\n";

my $gGeneRegDes;
my @gECgr;
my @tmpGRlist;
my @uniGRlist;
my %grPrinted;

my $methodDes;
my @methodList;
#my @tmpMelist;
my @uniMElist;
my %mePrinted;


my ($j, $k, $m, $n, $gReg);

$line = <GENE>;
while ($line=<GENE>) {
    chomp ($line);
    #print "$line\n"; #for script testing
    @tmp= (); 

    if  ($line =~ /^Gene/) {
	@tmp = split '"', $line;
	$g = $tmp[1];
	if ($geneName{$g}) {
	    $gName = $geneName{$g};
	} else {
	    $gName = "N.A.";
	}
	
	$i = 0;
	@gECgr = ();
	@methodList = ();
    } elsif ($line =~ /^Expression_cluster/) {
	@tmp = split '"', $line;
	$ec = $tmp[1];

	if ($ecGeneReg{$ec}) {
	    $gECgr[$i] = $ecGeneReg{$ec};
	    $methodList[$i] = $ecMethod{$ec};	    
	    #print "$ec -- $ecMethod{$ec}\n"; #for script testing
	    $i++;	
	}

    } elsif ($line eq "") {
	next unless ($i > 0);
	$gGeneRegDes = join ",", @gECgr;
	#print "$g: $gGeneRegDes\n";  #for script testing
	@tmpGRlist = split ",", $gGeneRegDes;

	$j = 0; #for GeneReg Regulated_by_genelist
	@uniGRlist = ();
	%grPrinted = ();	
	foreach $gReg (@tmpGRlist) {
	    if ($grPrinted{$gReg}) {
		#skip
	    } else {
		$grPrinted{$gReg} = 1;
		$uniGRlist[$j] = $gReg;
		$j++;
	    }
	}

	if ($j == 0) {
	    #skip the gene
	} elsif ($j == 1) {
	    print OUT "$g\t$gName\tRegulated by $uniGRlist[0]";
	}  elsif ($j == 2) {
	    print OUT "$g\t$gName\tRegulated by $uniGRlist[0] and $uniGRlist[1]";
	}  elsif  ($j >= 3) {
	    print OUT "$g\t$gName\tRegulated by $uniGRlist[0]";
	    $k = 1;
	    while ($k <= ($j-2)) {
		print OUT ", $uniGRlist[$k]";
		$k++;
	    }
	    print OUT " and $uniGRlist[$k]";
	}

	$m = 0; #for method list
	@uniMElist = ();
	%mePrinted = ();
	foreach $method (@methodList) {
	    if ($mePrinted{$method}) {
		#skip
	    } else {
		$mePrinted{$method} = 1;
		$uniMElist[$m] = $method;
		$m++;
	    }
	}

	if ($m == 0) {
	    #skip the method
	} elsif ($m == 1) {
	    print OUT " according to $uniMElist[0].\n";
	}  elsif ($m == 2) {
	    print OUT " according to $uniMElist[0] and $uniMElist[1].\n";
	}  elsif  ($m >= 3) {
	    print OUT " according to $uniMElist[0]";
	    $n = 1;
	    while ($n <= ($m-2)) {
		print OUT ", $uniMElist[$n]";
		$n++;
	    }
	    print OUT " and $uniMElist[$n].\n";
	}

    }
}
close (GENE);
close (OUT);

