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

#------------- build AO term and name hash--------------
my %aoName;
my ($ao, $aoname);
open (AO, "ace_files/AOname.ace") || die "cannot open AOname.ace!\n";
while ($line = <AO>) {
    chomp($line);
    @tmp=(); 
    if  ($line =~ /^Anatomy_term/) {
	@tmp = split '"', $line;
	$ao = $tmp[1];
    } elsif ($line =~ /^Term/){
	@tmp = split '"', $line;
	$aoname = $tmp[1];
	$aoName{$ao} = $aoname;
	#print "$ao -- $aoname\n"; #for script testing
    }
}
close (AO);

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

#----------------- build EC anatomy hash ----
my %ecAO;
my @ecList;
my ($i, $ecAOdes);

$AceFile = join "", "ace_files/", $specode, "ECaols.ace";
open (EC, "$AceFile") || die "cannot open $AceFile!\n";
#open (EC, "ace_files/bmaECaols.ace") || die "cannot open bmaECmethod.ace!\n";

$line = <EC>;
while ($line = <EC>) {
    chomp($line);
    @tmp= (); 
    if  ($line =~ /^Expression_cluster/) {
	@tmp = split '"', $line;
	$ec = $tmp[1];
	$i = 0;
	@ecList = ();
    } elsif (($line =~ /^Anatomy_term/)&&($line =~ /Enriched/)) {
	#print "$line\n"; # for script testing
	@tmp = split '"', $line;
	$ao = $tmp[1];
	if ($aoName{$ao}) {
	    $ecList[$i] = $aoName{$ao};
	    $i++;
	} else {
	    print "ERROR: no term name for $ao!\n";
	}
    } elsif ($line eq "") {
	#next unless ($i > 0);
	$ecAOdes = join ",", @ecList;
	if ($ecAOdes ne "") {
	    $ecAO{$ec} = $ecAOdes;
	    #print "$ec -- $ecAOdes\n"; #for script testing
        }
    }
}
close (EC);

#----------- print genomic expression summary -----------

my $OutputFile1 = join "", $specode, "ECsummary_anatomy.csv";
open (OUT, ">$OutputFile1") || die "cannot open $OutputFile1!\n";
print OUT "GeneID\tPublic Name\tAnatomic Expression Summary\n";

$AceFile = join "", "ace_files/", $specode, "GeneWithEC.ace";
open (GENE, "$AceFile") || die "cannot open $AceFile!\n";
#open (GENE, "ace_files/test.ace") || die "cannot open $AceFile!\n";

my $gAOdes;
my @gECao;
my @tmpAOlist;
my @uniAOlist;
my %aoPrinted;

my $methodDes;
my @methodList;
#my @tmpMelist;
my @uniMElist;
my %mePrinted;


my ($j, $k, $m, $n);

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
	@gECao = ();
	@methodList = ();
    } elsif ($line =~ /^Expression_cluster/) {
	@tmp = split '"', $line;
	$ec = $tmp[1];

	if ($ecAO{$ec}) {
	    $gECao[$i] = $ecAO{$ec};
	    $methodList[$i] = $ecMethod{$ec};	    
	    #print "$ec -- $ecMethod{$ec}\n"; #for script testing
	    $i++;	
	}

    } elsif ($line eq "") {
	next unless ($i > 0);
	$gAOdes = join ",", @gECao;
	#print "$g: $gAOdes\n";  #for script testing
	@tmpAOlist = split ",", $gAOdes;

	$j = 0; #for AO list
	@uniAOlist = ();
	%aoPrinted = ();	
	foreach $ao (@tmpAOlist) {
	#foreach $ao (@gECao) {
	    if ($aoPrinted{$ao}) {
		#skip
	    } else {
		$aoPrinted{$ao} = 1;
		$uniAOlist[$j] = $ao;
		$j++;
	    }
	}

	if ($j == 0) {
	    #skip the gene
	} elsif ($j == 1) {
	    print OUT "$g\t$gName\tEnriched in $uniAOlist[0]";
	}  elsif ($j == 2) {
	    print OUT "$g\t$gName\tEnriched in $uniAOlist[0] and $uniAOlist[1]";
	}  elsif  ($j >= 3) {
	    print OUT "$g\t$gName\tEnriched in $uniAOlist[0]";
	    $k = 1;
	    while ($k <= ($j-2)) {
		print OUT ", $uniAOlist[$k]";
		$k++;
	    }
	    print OUT " and $uniAOlist[$k]";
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

