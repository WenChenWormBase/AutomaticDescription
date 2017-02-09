#!/bin/csh
#-------------prepare expression cluster ace files ---------------------
setenv ACEDB /home/citace/WS/acedb/
## from Wen
/usr/local/bin/tace -tsuser 'wen' <<END_TACE
QUERY FIND Anatomy_term;
show -a -t Term -f ace_files/AOname.ace
QUERY FIND Life_stage;
show -a -t Public_name -f ace_files/LSname.ace

QUERY FIND Molecule;
show -a -t Public_name -f ace_files/MOLname.ace

QUERY FIND Expression_cluster; Species = *malayi
show -a -t Associated_with -f ace_files/bmaECaols.ace
show -a -t Regulation -f ace_files/bmaECreg.ace
show -a -t Attribute_of -f ace_files/bmaECmethod.ace
QUERY FIND Expression_cluster; Species = *malayi; follow Gene
show -a -t Expression_cluster -f ace_files/bmaGeneWithEC.ace

QUERY FIND Expression_cluster; Species = *pacificus
show -a -t Associated_with -f ace_files/ppaECaols.ace
show -a -t Regulation -f ace_files/ppaECreg.ace
show -a -t Attribute_of -f ace_files/ppaECmethod.ace
QUERY FIND Expression_cluster; Species = *pacificus; follow Gene
show -a -t Expression_cluster -f ace_files/ppaGeneWithEC.ace

QUERY FIND Expression_cluster; Species = *elegans
show -a -t Associated_with -f ace_files/ceECaols.ace
show -a -t Regulation -f ace_files/ceECreg.ace
show -a -t Attribute_of -f ace_files/ceECmethod.ace
QUERY FIND Expression_cluster; Species = *elegans; follow Gene
show -a -t Expression_cluster -f ace_files/ceGeneWithEC.ace

quit
END_TACE
