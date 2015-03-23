#!/bin/bash

## Store arguments from bash command line in special array
args=("$@")
## use $# variable to print out
## number of arguments passed to the bash script
#echo Number of arguments passed: $# 
##echo arguments to the shell
#echo arg0 ${args[0]}
#echo arg1 ${args[1]}

export HMMBINDIR="/usr/local/hmmer-3.0/hmmer-3.0-linux-intel-x86_64/binaries/"
export PFAMDIR="/usr/local/hmmer-3.0/"
export PHYLIPDIR="/u1/NGS/Tools/phylip-3.69/exe/"

if [ $# -eq 3 ] ; then
    s=${args[0]}
    PFAM=${args[1]}
    BINDIR=${args[2]}    
else 
    if [ $# -eq 2 ] ; then
	s=${args[0]}
	PFAM=${args[1]}
	echo "You give only 2 argument! Using default value for BINDIR!"
    else 
	echo "You didn't give any arguments! Using default values!"
	s="GMC_oxred_C.matches_both.reduced"
	PFAM="PF05199"
    fi
    
    ##BINDIR="/home/momerja/proj/enox/bin/"
    
    ## for clarity the functions from enox were copied to:
    BINDIR="/mnt/msa1000-2/momerja/proj/pentoval/bin"
fi

echo family $s
echo family $PFAM
echo BINDIR $BINDIR


if [ -e "$s.fasta" ]; then
    echo "$s.fasta found!"

    if [ -e "$s.mod.fasta" ]; then
	echo "Gene names already shortened"
    else
	## Needed for the phylip step! Names restored in the end of the script
	$BINDIR/shortenAndSaveNames.perl -i $s.fasta -o $s.mod.fasta -n $s.shortnames.tab
    fi
    
    if [ -e "$s.$PFAM.stockholm" ]; then
	echo "HMMalign already run"
    else
	echo "Running HMMalign to $PFAM"
	export PFAMnum=`grep $PFAM /usr/local/hmmer-3.0/Pfam-A.hmm | perl -pe "s/ACC   //g"`
	$HMMBINDIR/hmmfetch  $PFAMDIR/Pfam-A.hmm $PFAMnum > $PFAMnum.hmm
	$HMMBINDIR/hmmalign $PFAMnum.hmm $s.mod.fasta > $s.$PFAM.stockholm
    fi
    
    if [ -e "$s.$PFAM.nex" ]; then
	echo "Alignment format already converted to nex"
    else
	## For creating the final output to be visualized in Geneious
	$BINDIR/convertAlignFormat_stockholm2nex.perl -i $s.$PFAM.stockholm
    fi
    
    if [ -e "$s.$PFAM.phy" ]; then
	echo "Alignment format already converted to phylip"
    else
	## For creating a tree using PHYLIP
	$BINDIR/convertAlignFormat_stockholm2phy.perl -i $s.$PFAM.stockholm
    fi
    
    if [ -e "$s.$PFAM.NJtree.nex" ]; then
	echo "Phylip already been run"
    else
	echo "Running PHYLIP"
	if [ -e infile ]; then
	    rm infile
	fi
	if [ -e outfile ]; then
	    rm outfile
	fi
	if [ -e outtree ]; then
	    rm outtree
	fi

	cp $s.$PFAM.phy infile
	$PHYLIPDIR/protdist > protdist.log <<EOF
Y
EOF
	cp outfile $s.$PFAM.protdist.outfile
	mv protdist.log $s.$PFAM.protdist.log
	mv outfile infile
	$PHYLIPDIR/neighbor > neighbor.log  <<EOF
Y
EOF
	cp outtree $s.$PFAM.neighbor.outtree
	mv neighbor.log $s.$PFAM.neighbor.log

	if [ -e "outtree" ]; then
	    echo "#NEXUS" > temp.nex
	    echo "" >> temp.nex
	    echo "begin trees;" >> temp.nex
	    echo "tree tnt_1 = [&U]" >> temp.nex
	    cat temp.nex outtree > $s.$PFAM.NJtree.nex
	    echo "" >> $s.$PFAM.NJtree.nex
	    echo "end;" >> $s.$PFAM.NJtree.nex

	    mv outfile $s.$PFAM.NJtree.visualization.txt
	fi
    fi
    
    ##
    ## File for geneious
    
    if [ -e "$s.$PFAM.alntree.nex" ]; then
	echo "$s.$PFAM.alntree.nex already exits"
    else 
	## Still doesn't work as it should. Unable to read this into Geneious, probably somethings 
        ## $s.$PFAM.nex and $s.$PFAM.NJtree.nex will work. Below these are also coverted back to full
	## length names

	    echo "Creating final alignment file"
	    cat $s.$PFAM.nex $s.$PFAM.NJtree.nex > temp.nex
	    $BINDIR/shortenAndSaveNames_Reverse.perl -i temp.nex -n $s.shortnames.tab -o $s.$PFAM.alntree.nex
	    rm temp.nex
    fi
    if [ -e "$s.$PFAM.NJtree.FullNames.nex" ]; then
	echo "$s.$PFAM.NJtree.FullNames.nex already exits"
    else 
	    echo "Restoring names"
	    $BINDIR/shortenAndSaveNames_Reverse.perl -i $s.$PFAM.NJtree.nex -n $s.shortnames.tab -o $s.$PFAM.NJtree.FullNames.nex
    fi
    if [ -e "$s.$PFAM.FullNames.nex" ]; then
	echo "$s.$PFAM.FullNames.nex already exits"
    else 
	    echo "Restoring names"
	    $BINDIR/shortenAndSaveNames_Reverse.perl -i $s.$PFAM.nex -n $s.shortnames.tab -o $s.$PFAM.FullNames.nex
    fi
else
    echo "Error! Coudn't find $s.fasta!"
fi

