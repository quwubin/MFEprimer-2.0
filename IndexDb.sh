#!/bin/bash

# Index database for MFEprimer-2.0
# Wubin Qu <quwubin@gmail.com>

MFEHOME=$(dirname $(readlink -f "$0"))

if [ $# == 2 ]
then
    fasta_file=$1
    k=$2
elif [ $# == 1 ]
then
    fasta_file=$1
    k=9
else
    echo Usage:  
    echo
    echo     $(basename $0) Fasta_file K_value
    echo
    echo Example:  
    echo
    echo     $(basename $0) Human.fasta 9
    echo
    exit
fi

echo "Begin indexing ..."

$MFEHOME/chilli/UniFastaFormat.py -i $fasta_file

echo "Step 1/3: UniFasta done."

if [ `getconf LONG_BIT` == 64 ]
then
    $MFEHOME/bin/64bit/faToTwoBit $fasta_file.unifasta $fasta_file.2bit
else
    $MFEHOME/bin/32bit/faToTwoBit $fasta_file.unifasta $fasta_file.2bit
fi


echo "Step 2/3: faToTwoBit done."


echo "Step 3/3: Index begin ..."

$MFEHOME/chilli/mfe_index_db.py -f $fasta_file.unifasta -k $k
chmod 777 $fasta_file.uni

echo "Step 3/3: Index done"

rm $fasta_file.unifasta
