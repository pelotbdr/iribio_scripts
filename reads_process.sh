#!/bin/bash

if [ "$#" -ne 2 ]; then
	echo "Arguments incorrects, need two elements: folder name(1st) and reference file name(2nd)"
	exit 1
fi

cd $1
if ls ./F* >/dev/null  2>&1; then
	echo "Unzipping files..."
	gunzip -q *
	echo "Regrouping fastq files..."
	cat F*.fastq > complet.fastq
	echo "Deleting redundant fastq files..."
	rm F*.fastq
else
	if ls ./complet.fastq > /dev/null 2>&1; then
		echo "Regrouped fastq file already produced"
	else
		echo "No files found"
		exit 1
	fi
fi

if ls ./complet.sam > /dev/null 2>&1; then
	echo "Mapping file already present, delete the complet.sam file if you wish to redo the mapping"
	exit 1
else
	echo "Mapping..."
	minimap2 -ax map-ont "../"$2 complet.fastq > complet.sam
fi

echo "Finished!"
