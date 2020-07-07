#!/usr/bin/bash
#SBATCH -p short

VERSION=46
PREFIX=https://fungidb.org/common/downloads/release-${VERSION}
FILE=CimmitisRS_Broad.mRNA.fasta
if [ ! -f data/$FILE ]; then
	curl -o data/$FILE $PREFIX/CimmitisRS/fasta/data/FungiDB-${VERSION}_CimmitisRS_AnnotatedTranscripts.fasta
fi
