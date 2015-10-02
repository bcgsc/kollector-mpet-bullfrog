#!/usr/bin/make -Rrf
SHELL?=/bin/bash -o pipefail

#------------------------------------------------------------
# Description:
#
# This script aids in the evaluation of contigs from
# targeted assembly of MPET fragments.
#
# The assembled contigs typically extend beyond the ends
# of the reference Moleculo read.  As a result, aligning
# the contigs to the Moleculo reads results in clipping
# and chimeric alignments that look like misassemblies.
#
# In order to be able to detect *real* misassemblies,
# this script truncates the contigs to the lengths of their
# reference Moleculo reads.
#------------------------------------------------------------

#------------------------------------------------------------
# environment 
#------------------------------------------------------------

# make sure 'sort' and 'join' use the same sort order
LC_ALL=C

#------------------------------------------------------------
# params
#------------------------------------------------------------

# threads
j=1

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.DELETE_ON_ERROR:
.PHONY:
default: check-params $(name).truncated.fa.gz

check-params:
ifndef mp
	$(error missing required param 'mp' (mate pair reads))
endif
ifndef moleculo
	$(error missing required param 'moleculo' (moleculo reads))
endif
ifndef contigs
	$(error missing required param 'contigs' (assembled MPET fragments))
endif
ifndef name
	$(error missing required param 'name' (output file prefix))
endif

#------------------------------------------------------------
# alignment rules
#------------------------------------------------------------

# index FASTA file for alignment
%.bwt: %
	bwa index $<

# align MPET to moleculo
$(name).mpet-to-moleculo.sam.gz: $(moleculo).bwt $(mp)
	bwa mem -t$j $(moleculo) $(mp) | gzip > $@

# align MPET to contigs
$(name).mpet-to-contig.sam.gz: $(contigs).bwt $(mp)
	bwa mem -t$j $(contigs) $(mp) | gzip > $@

# get MPET "flank lengths" -- distances from left/right ends of MPET fragments
# to left/right ends of Moleculo reads
$(name).mpet-flanks.txt: $(name).mpet-to-moleculo.sam.gz
	zcat $< | sampe-flanks > $@

# get subseq coordinates to truncate each contig to the length of its reference
# Moleculo seq
$(name).truncated.bed: $(name).mpet-to-contig.sam.gz $(name).mpet-flanks.txt
	zcat $< | sampe-flanks -f $(name).mpet-flanks.txt > $@

# truncate each contig to length of its reference Moleculo seq
$(name).truncated.fa.gz: $(contigs) $(name).truncated.bed
	seqtk subseq $^ | gzip > $@
