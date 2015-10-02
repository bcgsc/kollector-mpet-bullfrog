#!/usr/bin/make -Rrf
SHELL?=/bin/bash -o pipefail

#------------------------------------------------------------
# params
#------------------------------------------------------------

# threads
j?=1

# random seed
seed?=1

#------------------------------------------------------------
# meta rules
#------------------------------------------------------------

.DELETE_ON_ERROR:
.PHONY: check-params

default: check-params report-params $(name)-assembly.fa.gz \
	$(name).misassembled-ids.txt

check-params:
ifndef moleculo
	$(error missing required param `moleculo` (reference sequences))
endif
ifndef mpet
	$(error missing required param `mpet` (MPET FASTQ pair))
endif
ifndef name
	$(error missing required param `name` (output file prefix))
endif
ifndef pet
	$(error missing required param `pet` (PET FASTQ pairs))
endif
ifndef repeat_filter
	$(error missing required param `repeat_filter` (BioBloom filter containing repeat k-mers))
endif
ifndef seed
	$(error missing required param `seed` (random seed for sampling MPET pairs))
endif

#------------------------------------------------------------
# calculations
#------------------------------------------------------------

# genome size
# (5 Gbp genome * 0.3X Moleculo coverage = 1.5 Gbp)
genome_size?=1500000000

# mean MPET fragment size (5 kbp)
# NOTE: the second MPET fragment size peak at 8 kbp is not well
# captured by the Moleculo reads, because Moleculo N50 = 3.7 kbp.
mp_frag_size?=5000

# target coverage for assembled MPET fragments;
# used to calculate number of input MPET pairs
target_genome_cov?=0.01

# target read coverage for each MPET fragment
target_read_cov?=100

# approx sequencing error rate
error_rate?=0.005

# number of MPET fragments to assemble
num_mp:=$(shell perl -MPOSIX -e \
	'print floor($(target_genome_cov)*$(genome_size)/$(mp_frag_size)/3)')

# max kmers to recruit, based on formula:
#
#    n = G + Gec
#
# n: number of distinct kmers
# G: genome size (sum of MPET fragment sizes * 3)
# e: sequencing error rate
# c: read coverage
max_kmers:=$(shell perl -MPOSIX -e \
	'print ceil($(num_mp)*$(mp_frag_size)*3*(1 + $(error_rate)*$(target_read_cov)))')

# kollector-mpet options (see kollector-mpet --help)
kollector_opt:=-l30 -L30 -K25 -j$j -o $(name) -R $(name).report.txt \
	-m20 -r $(repeat_filter) -n$(max_kmers)

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

# report calculated parameters
report-params:
	@echo "number of MPET pairs =" $(num_mp)
	@echo "max k-mers to recruit =" $(max_kmers)

# sample input MPET reads (read 1)
$(name).mp_read1.fq.gz: $(mpet)
	seqtk sample -s $(seed) $(word 1, $(mpet)) $(num_mp) | gzip > $@

# sample input MPET reads (read 2)
$(name).mp_read2.fq.gz: $(mpet)
	seqtk sample -s $(seed) $(word 2, $(mpet)) $(num_mp) | gzip > $@

# do targeted assembly of MPET frags
$(name)-assembly.fa.gz: $(name).mp_read1.fq.gz $(name).mp_read2.fq.gz
	kollector-mpet $(kollector_opt) $^ $(pet)

# truncate assembled contigs to Moleculo lengths
# (for evaluation purposes)
$(name).truncated.fa.gz: $(name)-assembly.fa.gz
	chop-to-moleculo.mk \
		mp="$(name).mp_read1.fq.gz $(name).mp_read2.fq.gz" \
		moleculo=$(moleculo) contigs=$< name=$(name)

# evaluate truncated contigs
$(name).misassembled-ids.txt: $(name).truncated.fa.gz
	kollector-eval.mk ref=$(moleculo) assembly=$< name=$(name)
