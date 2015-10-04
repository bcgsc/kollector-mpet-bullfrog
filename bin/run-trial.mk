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
ifndef assembly_k
	$(error missing required param `assembly_k` (ABySS k-mer size))
endif
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
# (sum length of Moleculo reads = 1.6Gbp)
ref_size?=1600000000

# mean MPET fragment size = 3.2 kbp
# (MPET frags >= 1000bp that map pairwise to Moleculo reads)
mp_frag_size?=3200

# target coverage for assembled MPET fragments;
# used to calculate number of input MPET pairs
mpet_cov?=0.01

# target read coverage for each MPET fragment
kmer_cov?=350

# number of MPET fragments to assemble
num_mp?=$(shell perl -MPOSIX -e \
	'print floor($(mpet_cov)*$(ref_size)/$(mp_frag_size)/3)')

# max kmers to recruit
max_kmers?=$(shell perl -MPOSIX -e \
	'print ceil($(num_mp)*$(mp_frag_size)*3*$(kmer_cov))')

# kollector-mpet options (see kollector-mpet --help)
kollector_opt:=-l30 -L30 -K25 -j$j -o $(name) -R $(name).report.txt \
	-m20 -k$(assembly_k) -r $(repeat_filter) -n$(max_kmers)

#------------------------------------------------------------
# main rules
#------------------------------------------------------------

# report calculated parameters
report-params:
	@echo "Starting targeted assembly:"
	@echo "  Number of MPET pairs =" $(num_mp)
	@echo "  Max k-mers to recruit =" $(max_kmers)

# sample input MPET reads (read 1)
$(name).mp_read1.fq.gz: $(mpet)
	seqtk sample -s $(seed) $(word 1, $(mpet)) $(num_mp) | gzip > $@

# sample input MPET reads (read 2)
$(name).mp_read2.fq.gz: $(mpet)
	seqtk sample -s $(seed) $(word 2, $(mpet)) $(num_mp) | gzip > $@

# do targeted assembly of MPET frags
$(name)-assembly.fa.gz: $(name).mp_read1.fq.gz $(name).mp_read2.fq.gz
	kollector-mpet $(kollector_opt) $(KOLLECTOR_OPT) $^ $(pet)

# truncate assembled contigs to Moleculo lengths
# (for evaluation purposes)
$(name).truncated.fa.gz: $(name).fa.gz
	chop-to-moleculo.mk \
		mp="$(name).mp_read1.fq.gz $(name).mp_read2.fq.gz" \
		moleculo=$(moleculo) contigs=$< name=$(name)
	chop-to-moleculo.mk name=$(name) clean

# evaluate truncated contigs
$(name).misassembled-ids.txt: $(name).truncated.fa.gz
	kollector-eval.mk ref=$(moleculo) assembly=$< name=$(name)
