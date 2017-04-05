#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import dendropy

dna1 = dendropy.DnaCharacterMatrix.get(file=open("test.aln"), schema="fasta")

print(dendropy.calculate.popgenstat.average_number_of_pairwise_differences(dna1, ignore_uncertain=True))
print(dendropy.calculate.popgenstat.num_segregating_sites(dna1, ignore_uncertain=True))
print(dendropy.calculate.popgenstat.tajimas_d(dna1, ignore_uncertain=True))


dna2 = dendropy.DnaCharacterMatrix.get(file=open("two_sample.aln"), schema="fasta")
print(dendropy.calculate.popgenstat._average_number_of_pairwise_differences(dna2.sequences(), dna2.default_state_alphabet, True))

