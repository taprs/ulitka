``ulitka``: the Unfancy Liftover Toolkit in Awk <img src="logo/logo.svg" align="right" width="25%">
==========

Does what it says: `ulitka` handles **all stages of genomic coordinates liftover** from alignments of anchors to making/manipulating chain files and modifying coordinates in VCF and BED files. No sophisticated indel handling, thus only usable for coordinate conversion and SNP liftover. Sometimes that's the only thing needed, right?

If no premade chain file is available, `ulitka` can make it. It takes knowing where orthologous loci (anchors) are found in two genomes. Then, chain files describing the coordinate conversion are made from FASTA alignments of these regions. 

It comes with convenient **chain file manipulation** commands. Specifically, chains can be "stretched out" into a tabixable one-chain-one-line `lchain` format and "reversed" by swapping reference and query genomes (effectively enabling tabix queries across both coordinate systems).

# Imaginary example

```bash
# Make alignment: put at least two anchor sequences into a FASTA file, new (query) genome first
# Coordinates and orientation are detected from standard faidx naming
ulitka faidx new.fa chrA:1-10 > aln.fa
ulitka faidx --reverse-complement old.fa chrB:11-20 >> aln.fa

# Align with your favorite tool, e.g. mafft aln.fa > aln_fin.fa

# Make the chain file from any number of alignments
ulitka fa2chain -1 old.fa -2 new.fa -f aln_fin.fa > aln.chain

# Lift over VCF, modify header to contain correct ##contig definitions
ulitka easy -f new.fa -c aln.chain -v file.vcf.gz | bgzip > file_liftover.vcf.gz 
```

# Installation

Prerequisites are having the following programs in `$PATH` environmental variable:

 - `gawk`
 - `tabix`
 - `bgzip`
 - GNU `parallel` (to enable parallel execution in `ulitka easy`)

`ulitka` itself is zero-install, being entirely written in GNU AWK. Just clone the repo and make `ulitka` executable:

```bash
git clone https://github.com/taprs/ulitka.git
chmod +x ulitka/ulitka
```

