## Getting Started

```sh
# Download sample PacBio from the PBcR website
wget -O- http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz | tar zxf -
ln -s selfSampleData/pacbio_filtered.fastq reads.fq
# Install minimap and miniasm (requiring gcc and zlib)
git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)
git clone https://github.com/lh3/miniasm  && (cd miniasm  && make)
# Overlap for PacBio reads (or use "-x ava-ont" for nanopore read overlapping)
minimap2/minimap2 -x ava-pb -t8 pb-reads.fq pb-reads.fq | gzip -1 > reads.paf.gz
# Layout
miniasm/miniasm -f reads.fq reads.paf.gz > reads.gfa
```

## Introduction

Miniasm is a very fast OLC-based *de novo* assembler for noisy long reads. It
takes all-vs-all read self-mappings (typically by [minimap][minimap]) as input
and outputs an assembly graph in the [GFA][gfa] format. Different from
mainstream assemblers, miniasm does not have a consensus step. It simply
concatenates pieces of read sequences to generate the final [unitig][unitig]
sequences. Thus the per-base error rate is similar to the raw input reads.

So far miniasm is in early development stage. It has only been tested on
a dozen of PacBio and Oxford Nanopore (ONT) bacterial data sets. Including the
mapping step, it takes about 3 minutes to assemble a bacterial genome. Under
the default setting, miniasm assembles 9 out of 12 PacBio datasets and 3 out of
4 ONT datasets into a single contig. The 12 PacBio data sets are [PacBio E.
coli sample][PB-151103], [ERS473430][ERS473430], [ERS544009][ERS544009],
[ERS554120][ERS554120], [ERS605484][ERS605484], [ERS617393][ERS617393],
[ERS646601][ERS646601], [ERS659581][ERS659581], [ERS670327][ERS670327],
[ERS685285][ERS685285], [ERS743109][ERS743109] and a [deprecated PacBio E.
coli data set][PB-deprecated]. ONT data are acquired from the [Loman
Lab][loman-ont].

For a *C. elegans* [PacBio data set][ce] (only 40X are used, not the whole
dataset), miniasm finishes the assembly, including reads overlapping, in ~10
minutes with 16 CPUs. The total assembly size is 105Mb; the N50 is 1.94Mb. In
comparison, the [HGAP3][hgap] produces a 104Mb assembly with N50 1.61Mb. [This
dotter plot][ce-img] gives a global view of the miniasm assembly (on the X
axis) and the HGAP3 assembly (on Y). They are broadly comparable. Of course,
the HGAP3 consensus sequences are much more accurate. In addition, on the whole
data set (assembled in ~30 min), the miniasm N50 is reduced to 1.79Mb. Miniasm
still needs improvements.

Miniasm confirms that at least for high-coverage bacterial genomes, it is
possible to generate long contigs from raw PacBio or ONT reads without error
correction. It also shows that [minimap][minimap] can be used as a read
overlapper, even though it is probably not as sensitive as the more
sophisticated overlapers such as [MHAP][mhap] and [DALIGNER][daligner].
Coupled with long-read error correctors and consensus tools, miniasm
may also be useful to produce high-quality assemblies.

## Algorithm Overview

1. Crude read selection. For each read, find the longest contiguous region
   covered by three good mappings. Get an approximate estimate of read
   coverage.

2. Fine read selection. Use the coverage information to find the good regions
   again but with more stringent thresholds. Discard contained reads.

3. Generate a [string graph][sg]. Prune tips, drop weak overlaps and collapse
   short bubbles. These procedures are similar to those implemented in
   short-read assemblers.

4. Merge unambiguous overlaps to produce unitig sequences.

## Limitations

1. Consensus base quality is similar to input reads (may be fixed with a
   consensus tool).

2. Only tested on a dozen of high-coverage PacBio/ONT data sets (more testing
   needed).

3. Prone to collapse repeats or segmental duplications longer than input reads
   (hard to fix without error correction).



[unitig]: http://wgs-assembler.sourceforge.net/wiki/index.php/Celera_Assembler_Terminology
[minimap]: https://github.com/lh3/minimap
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[gfa]: https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md
[ERS473430]: http://www.ebi.ac.uk/ena/data/view/ERS473430
[ERS544009]: http://www.ebi.ac.uk/ena/data/view/ERS544009
[ERS554120]: http://www.ebi.ac.uk/ena/data/view/ERS554120
[ERS605484]: http://www.ebi.ac.uk/ena/data/view/ERS605484
[ERS617393]: http://www.ebi.ac.uk/ena/data/view/ERS617393
[ERS646601]: http://www.ebi.ac.uk/ena/data/view/ERS646601
[ERS659581]: http://www.ebi.ac.uk/ena/data/view/ERS659581
[ERS670327]: http://www.ebi.ac.uk/ena/data/view/ERS670327
[ERS685285]: http://www.ebi.ac.uk/ena/data/view/ERS685285
[ERS743109]: http://www.ebi.ac.uk/ena/data/view/ERS743109
[PB-151103]: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly
[PB-deprecated]: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-20kb-Size-Selected-Library-with-P6-C4/ce0533c1d2a957488594f0b29da61ffa3e4627e8
[ce]: https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
[mhap]: https://github.com/marbl/MHAP
[daligner]: https://github.com/thegenemyers/DALIGNER
[sg]: http://bioinformatics.oxfordjournals.org/content/21/suppl_2/ii79.abstract
[loman-ont]: http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/
[hgap]: https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP
[ce-img]: http://lh3lh3.users.sourceforge.net/download/ce-miniasm.png
