## Getting Started

```sh
# Install minimap and miniasm
git clone https://github.com/lh3/minimap && (cd minimap && make)
git clone https://github.com/lh3/miniasm && (cd miniasm && make)
# Overlapping
minimap/minimap -Sw5 -L100 -m0 -t8 reads.fa reads.fa | gzip -1 > reads.paf.gz
# Assembly
miniasm/miniasm -f reads.fa reads.paf.gz > reads.gfa
```

## Introduction

Miniasm is a very fast OLC-based *de novo* assembler for noisy long reads. It
takes all-vs-all read self-mappings (typically by [minimap][minimap]) as input
and outputs an assembly graph in the [GFA][gfa] format. Different from
mainstream assemblers, miniasm does not have a consensus step. It simply
concatenates pieces of read sequences to generate the final [unitig][unitig]
sequences. Thus the per-base error rate is similar to the raw input reads.

So far miniasm is in very early development stage. It has only been tested on
twelve bacterial genomes. Including the mapping step, it takes about 3 minutes
to assmble a bacterial genome. Under the default setting, miniasm assembles 5
out of 12 datasets into a single contig. The 12 data sets are [PacBio E. coli
sample][PB-151103], [ERS473430][ERS473430], [ERS544009][ERS544009],
[ERS554120][ERS554120], [ERS605484][ERS605484], [ERS617393][ERS617393],
[ERS646601][ERS646601], [ERS659581][ERS659581], [ERS670327][ERS670327],
[ERS685285][ERS685285], [ERS743109][ERS743109] and a [deprecated PacBio E. coli
sample][PB-deprecated].



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
[PB-deprecated]: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-20kb-Size-Selected-Library-with-P6-C4
