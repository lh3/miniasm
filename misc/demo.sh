#
# Dependencies: awk, wget, git, gcc and zlib
#

# Download Oxford Nanopore 2D reads (from the Loman lab)
wget -O reads.fa http://nanopore.climb-radosgw01.bham.ac.uk/MAP006-1_2D_pass.fasta

# Install minimap and miniasm (requiring gcc and zlib)
git clone https://github.com/lh3/minimap && (cd minimap && make)
git clone https://github.com/lh3/miniasm && (cd miniasm && make)

# Overlap
minimap/minimap -Sw5 -L100 -m0 -t8 reads.fa reads.fa | gzip -1 > reads.paf.gz

# Layout
miniasm/miniasm -f reads.fa reads.paf.gz > utg.gfa

# Convert to FASTA
awk '/^S/{print ">"$2"\n"$3}' utg.gfa > utg.fa

# Download E. coli K-12 sequence
wget -O NC_000913.fa 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&db=nuccore&dopt=fasta&val=556503834'

# Map assembly to ref
minimap/minimap NC_000913.fa utg.fa | miniasm/minidot - > utg.eps
