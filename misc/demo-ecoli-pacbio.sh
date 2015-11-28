#
# Dependencies: awk, wget, git, gcc and zlib
#

# Download sample PacBio from the PBcR website
wget -O- http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData.tar.gz | tar zxf -
ln -s selfSampleData/pacbio_filtered.fastq reads.fq

# Install minimap and miniasm (requiring gcc and zlib)
git clone https://github.com/lh3/minimap && (cd minimap && make)
git clone https://github.com/lh3/miniasm && (cd miniasm && make)

# Overlap
minimap/minimap -Sw5 -L100 -m0 -t8 reads.fq reads.fq | gzip -1 > reads.paf.gz

# Layout
miniasm/miniasm -f reads.fq reads.paf.gz > utg.gfa

# Convert to FASTA
awk '/^S/{print ">"$2"\n"$3}' utg.gfa > utg.fa

# Download E. coli K-12 sequence
wget -O NC_000913.fa 'http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?sendto=on&db=nuccore&dopt=fasta&val=556503834'

# Map assembly to ref
minimap/minimap NC_000913.fa utg.fa | miniasm/minidot - > utg.eps
