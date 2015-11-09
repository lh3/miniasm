prefix=ce-40X

# list of read files
cat > $prefix.files <<EOF
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0001/Analysis_Results/m140928_184123_42139_c100719602550000001823155305141590_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0001/Analysis_Results/m140928_184123_42139_c100719602550000001823155305141590_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0001/Analysis_Results/m140928_184123_42139_c100719602550000001823155305141590_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0002/Analysis_Results/m140928_230547_42139_c100719602550000001823155305141591_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0002/Analysis_Results/m140928_230547_42139_c100719602550000001823155305141591_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0002/Analysis_Results/m140928_230547_42139_c100719602550000001823155305141591_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0003/Analysis_Results/m140929_033247_42139_c100719602550000001823155305141592_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0003/Analysis_Results/m140929_033247_42139_c100719602550000001823155305141592_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0003/Analysis_Results/m140929_033247_42139_c100719602550000001823155305141592_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0004/Analysis_Results/m140929_075857_42139_c100719602550000001823155305141593_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0004/Analysis_Results/m140929_075857_42139_c100719602550000001823155305141593_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0004/Analysis_Results/m140929_075857_42139_c100719602550000001823155305141593_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0005/Analysis_Results/m140929_122654_42139_c100719602550000001823155305141594_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0005/Analysis_Results/m140929_122654_42139_c100719602550000001823155305141594_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590970/0005/Analysis_Results/m140929_122654_42139_c100719602550000001823155305141594_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0001/Analysis_Results/m140928_191128_sidney_c100699772550000001823139903261590_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0001/Analysis_Results/m140928_191128_sidney_c100699772550000001823139903261590_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0001/Analysis_Results/m140928_191128_sidney_c100699772550000001823139903261590_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0002/Analysis_Results/m140928_233028_sidney_c100699772550000001823139903261591_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0002/Analysis_Results/m140928_233028_sidney_c100699772550000001823139903261591_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0002/Analysis_Results/m140928_233028_sidney_c100699772550000001823139903261591_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0003/Analysis_Results/m140929_034941_sidney_c100699772550000001823139903261592_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0003/Analysis_Results/m140929_034941_sidney_c100699772550000001823139903261592_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0003/Analysis_Results/m140929_034941_sidney_c100699772550000001823139903261592_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0004/Analysis_Results/m140929_080908_sidney_c100699772550000001823139903261593_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0004/Analysis_Results/m140929_080908_sidney_c100699772550000001823139903261593_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0004/Analysis_Results/m140929_080908_sidney_c100699772550000001823139903261593_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0005/Analysis_Results/m140929_122826_sidney_c100699772550000001823139903261594_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0005/Analysis_Results/m140929_122826_sidney_c100699772550000001823139903261594_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0005/Analysis_Results/m140929_122826_sidney_c100699772550000001823139903261594_s1_p0.3.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0006/Analysis_Results/m140929_164720_sidney_c100699772550000001823139903261595_s1_p0.1.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0006/Analysis_Results/m140929_164720_sidney_c100699772550000001823139903261595_s1_p0.2.subreads.fasta
http://datasets.pacb.com.s3.amazonaws.com/2014/c_elegans/40X/raw_data/2590971/0006/Analysis_Results/m140929_164720_sidney_c100699772550000001823139903261595_s1_p0.3.subreads.fasta
EOF

# download read file
if [ ! -f $prefix.fa.gz ]; then
	wget -O- -qi $prefix.files | gzip -1 > $prefix.fa.gz
fi

# Install minimap and miniasm (requiring gcc and zlib)
git clone https://github.com/lh3/minimap && (cd minimap && make)
git clone https://github.com/lh3/miniasm && (cd miniasm && make)

# Overlap
minimap/minimap -Sw5 -L100 -m0 -t8 -I6G $prefix.fa.gz $prefix.fa.gz 2> $prefix.paf.gz.log | gzip -1 > $prefix.paf.gz

# Layout
miniasm/miniasm -f $prefix.fa.gz $prefix.paf.gz > $prefix.gfa 2> $prefix.gfa.log

# Convert to FASTA
awk '/^S/{print ">"$2"\n"$3}' $prefix.gfa > $prefix.utg.fa
