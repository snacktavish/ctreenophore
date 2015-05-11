##IDEA - mapping reads straight from SRA or other NGS files, to alignments and feeding them into the pipeline.

##Notes from Ido Ebersberger
## interrogate reads using hmmer or hmmerscan.
## Then assemble like they are tiny genomes using tirnity or something else!

## Shouldn't work for really distnt stuff:
## then Blast? reads to different *profiles*? and assemble.
## Need to think about what is appropriate at high divergence - translating reads to amino acids not scaleable.

##BUT : once potenially correct, eveni f quite messy, reads are targeted, THEN you can build sequences like they are genomes. Which they are.

#Option A
#Use aTRAM to generate seqs from reads.
for readarch in ${SRADIR}*.sra
  do
	a=$(basename $readarch); 
	stub=${a%.*}; 
	fastq-dump --split-files nam
	sed -e 's/length=\([0-9]*\)$/length=\1\\1/' ${stub}_1.fastq > ${stub}c_1.fastq 
	sed -e 's/length=\([0-9]*\)$/length=\1\\2/' ${stub}_2.fastq > ${stub}c_2.fastq 
	cat ${stub}c_1.fastq ${stub}c_2.fastq > ${stub}.fastq
  done

##Option 1 - pagan
##Is a logical option but is VERY slow
##I this case 'testNGSfasta' is just concateneated paired end reads
#pagan --ref-seqfile schirrmeister2011_1220_clean.aln -t RAxML_bestTree.schirrmeister2011_1220_clean --queryfile tesdtNGScombo.fastq --outfile read_alignment --fast-placement

##Option 2 - locus as reference
##FAST, but is locus is distant from reads none will map.
#bwa index -a bwtsw sequence.fasta 
#bwa mem sequence.fasta tesdtNGScombo.fastq > test.sam
#samtools view -bT sequence.fasta test.sam > test.bam
#samtools sort test.bam test_sorted
#samtools index test_sorted.bam test_sorted.bai
#samtools mpileup -uf sequence.fasta test_sorted.bam | bcftools view -bvcg - > var.raw.bcf
#samtools mpileup -uf sequence.fasta test_sorted.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq
#fastq_to_fasta -i cns.fq -o cns.fa

##Option 3
##Use ancestral seqs in full index?!
##Requires making some mildly annoying genome mapper variant file, so I haven't tried it.
#gmindex -i nogap.fas -x gmindextest -t tgmindextest
#genomemapper -i nogap.fas -x gmindextest -t tgmindextest -q tesdtNGScombo.fastq -f bed -o gmtest.bed


##Option 4.
##Use nhmmer to mathch reads, and then align only those using a slower mthod, i.e. Pagan
##Question: how to pull out matched reads? Grep seems silly...
#hmmbuild --dna -n schirrmhmm schirrmeister2011_1220_clean.aln 
#nhmmer --dna  --tblout hmmrhits schirrmhmm allreads.fas > hmm.out
#cut -f1 -d' ' hmmrhits | uniq | tail -n+3 > matches
#for item in `cat matches`; do grep -A 1 "$item " allreads.fas >> matched_reads; done
#pagan --ref-seqfile schirrmeister2011_1220_clean.aln -t RAxML_bestTree.schirrmeister2011_1220_clean --queryfile matched_reads --outfile read_alignment


##Option 5. 
##Might be best, but is a bit slow, and doesn't grab the distant stuff? Does output alignment tho
#hmmbuild --dna -n schirrmhmm schirrmeister2011_1220_clean.aln 
#hmmalign --dna --trim -o hmmalignment_reads.sth schirrmhmm matched_reads


##OPtion 6 using a short read aligner that can handle higher divergence
#./stampy.py -G hg18 ../ctreenophore/updating/sequence.fasta
#./stampy.py -g hg18 -H hg18
#./stampy.py -g hg18 -h hg18 -M ../ctreenophore/updating/matched_reads_1.fq  ../ctreenophore/updating/matched_reads_2.fq


#more on hmmalign:
#hmmsearch --tblout testsearch.out ../../ctreenophore/updating/schirrmhmm SRR610375.fasta > tbl.tes
#grep -e "     SRR6" tbl.test | tr -s ' ' | cut -d ' ' -f2 >matches
'''
small python prog:


matchset=set()
matches=open("matches")
for lin in matches:
	matchset.add(lin.strip())


matchout=open('matchout.fasta','w')

fi=open("../../fasta/SRR610375.fasta")
i=0
for lin in fi:
  i+=1
  if i%2==1:
  	lihead=lin
  else:
    if lihead.split()[0][1:] in matchset:
        matchout.write(lihead)
        matchout.write(lin)

matchout.close()
  


 2108  head test.fasta 
 2109  grep -e "     SRR6" tbl.test | tr -s ' ' | cut -d ' ' -f2 >matches
 2110  ls
 2111  head matches 
 2112  grep -e "     SRR6" tbl.test | tr -s ' ' | cut -d ' ' -f2 | uniq > matches
 2113  head test.fasta 
 2114  python
 2115  hmmalign --dna --mapali ../../ctreenophore/updating/schirrmeister2011_1220_clean.aln  ../../ctreenophore/updating/schirrmhmm matchout.fasta > searchout.aln
 2116  seaview searchout.aln 
 2117  sudo apt-get install minia
 2118  ./minia -in matchout.fasta -kmer-size 31 -abundance-min 3 -out output_minia
 2119  minia -in matchout.fasta -kmer-size 31 -abundance-min 3 -out output_minia
 2120  vi matchout.fasta 
 2121  minia -in test.fasta -kmer-size 31 -abundance-min 3 -out output_minia
 2122  minia -in test.fasta -kmer-size 3 -abundance-min 3 -out output_minia
 2123  minia iin test.fasta
 2124  minia test.fasta 31 1000 testmin
 2125  minia test.fasta 31 3 1000 testmin
 2126  cat matchout.fasta ../alttest/ref_genes.fa > goosed.fa
 2127  minia goosed.fasta 31 3 1000 testmin
 2128  minia matchout.fasta 31 3 1000 testmin
 2129  ls -ltr
 2130  vi testmin.contigs.fa 
 2131  perl ./stockholm2fasta.pl 
 2132  perl ./stockholm2fasta.pl testsearch.out 
 2133  ls -ltr
 2134  perl ./stockholm2fasta.pl searchout.aln > searchout.fasta
 2135  seaview searchout.fasta 
 2136  history
'''


bowtie2 -p "${PROCESSORS}" -N 1 --local -x "${MAINFOLDER}"/ref_genes -1 ${FILE} -2 $( echo ${FILE}|sed 's/R1/R2/' ) > >(tee ${FILE/.fastq/}_stdout.log) 2> >(tee ${FILE/.fastq/}_stderr.log >&2) | samtools view -Su -F 4 - | samtools sort - "$FOLDER"_loci/$(basename ${FILE/.fastq})  


bowtie2 -x full_aln/ref_genes -1 full_aln/datafiles/SRR610374_R1.fastq -2 full_aln/datafiles/SRR610374_R2.fastq -S tes_full.sam
samtools mpileup -uf REFERENCE.FA QUERY.BAM | bcftools view -cg - | vcfutils.pl vcf2fq > OUTPUT.FASTQ
