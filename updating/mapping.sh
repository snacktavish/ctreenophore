##IDEA - mapping reads straight from SRA or other NGS files, to alignments and feeding them into the pipeline.


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
