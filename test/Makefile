
SNIPPA = ../bin/snippa.pl
REF = wildtype.fna
R1 = R1.fq
R2 = R2.fq
OUTDIR = out
BAM = $(OUTDIR)/reads.bam
GZIP = pigz

all: sim snippa

snippa:
	$(SNIPPA) makeref -r $(REF) -o $(OUTDIR)
	$(SNIPPA) align -r $(OUTDIR) -1 $(R1).gz -2 $(R2).gz -o $(OUTDIR)
	samtools flagstat $(BAM)
	$(SNIPPA) call -r $(OUTDIR) -b $(BAM) -o $(OUTDIR) -m bcftools
	
sim: clean
	wgsim wildtype.fna -1 150 -2 150 -N 5000 $(R1) $(R2) > snps.tab
	$(GZIP) -f $(R1)
	$(GZIP) -f $(R2)

clean:
	$(RM) -fr out
	$(RM) -f $(R1) $(R2) snps.tab
