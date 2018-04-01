
SNIPPA = ../bin/snippa.pl
REF = wildtype.fna
DWGSIM = sim
R1 = $(DWGSIM).bwa.read1.fastq
R2 = $(DWGSIM).bwa.read2.fastq
OUTDIR = out
BAM = $(OUTDIR)/reads.bam
GZIP = pigz
VCF = snps.vcf
GFF = $(REF).gff
GBK = $(REF).gbk
#CALLER = bcftools
CALLER = freebayes
SEED = 42

all: sim snippa stats

snippa: $(GFF)
	$(SNIPPA) makeref -r $(REF) -o $(OUTDIR)
	$(SNIPPA) align -r $(OUTDIR) -1 $(R1) -2 $(R2) -o $(OUTDIR) 
	samtools flagstat $(BAM)
	$(SNIPPA) call -r $(OUTDIR) -b $(BAM) -o $(OUTDIR) -m $(CALLER)
	$(SNIPPA) snpeff -r $(OUTDIR) --vcf $(OUTDIR)/$(VCF) -o $(OUTDIR) --gff $(GFF)

$(GFF): $(REF)
	prokka --force --fast --outdir prokka --prefix anno $(REF) 
	mkdir -p $(OUTDIR)
	../bin/genbank-to-bcftools_gff.pl prokka/anno.gbk > $(GFF)
#	rm -fr prokka

stats:
	vt peek -r $(REF) $(DWGSIM).mutations.vcf
	vt peek -r $(REF) $(OUTDIR)/$(VCF)
#	@echo -n "Variants  : " && cat snps.tab | wc -l
#	@echo -n "SNPSs     : " && cut -f 3,4 snps.tab | grep -c -v -- '-'
#	@echo -n "Indels    : " && cut -f 3 snps.tab | grep -c -- '-'
#	@echo -n "Deletions : " && cut -f 4 snps.tab | grep -c -- '-'
	
sim: clean
	#wgsim wildtype.fna -S $(SEED) -h -s 100 -1 150 -2 150 -N 5000 $(R1) $(R2) > snps.tab
	dwgsim wildtype.fna -q B -z $(SEED) -C 20 -H -s 100 -1 150 -2 150 $(DWGSIM)
#	$(GZIP) -f $(R1)
#	$(GZIP) -f $(R2)

clean:
	$(RM) -fr out
	$(RM) -f $(DWGSIM).*
	
	
