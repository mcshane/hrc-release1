all: htslib/tabix htslib/bgzip samtools/samtools bcftools/bcftools

htslib/tabix:
	cd htslib && $(MAKE)

htslib/bgzip:
	cd htslib && $(MAKE)

samtools/samtools:
	cd samtools && $(MAKE)

bcftools/bcftools:
	cd bcftools && $(MAKE)
