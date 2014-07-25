# HRC Likelihood Generation

This is a README for the SNP likelihood calculation for the [Haplotype Reference Consortium](http://haplotype-reference-consortium.org) (Release 1) using the [htslib](https://github.com/samtools/htslib) versions of [samtools](https://github.com/samtools/samtools)/[bcftools](https://github.com/samtools/bcftools).

1. Download and install this tagged (rc8) version of samtools/bcftools/htslib

        git clone --recursive git://github.com/mcshane/hrc-release1.git hrc-release1
        cd hrc-release1
        make

2. Run samtools mpileup on your cohort (or multi-cohort) BAM files to generate BCF2 files. These should be the BAM files associated with all the samples submitted to the project in the haplotype VCFs. The VCF `wgs.union.sites.breakMultiAC2.vcf.gz` (`md5sum : 78952ae50b90cb4707579cee7a3adaa5`) in the union site list currently only available within the project.

    The command can be split into multiple jobs for non-overlapping chunks across the genome. For a particular chunk, 

        samtools/samtools mpileup -IE -C50 -d$MAXDEPTH -t DP,DP4 -l wgs.union.sites.breakMultiAC2.vcf.gz -g -r $CHROM:$FROM-$TO -f $REF -b $BAMLIST > $CHROM:$FROM-$TO.$COHORT.bcf
        bcftools/bcftools index $CHROM:$FROM-$TO.$COHORT.bcf

    where

        $CHROM=chromosome
        $FROM=start position
        $END=end position
        $REF=GRCh37 reference fasta file used to align the BAM files, e.g. hs37d5.fa or human_g1k_v37.fasta
        $MAXDEPTH=Maximum depth per-BAM file. Should be adjusted depending on your coverage and if you have BAMs merged across multiple samples. Generally should be set to at least 10x the depth in your highest coverage input BAM.
        $BAMLIST=file listing location of input BAM files. One file per-line

    Positions are 1-based, so chunking `chr20` into 1Mbp chunks would look like `20:1-1000000`, `20:1000001-2000000`, ..., `20:63000001-63025520`. File `1Mbp-chunks.txt` could be used to generate the commands:

        awk '{print "samtools/samtools mpileup -IE -C50 -d$MAXDEPTH -t DP,DP4 -l wgs.union.sites.breakMultiAC2.vcf.gz -g -r $0 -f $REF -b $BAMLIST > $0.$COHORT.bcf"}' 1Mbp-chunks.txt

    If you want to generate BCF one sample at a time, then chunking is not necessary and the command would be:

        samtools/samtools mpileup -IE -C50 -d$MAXDEPTH -t DP,DP4 -l wgs.union.sites.breakMultiAC2.vcf.gz -g -f $REF $SAMPLE_BAM > $SAMPLE.bcf
        bcftools/bcftools index $SAMPLE.bcf

    If you generate per-sample BCF, skip step (3).

3. If you chunked in step (2), concat your chunks together with `bcftools concat`:

        bcftools/bcftools concat -Ob -f $CHROM.$COHORT.concat.list > $CHROM.$COHORT.bcf
        bcftools/bcftools index $CHROM.$COHORT.bcf

    where `$CHROM.$COHORT.concat.list` is a file listing the *coordinate-sorted* chunked BCF files (one file per-line) from step (2) for a particular chromosome and cohort. One can create this list with

        ls $CHROM:*.$COHORT.bcf | sort -V > $CHROM.$COHORT.concat.list

    If your list of files is not coordinate-sorted, then include option `-a` in the concat command.

4. Compare the sample names with the sample names in the haplotype VCFs that were submitted. The sample names in the BCF header can be queried with:

        bcftools/bcftools query -l $CHROM.$COHORT.bcf

    If the sample names are not the same, please fix here or provide a map of the sample IDs.

        bcftools/bcftools view -h $CHROM.$COHORT.bcf > old_header.vcf
        ** map sample names to make new header ** > new_header.vcf
        (cat new_header.vcf ; bcftools/bcftools view -H $CHROM.$COHORT.bcf) | bcftools/bcftools view -Ob > $CHROM.$COHORT.new_header.bcf

5. Submit your cohort likelihood BCFs.

6. _ONLY FOR CENTRAL ANALYSIS_. Once all likelihood BCFs for all cohorts are in, we will merge with the following:

        bcftools/bcftools norm -Ou -m+ wgs.union.sites.breakMultiAC2.vcf.gz | bcftools/bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' | htslib/bgzip -c > wgs.union.sites.AC2.alleles.gz
        htslib/tabix -p vcf wgs.union.sites.AC2.alleles.gz
        bcftools/bcftools merge -Ou -i DP:sum,I16:sum,QS:sum $CHROM.*.bcf | bcftools/bcftools call -Ou -mAC alleles -f GQ,GP -T wgs.union.sites.AC2.alleles.gz | bcftools/bcftools norm -Ob -m- > $CHROM.bcf
