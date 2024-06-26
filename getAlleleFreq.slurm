#!/bin/bash
#SBATCH --job-name=AlleleFreq
#SBATCH --mem=1GB
#SBATCH --time=5:00:00
#SBATCH --mail-user=rumika.mascarenhas@ucalgary.ca
#SBATCH --mail-type=END
#SBATCH --output=jobs.out
#SBATCH --error=jobs.err

GENOMIC_FILE="[SampleID]Blood_freebayes_filt_normalized_decompose.vcf"
AFFECTED_FILE="[SampleID]N_freebayes_filt_normalized_decompose.vcf"
UNAFFECTED_FILE="[SampleID]N_freebayes_filt_normalized_decompose.vcf"

GENOMIC_OUT="[SampleID]Blood.table"
AFFECTED_OUT="[SampleID]N.table"
UNAFFECTED_OUT="[SampleID]N.table"


#Comment this for samples with no blood/saliva
#/work/mtgraovac_lab/tools/gatk-4.1.8.1/gatk VariantsToTable \
#     -V /work/mtgraovac_lab/Rumika/SeqRound2/calling/freebayes/$GENOMIC_FILE \
#     -F CHROM -F POS -F REF -F ALT -F DP -GF AD \
#     -O $GENOMIC_OUT

/work/mtgraovac_lab/tools/gatk-4.1.8.1/gatk VariantsToTable \
      -V /work/mtgraovac_lab/Rumika/SeqRound2/calling/freebayes/$AFFECTED_FILE \
      -F CHROM -F POS -F REF -F ALT -F DP -GF AD \
      -O $AFFECTED_OUT

#Comment this for samples with no unaffected
/work/mtgraovac_lab/tools/gatk-4.1.8.1/gatk VariantsToTable \
   	-V /work/mtgraovac_lab/Rumika/SeqRound2/calling/freebayes/$UNAFFECTED_FILE \
    	-F CHROM -F POS -F REF -F ALT -F DP -GF AD \
    	-O $UNAFFECTED_OUT

#Format: python getAlleleFreq.py Genomic Affected Unaffected

python getAlleleFreq.py $GENOMIC_OUT $AFFECTED_OUT $UNAFFECTED_OUT

