#!/bin/bash


##########################################
## Bash Script for Onepot-Seq alignment ##
##########################################


INHOUSE_SCRIPT_DIR=/Data2/KMS/Single_cell_analysis_tools/in_house_script
# PICARD_DIR=/Data2/KMS/Single_cell_analysis_tools/Drop-seq_tools-1.12/3rdParty/picard
PICARD_DIR=/Data2/SDJ/script/Drop-seq_tools-1.12/3rdParty/picard
STAR_DIR=/Data/programs/aligners/STAR-2.7.6a/bin/Linux_x86_64
# DROP_SEQ_TOOL_DIR=/Data2/KMS/Single_cell_analysis_tools/Drop-seq_tools-1.12
DROP_SEQ_TOOL_DIR=/Data2/SDJ/script/Drop-seq_tools-1.12



### Working Directory
WD=$1

### Script Directory
SCRIPT_DIR=$2

### Sample prefix
SAMPLE_PREFIX=$3

### Output number of cells
CELL=$4

### Reference identity 
### mix (human mouse mix) or hg (hg19)
Refident=$5

### Read1 fastq file
cDNA=$6

### Read2 fastq file
Barcode=$7



OUT=$SAMPLE_PREFIX
INHOUSE_SCRIPT_DIR=$SCRIPT_DIR/INHOUSE
PICARD_DIR=$SCRIPT_DIR/PICARD
STAR_DIR=$SCRIPT_DIR/STAR
DROP_SEQ_TOOL_DIR=$SCRIPT_DIR/DROPSEQ
REF_DIR=$SCRIPT_DIR/REF




mkdir $WD/pre-process

echo "Step 1 : Make fastq to BAM"
java -jar -Xmx4g $PICARD_DIR/picard.jar FastqToSam \
F1=$Barcode F2=$cDNA O=$WD/pre-process/"$OUT"_unaligned_data.bam SM=$OUT

echo "Step 2-1 : Tag BAM with Cell barcode"
$DROP_SEQ_TOOL_DIR/TagBamWithReadSequenceExtended \
INPUT=$WD/pre-process/"$OUT"_unaligned_data.bam \
OUTPUT=$WD/pre-process/"$OUT"_unaligned_tagged_Cell.bam \
SUMMARY=$WD/pre-process/"$OUT"_unaligned_tagged_Cellular.bam_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1

echo "Step 2-2 : Tag BAM with Molecular barcode"
$DROP_SEQ_TOOL_DIR/TagBamWithReadSequenceExtended \
INPUT=$WD/pre-process/"$OUT"_unaligned_tagged_Cell.bam \
OUTPUT=$WD/pre-process/"$OUT"_unaligned_tagged_CellMolecular.bam \
SUMMARY=$WD/pre-process/"$OUT"_unaligned_tagged_Molecular.bam_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=TRUE \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1

echo "Step 2-3 : Filter BAM"
$DROP_SEQ_TOOL_DIR/FilterBAM \
TAG_REJECT=XQ \
INPUT=$WD/pre-process/"$OUT"_unaligned_tagged_CellMolecular.bam \
OUTPUT=$WD/pre-process/"$OUT"_unaligned_tagged_filtered.bam

echo "Step 3-1 : Trim offtarget"
$DROP_SEQ_TOOL_DIR/TrimStartingSequence \
INPUT=$WD/pre-process/"$OUT"_unaligned_tagged_filtered.bam \
OUTPUT=$WD/pre-process/"$OUT"_unaligned_tagged_trimmed_smart.bam \
OUTPUT_SUMMARY=$WD/pre-process/"$OUT"_adapter_trimming_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5

echo "Step 3-2 : Trim polyA"
$DROP_SEQ_TOOL_DIR/PolyATrimmer \
INPUT=$WD/pre-process/"$OUT"_unaligned_tagged_trimmed_smart.bam \
OUTPUT=$WD/pre-process/"$OUT"_unaligned_mc_tagged_polyA_filtered.bam \
OUTPUT_SUMMARY=$WD/pre-process/"$OUT"_polyA_trimming_report.txt \
MISMATCHES=0 \
NUM_BASES=6

echo "Step 4 :  Make BAM to fastq"
mkdir $WD/star_align

java -Xmx4g -jar $PICARD_DIR/picard.jar SamToFastq \
INPUT=$WD/pre-process/"$OUT"_unaligned_mc_tagged_polyA_filtered.bam \
F=$WD/star_align/"$OUT"_unaligned_mc_tagged_polyA_filtered.fastq

if [ "$Refident" = "hg" ]
then
	echo "Step 6 : Align using STAR"
	$STAR_DIR/STAR \
	--genomeDir $REF_DIR/hg_ref \
	--readFilesIn $WD/star_align/"$OUT"_unaligned_mc_tagged_polyA_filtered.fastq
		
	mv Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab $WD/star_align

	java -Xmx4g -jar $PICARD_DIR/picard.jar SamFormatConverter \
	I=$WD/star_align/Aligned.out.sam \
	O=$WD/star_align/Aligned.out.bam

	echo "Step 7 : Sort SAM"

	java -Xmx4g -jar $PICARD_DIR/picard.jar SortSam \
	I=$WD/star_align/Aligned.out.bam \
	O=$WD/star_align/accepted_hits_sorted.bam \
	SO=queryname

	echo "Step 8 : Merge aligned data with barcode data"
	java -Xmx4g -jar $PICARD_DIR/picard.jar MergeBamAlignment \
	REFERENCE_SEQUENCE=$REF_DIR/hg_ref/GSM1629193_hg19_ERCC.fasta \
	UNMAPPED_BAM=$WD/pre-process/"$OUT"_unaligned_mc_tagged_polyA_filtered.bam \
	ALIGNED_BAM=$WD/star_align/accepted_hits_sorted.bam \
	OUTPUT=$WD/star_align/"$OUT"_merged.bam \
	INCLUDE_SECONDARY_ALIGNMENTS=false \
	PAIRED_RUN=false


	echo "Step 9 : Tag read with Exon"
	$DROP_SEQ_TOOL_DIR/TagReadWithGeneExon \
	I=$WD/star_align/"$OUT"_merged.bam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged.bam \
	ANNOTATIONS_FILE=$REF_DIR/hg_ref/GSM1629193_hg19_ERCC.gtf \
	TAG=GE


	echo "Step 10 : Detect bead synthesis errors"
	$DROP_SEQ_TOOL_DIR/DetectBeadSynthesisErrors \
	I=$WD/star_align/"$OUT"_star_gene_exon_tagged.bam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.bam \
	OUTPUT_STATS=$WD/star_align/"$OUT"_star_gene_exon_tagged.synthesis_stats.txt \
	SUMMARY=$WD/star_align/"$OUT"_star_gene_exon_tagged.synthesis_stats.summary.txt \
	NUM_BARCODES=$CELL \
	PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
	
	echo "Step 12 : Make Digital Gene Expression Matrix"
	mkdir $WD/dge
	java -jar -Xmx8g $DROP_SEQ_TOOL_DIR/jar/dropseq.jar \
	DigitalExpression I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.bam \
	O=$WD/dge/"$OUT"_hg_dge.txt \
	EDIT_DISTANCE=1 \
	NUM_CORE_BARCODES=$CELL \
	SUMMARY=$WD/dge/"$OUT"_hg_dge.summary.txt
	
	
	

elif [ "$Refident" = "mix" ]
then
	echo "Step 6 : Align using STAR"
	$STAR_DIR/STAR \
	--genomeDir $REF_DIR/hg_mm_ref \
	--readFilesIn $WD/star_align/"$OUT"_unaligned_mc_tagged_polyA_filtered.fastq
	
	mv Aligned.out.sam Log.final.out Log.out Log.progress.out SJ.out.tab $WD/star_align

	java -Xmx4g -jar $PICARD_DIR/picard.jar SamFormatConverter \
	I=$WD/star_align/Aligned.out.sam \
	O=$WD/star_align/Aligned.out.bam

	echo "Step 7 : Sort SAM"

	java -Xmx4g -jar $PICARD_DIR/picard.jar SortSam \
	I=$WD/star_align/Aligned.out.bam \
	O=$WD/star_align/accepted_hits_sorted.bam \
	SO=queryname

	echo "Step 8 : Merge aligned data with barcode data"
	java -Xmx4g -jar $PICARD_DIR/picard.jar MergeBamAlignment \
	REFERENCE_SEQUENCE=$REF_DIR/hg_mm_ref/hg19_mm10_transgenes.fa \
	UNMAPPED_BAM=$WD/pre-process/"$OUT"_unaligned_mc_tagged_polyA_filtered.bam \
	ALIGNED_BAM=$WD/star_align/accepted_hits_sorted.bam \
	OUTPUT=$WD/star_align/"$OUT"_merged.bam \
	INCLUDE_SECONDARY_ALIGNMENTS=false \
	PAIRED_RUN=false


	echo "Step 9 : Tag read with Exon"
	$DROP_SEQ_TOOL_DIR/TagReadWithGeneExon \
	I=$WD/star_align/"$OUT"_merged.bam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged.bam \
	ANNOTATIONS_FILE=$REF_DIR/hg_mm_ref/hg19_mm10_transgenes.gtf \
	TAG=GE

	echo "Step 10 : Detect bead synthesis errors"
	$DROP_SEQ_TOOL_DIR/DetectBeadSynthesisErrors \
	I=$WD/star_align/"$OUT"_star_gene_exon_tagged.bam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.bam \
	OUTPUT_STATS=$WD/star_align/"$OUT"_star_gene_exon_tagged.synthesis_stats.txt \
	SUMMARY=$WD/star_align/"$OUT"_star_gene_exon_tagged.synthesis_stats.summary.txt \
	NUM_BARCODES=$CELL \
	PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

	java -Xmx4g -jar $PICARD_DIR/picard.jar SamFormatConverter \
	I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.bam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.sam

	echo "Step 11 : Separate Human data and Mouse data"
	$INHOUSE_SCRIPT_DIR/Separate_sam.py \
	$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.sam \
	$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_hg.sam \
	$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_mm.sam \
	$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_separate_summary.txt

	java -Xmx4g -jar $PICARD_DIR/picard.jar SamFormatConverter \
	I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_hg.sam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_hg.bam

	java -Xmx4g -jar $PICARD_DIR/picard.jar SamFormatConverter \
	I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_mm.sam \
	O=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_mm.bam

	rm $WD/star_align/*.sam

	echo "Step 12 : Make Digital Gene Expression Matrix"
	mkdir $WD/dge
	java -jar -Xmx8g $DROP_SEQ_TOOL_DIR/jar/dropseq.jar \
	DigitalExpression I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.bam \
	O=$WD/dge/"$OUT"_total_dge.txt \
	EDIT_DISTANCE=1 \
	NUM_CORE_BARCODES=$CELL \
	SUMMARY=$WD/dge/"$OUT"_total_dge.summary.txt

	rm $WD/dge/topcellbcd.txt
	i=0
	while IFS=$'\t' read -r -a value
	do
		((i++))
		test $i -eq 1 && continue
		test $i -eq 2 && continue
		test $i -eq 3 && continue
		echo "${value[0]}">>$WD/dge/topcellbcd.txt
	done < $WD/dge/"$OUT"_total_dge.summary.txt
		
		
	java -jar -Xmx8g $DROP_SEQ_TOOL_DIR/jar/dropseq.jar \
	DigitalExpression I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_hg.bam \
	O=$WD/dge/"$OUT"_human_dge.txt \
	EDIT_DISTANCE=1 \
	CELL_BC_FILE=$WD/dge/topcellbcd.txt \
	SUMMARY=$WD/dge/"$OUT"_human.summary.txt


	java -jar -Xmx8g $DROP_SEQ_TOOL_DIR/jar/dropseq.jar \
	DigitalExpression I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean_mm.bam \
	O=$WD/dge/"$OUT"_mouse_dge.txt \
	EDIT_DISTANCE=1 \
	CELL_BC_FILE=$WD/dge/topcellbcd.txt \
	SUMMARY=$WD/dge/"$OUT"_mouse.summary.txt

	echo "Step 13 : Data combine using summary"
	$INHOUSE_SCRIPT_DIR/data_combine_using_summary.py \
	$WD/dge/"$OUT"_human.summary.txt \
	$WD/dge/"$OUT"_mouse.summary.txt \
	$WD/dge/"$OUT"_mix2_dge.summary.txt
	

else
	echo "Error : unspecified species identity"
fi

	
echo "Step 14 : BAM Tag Histogram"
$DROP_SEQ_TOOL_DIR/BAMTagHistogram \
I=$WD/star_align/"$OUT"_star_gene_exon_tagged_clean.bam \
O=$WD/dge/"$OUT"_cell_readcounts.txt \
TAG=XC

echo "Step 15 : Making summary files"

mkdir $WD/summary

$INHOUSE_SCRIPT_DIR/used_ratio.py \
$WD/dge/"$OUT"_cell_readcounts.txt \
$CELL \
$WD/summary/"$OUT"_used_read_ratio.txt


echo "Step 16 : Delete files"
rm -r $WD/pre-process \
$WD/star_align/"$OUT"_star_gene_exon_tagged.bam \
$WD/star_align/accepted_hits_sorted.bam \
$WD/star_align/"$OUT"_unaligned_mc_tagged_polyA_filtered.fastq \
$WD/star_align/"$OUT"_merged.bam \
$WD/star_align/Aligned.out.sam

