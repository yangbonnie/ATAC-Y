# https://www.jianshu.com/p/09e05bcd6981
##1.filter adapter
cd /md01/yangjw28/data/Fibroblasts/9sample20220523
input_dir=/md01/yangjw28/data/Fibroblasts/9sample20220523/rawdata/220510_A00682_0759_BHGF2TDSX3/
output_dir=/md01/yangjw28/data/Fibroblasts/9sample20220523/fastp/
for i in {D1_L2_804F0002203,D3_L2_805F0002203,N3_L2_808F0004094,N4_L2_809F0004094,K5_L2_801F0002203,K7_L2_802F0002203,K9_L2_803F0002203};
do
      fastp --detect_adapter_for_pe \
	 -i ${input_dir}$i.R1.fastq.gz \
	 -I ${input_dir}$i.R2.fastq.gz \
	 -o ${output_dir}$i.R1.fastq \
	 -O ${output_dir}$i.R2.fastq
done

##2.mapping
#for file in `ls  ../Raw_data/merge_file/*_R1.fq.gz.trimmed.gz`
#do
#file2=`echo $file |sed 's/R1/R2/g'` 
#echo $file
#echo $file2
#name=`basename $file`
#output=`echo $name| awk -F_ '{print $1}'`
#echo $output
#bowtie2   -p 15  -x /home/jinxu/DB/mmu10/mm10_allchr_bowtie2Index/mm10 -1 $file -2 $file2 -S  $output.sam 1>$output.mapping.log 2>$output.mapping.err 
#done
cd /md01/yangjw28/data/Fibroblasts/9sample20220523
index_dir=/public/home/yangjw28/genomeDB/human-hg38/
input_dir=/md01/yangjw28/data/Fibroblasts/9sample20220523/fastp/
output_dir=/md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/

for i in {D1_L2_804F0002203,D3_L2_805F0002203,D4,N2,N3_L2_808F0004094,N4_L2_809F0004094,K5_L2_801F0002203,K7_L2_802F0002203,K9_L2_803F0002203};
do
	/md01/zhangwx/biosoft/bowtie2/bowtie2 -p 16 -k 10\
	--very-sensitive \
	-x ${index_dir}hg38 \
	-1 ${input_dir}$i.R1.fastq \
	-2 ${input_dir}$i.R2.fastq | /md01/shipy3/software/samtools-1.12/samtools view -F 4 -u | /md01/shipy3/software/samtools-1.12/samtools sort -@ 6 \
	-o ${output_dir}$i.sort.bam
done

##3.paired,rm chrM,sort,rmduplication
#for file in `ls *.sam`
#do
#echo $file
#file2=`echo $file|sed 's/\.sam//g'`
#echo $file2
#{
#awk '$3!="chrM"' $file | samtools view -S -b -F 4  -   -o $file2.bam  
#samtools sort $file2.bam $file2.sort
#&
#done 

#后台运行
cd /md01/yangjw28/data/Fibroblasts/9sample20220523/mapping
input_dir=/md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/rm_chrM_duplic/
for file in `ls *.bam`
do
	{	
		#samtools view -h -@ 9 $file |awk '$3!="chrM"'|samtools view -@ 9 -S -b -F 4 > /md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/rm_chrM3/$file.rmchrM3.bam   
		samtools sort  -@ 6 /md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/rm_chrM3/$file.rmchrM3.bam > /md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/rm_chrM3/$file.rmchrM3.bam.sort	
		#samtools view -h -f 2 -q 30 -@ 9 $file |grep -v chrM| samtools sort -@ 6 -o - > ${input_dir}$file.sort.bam
		java -jar /md01/zhangwx/miniconda3/bin/picard.jar MarkDuplicates -I ${input_dir}$file.rmchrM3.bam.sort \
					-O ${input_dir}$file.rmchrM3.bam.sort.markdup.bam -M ${input_dir}$file.rmchrM3.bam.sort.markdup.txt -REMOVE_DUPLICATES true 
	}& 
done

#wait



##4.merge.bam has beedn sorted and rm chrM 
cd /md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/rm_chrM3
for file in `ls *.sort`
do
{
genomeCoverageBed -bg -split -i qbam $file -g /public/home/yangjw28/genomeDB/human-hg38/hg38.chrom.sizes > $file.bedGraph 
/public/home/yangjw28/software/norm_bedGraph.pl $file.bedGraph $file.bedGraph.norm
/public/home/yangjw28/software/bedGraphToBigWig $file.bedGraph.norm   /public/home/yangjw28/genomeDB/human-hg38/hg38.chrom.sizes   $file.bedGraph.norm.bw 
}& 
done

#wait

##5.call peaks(MACS2)
cd /md01/yangjw28/data/Fibroblasts/9sample20220523/mapping/rm_chrM3
for file in `ls *.markdup.bam`
do
{
#bam2bed  $file # using good paired-end reads only.
bedtools bamtobed -i $file > $file.bed
bedfile=`echo $file|sed 's/bam/bed/g'`
outfile=`echo $file |sed 's/\.merge\.bam//g'`
path=`pwd`
macs2  callpeak -t  $file.bed -f BED  -g hs --outdir $path  -q 0.01 -n $outfile --nomodel  --shift 0 > $file.macs2.log 
} 
done


