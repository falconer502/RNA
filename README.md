# RNA
RNA-SEQ
read 密度图
1、首先对tophat或hisat软件比对产生的sam文件进行排序，并提取比对到每个scaffold的sam文件
将sam文件转换为bam文件
··samtools view -bS abc.sam > abc.bam
2、将bam文件进行排序
	samtools sort abc.bam abc.sort.bam
3、提取bam文件中比对到caffold1上的比对结果，并保存到sam文件格式
···$ samtools view abc.sort.bam scaffold1 > scaffold1.sam
4、运行python脚本进行画图

cumulative distribution curve 图绘制
将比对的生成的sam文件转换为bam文件，然后按名称进行排序，再转换为sam文件
samtools sort -n file.bam #sort bam by name
samtools view -h bamfile.bam > samfile.sam
1、首先提取exon、intergenic、intron、gene等gtf文件
python3 ex_gtf.py Artha....gff3 生成exon.gtf intron.gtf gene.gtf intergenic.gtf
2、利用htseq-count 对每种类型进行计数
htseq-count  ./t01.sorted.sam ./exon.gtf -i ID  -t CDS >CDS_count
3、利用脚本进行画图
python plot_curve.py inter_count intron_count CDS_count.. ..
