import subprocess
import os,shutil,getopt,sys
def write_log(log):
	f=open('log','a')
	print(log,file = f)
	f.close()
	
def exec(*arg) :   #write logfile 
	out=subprocess.Popen(arg[0],shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)#arg[0] command ,arg[1] comment
	write_log(arg[1])
	while out.poll()==None:       #check whether the process is commplete, if  have completed return 0
		write_log(out.stdout.read().decode('utf-8'))
		write_log("poll code is"+str(out.poll()))
	write_log("#######################################")

def mkdir(dirname):
	dirname=dirname.strip()
	isExists=os.path.exists(dirname)
	if not isExists:
		os.mkdir(dirname)
		return True
	else:
		print(dirname+'already exists!')
		#shutil.rmtree(dirname)
		#os.mkdir(dirname)
def extract_es(gtf):                   #create splicesites fileT
	
	cmd1='extract_exons.py '+gtf+'> exons.txt'
	cmd2='extract_splice_sites.py '+gtf+'>splicesites.txt'
	exec(cmd1,'creat exons.txt')
	exec(cmd2,'creat splicesites.txt')

def build_index(genome,gtf):
	extract_es(gtf)
	mkdir('index')
	cmd='hisat2-build -f %s --ss splicesites.txt --exon exons.txt -p 20 ./index/genome.index' %genome
	exec(cmd,'build index with hisat2')

#exec("ls -s","pirnt dir....")
def align_genome(genome,fa_path):   #align genome,input path of data files
	sequenced_file=os.listdir(fa_path)													#path of input file
	print(sequenced_file)  #command of genome alginment
	for i in sequenced_file:
		if os.path.isdir(fa_path+i)==True:
			file1,file2=os.listdir(fa_path+i)
			file1=fa_path+i+'/'+file1
			file2=fa_path+i+'/'+file2
			cmd='hisat2 -q -x ./index/genome.index -1 %s -2 %s -S ./temp/%s.sam -p 20 -t'%(file1,file2,i)
			exec(cmd ,"%s is algining against genome  using hisat2"%i) #execute commad
def stringtie(gtf,filepath):
	datalist_tmp=[i+'.sam' for i in datalist]
	write_log(datalist_tmp)
	if set(datalist_tmp).issubset(set(os.listdir(filepath))) and 'merged.gtf' in gtf :
		for i in datalist_tmp:
			name=i.split('.')[0]
			print(i)
			exec('stringtie ./temp/%s.sorted.bam -b ./results/%s -e -G ./merge/merged.gtf -C ./results_gtf/%s.cov_ref.gtf -p 20 -o ./results_gtf/%s.stringtie.out.gtf'%(i,name,name,name),'%s transcript assembly and qualntification for RNA-SEQ using stingtie'%(name))
			# stringite after emerging  gtf  #t01.sam.sorted.bam  -b t01.sam -G ./temp/merged.gtf -c t01.cov_ref.gtf -o t01.sam.stringtie.out.gtf
	else:	
		for i in datalist_tmp:
			name=i
			exec('samtools view -Su ./temp/%s | samtools sort > ./temp/%s.sorted.bam'%(name,name),'converted %s to bam format using samtools'%(name))
			exec('stringtie ./temp/%s.sorted.bam -G %s -p 20 -o ./temp/stringtie/%s.stringtie.out.gtf'%(name,gtf,name),'%s transcript assembly and qualntification for RNA-SEQ using stingtie'%(name))
			#t01.sam.sorted.bam
def merge(filepath,gtf):
	list=os.listdir(filepath)
	list=[filepath+i for i in list]
	listfile=' '.join(list)
	write_log(listfile)
	exec('stringtie --merge %s -G %s -o ./merge/merged.gtf'%(listfile,gtf),'merge all gtf.......')
	exec('rm %s'%(listfile),'delete gtf file.........')
#def ballgown():



def usuage():
	write_log('''
		g: input genome file
		f: input cleandata file path
		t: input gtf file''')



if __name__=="__main__":
	global datalist
	try:
		opts,args=getopt.getopt(sys.argv[1:],"hg:t:f:",["help","genome=","gtf=","filepath="])
		print("===============opt=====================")
		for name,value in opts:
			if name in ('-g',"--genome"):
				genome=value
				print('%-20s: %-10s'%('genome file is',value))
			if name in ('-t','-gtf'):
				gtf=value
				print('%-20s: %-10s'%('gtf file is',value))
			if name in ('-f','--filepath'):
				filepath=value
				datalist=os.listdir(filepath)
				print('%-20s: %-10s'%('fasta filepath are',value))
			if name in ('-h','-help'):
				usuage()
				sys.exit()
		write_log("===============opt=====================")
		write_log('start analysing..............')	
		
	except getopt.GetoptError:
		write_log("getopt error")
		usuage()
		sys.exit()
	
	
	write_log("setp1: creat temp directory.....")
	mkdir('temp')
	mkdir('results_gtf')
	mkdir('./temp/stringtie')
	mkdir('merge')
	mkdir('results')    #create a temporary  folder
	write_log('setp2: build index..............')
	build_index(genome,gtf)
	write_log('step3: algin genome.............')
	align_genome(genome,filepath) # align genome and remove the result to temp)
	write_log('step4:stringtie.................')
	stringtie(gtf,'./temp') 
	write_log('step5:merging...................')
	merge('./temp/stringtie/',gtf)
	write_log('step6:stringtie after merged....')
	stringtie('./merge/merged.gtf','./temp')
