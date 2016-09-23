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
                write_log(out.stderr.read().decode('utf-8'))
                print("poll code is"+str(out.poll()))
        print("#######################################")

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
        cmd='hisat2-build -f %s --ss splicesites.txt --exon exons.txt -p 7 ./index/setaria.index' %genome
        exec(cmd,'build index with hisat2')

#exec("ls -s","pirnt dir....")
def align_genome(genome,fa_path):   #align genome,input path of data files
        sequenced_file=os.listdir(fa_path)                                                                                                      #path of input file
        print(sequenced_file)  #command of genome alginment
        for i in sequenced_file:
                if os.path.isdir(fa_path+i)==True:
                        file1,file2=os.listdir(fa_path+i)
                        file1=fa_path+i+'/'+file1
                        file2=fa_path+i+'/'+file2
                        cmd='hisat2 -q -x ./index/setaria.index -1 %s -2 %s -S ./temp/%s.sam -t'%(file1,file2,i)
                        exec(cmd ,"%s algin genome  complete with hisat2"%i) #execute commad
def stringtie(gtf,filepath):
        for i in os.listdir(filepath):
                if 'sam' in i: #i==t01.sam t02.sam t01.sam.sorted.bam
                        name=i  #name==t01.sam.sorted.bam
                        if 'merged.gtf' in gtf :
                                if 'sorted.bam' in i:
                                        name=i.split('.')[0]
                                        print(i)
                                        exec('stringtie ./temp/%s -b ./results/%s -e -G ./merge/merged.gtf -C ./results/%s.cov_ref.gtf -p 7 -o ./results/%s.stringtie.out.gtf'%(i,name,name,name),'transcript assembly and qualntification for RNA-SEQ using stingtie ')
                                else:
                                        continue
                                        # stringite after emerging  gtf  #t01.sam.sorted.bam  -b t01.sam -G ./temp/merged.gtf -c t01.cov_ref.gtf -o t01.sam.stringtie.out.gtf
                        else:
                                exec('samtools view -Su ./temp/%s | samtools sort > ./temp/%s.sorted'%(name,name),'converted to bam format using samtools')
                                exec('stringtie ./temp/%s.sorted.bam -G %s -p 7 -o ./temp/stringtie/%s.stringtie.out.gtf'%(name,gtf,name),'transcript assembly and qualntification for RNA-SEQ using stingtie ')
                        #t01.sam.sorted.bam
def merge(filepath,gtf):
        list=os.listdir(filepath)
        list=[filepath+i for i in list]
        listfile=' '.join(list)
        print(listfile)
        exec('stringtie --merge %s -G %s -o ./merge/merged.gtf'%(listfile,gtf),'merge all gtf.......')
#def ballgown():



def usuage():
        print('''
                g: input genome file
                f: input cleandata file path
                t: input gtf file''')



if __name__=="__main__":

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
                                print('%-20s: %-10s'%('fasta filepath are',value))
                        if name in ('-h','-help'):
                                usuage()
                                sys.exit()
                print("===============opt=====================")
                print('start analysing..............')

        except getopt.GetoptError:
                print("getopt error")
                usuage()
                sys.exit()


        print("setp1: creat temp directory.....")
        mkdir('temp')
        mkdir('./temp/stringtie')
        mkdir('merge')
        mkdir('results')    #create a temporary  folder
        print('setp2: build index..............')
        build_index(genome,gtf)
        print('step3: algin genome.............')
        align_genome(genome,filepath) # align genome and remove the result to temp
        print('step4:stringtie.................')
        stringtie(gtf,'./temp')
        print('step5:merging...................')
        merge('./temp/stringtie/',gtf)
        print('step6:stringtie after merged....')
        stringtie('./merge/merged.gtf','./temp')
