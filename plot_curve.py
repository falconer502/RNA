import os,sys,random,math
import matplotlib.pyplot as plt
from matplotlib import style

files=sys.argv[1:]

def del_file(infile):
        gene_totalcount=0
        read_totalcount=0
        data=[]
	newdata=[]
	totalrpkm=0
        with open (infile,'r') as f:
                for line in f:
                        gene,count=line.split('\t')
			#print gene
			gene_len=gene.split(':')[1]
			
			if int(count)>0:
                        	gene_totalcount+=1
                        	read_totalcount+=int(count)
                        	data.append([gene,int(count),int(gene_len)]) 
		for i in data:
			gene,count,gene_len=i[:]
			
			rpkm=count/float(gene_totalcount*gene_len)*1000000000
			newdata.append([gene,rpkm])
			totalrpkm+=rpkm
			#print count,gene_totalcount, rpkm #gene,count,length
        return (newdata,gene_totalcount,totalrpkm)

def del_data(data,gene_totalcount,totalrpkm) :
        gene_percent=0.0
        rpkm_percent=0.0
        xlist=[]
        ylist=[]
	rpkmlist=[]
        data.sort(key=lambda x:-x[1])
        for d in data:
                gene_id,rpkm=d[:]
		rpkmlist.append(rpkm)
                gene_percent+=1/float(gene_totalcount)*100
                rpkm_percent+=rpkm/float(totalrpkm)
		
                xlist.append(gene_percent)
                ylist.append(rpkm_percent)
        return(xlist,ylist,rpkmlist)
color=['r','b','m','c','g','k','y']
col=0
rpkmlog=[]
for file in files:
	data,gene_totalcount,totalrpkm=del_file(file)
	xlist,ylist,rpkmlist=del_data(data,gene_totalcount,totalrpkm)
	rpkmlog.append([math.log(i,10) for i in rpkmlist])
	plt.figure(1)
	style.use('ggplot')
	plt.plot(xlist, ylist, color = color[col], label = file)
	#plt.figure(2)
	#plt.boxplot(rpkmlog)
	col=col+1
plt.figure(1)
plt.legend(loc='center right')
plt.xlim([-1,100])
plt.ylim([-0.02,1.05])
plt.figure(2)
bp=plt.boxplot(rpkmlog,patch_artist=True)
col=0
for box in bp['boxes']:
	box.set(color = color[col])

plt.show()
