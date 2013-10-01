# usage: program contig.fa ref_path
import sys, os
from Bio import SeqIO

if len(sys.argv) < 3:
	ref_loc = '/compbio/data/bwaIndex/ryang/hg19'
else:
	ref_loc = sys.argv[2]

mapq_cutoff = 1 
mapl_cutoff = 0.9

# mapping contigs to reference genome
#os.system('bwa bwasw -b 33 -r 17 -q 50 -w 500 -z 100 -t 6 '+ref_loc+' '+sys.argv[1]+' >tmp.sam')
os.system('bwa bwasw -t 6 '+ref_loc+' '+sys.argv[1]+' >tmp.sam')
os.system('samtools view -bS tmp.sam >tmp.bam')
os.system('samtools sort tmp.bam tmp.sorted')
os.system('samtools index tmp.sorted.bam')
os.system('bamToBed -i tmp.sorted.bam > tmp.bed')

# get mapped contigs set
os.system('sort -n -r -k 5,5 tmp.bed > tmp.sorted.bed')
ifp = open('tmp.sorted.bed')
mapped_name = set()
unique_name = set() 
ofp = open(sys.argv[1]+'.unique.bed','w')
for line in ifp:
	item = line.rstrip().split()
	name = item[3]
	if name in mapped_name:
		continue
	else:
		length1 = int(item[2]) - int(item[1])
		length2 = float(item[3].split('_')[-1])
		ratio = length1/length2
		if ratio >= mapl_cutoff:
			mapped_name.add(item[3])
			if int(item[4]) >=  mapq_cutoff:
				unique_name.add(item[3])
				print >> ofp, line.rstrip()
ifp.close()
ofp.close()

# output all three categories contigs: mapped, uniquely mapped and unmapped contigs
ofp1 = open(sys.argv[1]+'.unmapped.fa','w')
ofp2 = open(sys.argv[1]+'.unique.fa','w')
ofp3 = open(sys.argv[1]+'.mapped.fa','w')
for i,record in enumerate(SeqIO.parse(sys.argv[1],'fasta')):
	if record.id not in mapped_name:
		print >> ofp1, '>'+record.id
		print >> ofp1, record.seq
	if record.id in unique_name:
		print >> ofp2, '>'+record.id
		print >> ofp2, record.seq
	if record.id in mapped_name:
		print >> ofp3, '>'+record.id
		print >> ofp3, record.seq
ofp1.close()
ofp2.close()
ofp3.close()
os.system('rm tmp*')
print str(i+1)+'\tin total'
print str(len(mapped_name))+'\tmapped (identify>='+str(mapl_cutoff)+')'
print str(len(unique_name))+'\tuniquely mapped (mapQ>='+str(mapq_cutoff)+')'
