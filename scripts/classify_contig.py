# usage: program contig.fa peak.bed ref_path
import sys, os
from Bio import SeqIO

if len(sys.argv) < 4:
	ref_loc = '/compbio/data/bwaIndex/ryang/hg19'
else:
	ref_loc = sys.argv[3]
mapq_cutoff = 1 

# mapping contigs to reference genome
#os.system('bwa bwasw -b 33 -r 17 -q 50 -w 500 -z 100 -t 6 '+ref_loc+' '+sys.argv[1]+' >tmp.sam')
os.system('bwa bwasw -t 6 '+ref_loc+' '+sys.argv[1]+' >tmp.sam')
os.system('samtools view -bS tmp.sam >tmp.bam')
os.system('samtools sort tmp.bam tmp.sorted')
os.system('samtools index tmp.sorted.bam')
os.system('bamToBed -i tmp.sorted.bam > tmp.bed')

# get mapped contigs set
os.system('sort -n -r -k 5,5 tmp.bed > tmp.sorted.bed')
ofp = open('tmp.uniqname.bed','w')
ifp = open('tmp.sorted.bed')
mapped_name = set()
for line in ifp:
	item = line.rstrip().split()
	name = item[3]
	if name in mapped_name:
		continue
	else:
		mapped_name.add(item[3])
		print >> ofp, line.rstrip()
ifp.close()
ofp.close()


# get repeat mappped contigs set
ofp3 = open(sys.argv[1]+'.repeat.bed','w')
repeat_name = set()
ifp = open('tmp.uniqname.bed')
ofp = open('tmp.reliable.bed','w')
for line in ifp:
	item = line.rstrip().split()
	if int(item[4]) < mapq_cutoff:
		repeat_name.add(item[3])
		print >> ofp3, line.rstrip()
	else:
		print >> ofp, line.rstrip()
ifp.close()
ofp.close()
ofp3.close()


# calculate coverage between reliable mapped contigs and peaks
os.system('coverageBed -a '+sys.argv[2]+' -b tmp.reliable.bed >tmp.coverage.bed')
os.system('sort -n -r -k 5,5 tmp.coverage.bed > tmp.coverage.sorted.bed')

reliable_name = set() # newly found peak region from reliable mapped contigs
peak_name = set() # contigs has overlap with peaks from peak caller
ifp = open('tmp.coverage.sorted.bed')
ofp1 = open(sys.argv[1]+'.reliable_new.bed','w')
ofp2 = open(sys.argv[1]+'.reliable_peak.bed','w')
chr_name = ['chr'+str(i) for i in range(1,23)]
chr_name.extend(['chrX','chrY'])
for line in ifp:
	item = line.rstrip().split()
	if item[0] not in chr_name:
		continue
	elif item[6] == '0':
		reliable_name.add(item[3])
		print >> ofp1, line.rstrip()
	else:
		peak_name.add(item[3])
		print >> ofp2, line.rstrip()
ifp.close()
ofp1.close()
ofp2.close()

# output all three categories contigs: unmapped, repeat, and reliable new in fasta format
ofp1 = open(sys.argv[1]+'.unmapped.fa','w')
ofp2 = open(sys.argv[1]+'.repeat.fa','w')
ofp3 = open(sys.argv[1]+'.reliable_new.fa','w')
ofp4 = open(sys.argv[1]+'.reliable_peak.fa','w')
for record in SeqIO.parse(sys.argv[1],'fasta'):
	if record.id not in mapped_name:
		print >> ofp1, '>'+record.id
		print >> ofp1, record.seq
	if record.id in repeat_name:
		print >> ofp2, '>'+record.id
		print >> ofp2, record.seq
	if record.id in reliable_name:
		print >> ofp3, '>'+record.id
		print >> ofp3, record.seq
	if record.id in peak_name:
		print >> ofp4, '>'+record.id
		print >> ofp4, record.seq
ofp1.close()
ofp2.close()
ofp3.close()
ofp4.close()
os.system('rm tmp*')
os.system('cat '+sys.argv[1]+'.reliable_new.fa '+sys.argv[1]+'.repeat.fa >'+sys.argv[1]+'.mapped_new.fa')
os.system('cat '+sys.argv[1]+'.reliable_new.bed '+sys.argv[1]+'.repeat.bed >'+sys.argv[1]+'.mapped_new.tmp')
os.system('cut -f1-6 '+sys.argv[1]+'.mapped_new.tmp >'+sys.argv[1]+'.mapped_new.bed')
os.remove(sys.argv[1]+'.mapped_new.tmp')
