import sys,HTSeq
import os,getopt
import numpy as np
from Bio import SeqIO
import itertools
from scipy.stats import scoreatpercentile
import shutil, glob

def _boolrelextrema(data, comparator,
                  axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`.

    Relative extrema are calculated by finding locations where
    comparator(data[n],data[n+1:n+order+1]) = True.

    Parameters
    ----------
    data: ndarray
    comparator: function
        function to use to compare two data points.
        Should take 2 numbers as arguments
    axis: int, optional
        axis over which to select from `data`
    order: int, optional
        How many points on each side to require
        a `comparator`(n,n+x) = True.
    mode: string, optional
        How the edges of the vector are treated.
        'wrap' (wrap around) or 'clip' (treat overflow
        as the same as the last (or first) element).
        Default 'clip'. See numpy.take

    Returns
    -------
    extrema: ndarray
        Indices of the extrema, as boolean array
        of same shape as data. True for an extrema,
        False else.

    See also
    --------
    argrelmax,argrelmin

    Examples
    --------
    >>> testdata = np.array([1,2,3,2,1])
    >>> argrelextrema(testdata, np.greater, axis=0)
    array([False, False,  True, False, False], dtype=bool)
    """

    if((int(order) != order) or (order < 1)):
        raise ValueError('Order must be an int >= 1')

    datalen = data.shape[axis]
    locs = np.arange(0, datalen)

    results = np.ones(data.shape, dtype=bool)
    main = data.take(locs, axis=axis, mode=mode)
    for shift in xrange(1, order + 1):
        plus = data.take(locs + shift, axis=axis, mode=mode)
        minus = data.take(locs - shift, axis=axis, mode=mode)
        results &= comparator(main, plus)
        results &= comparator(main, minus)
        if(~results.any()):
            return results
    return results

def argrelextrema(data, comparator,
                  axis=0, order=1, mode='clip'):
    """
    Calculate the relative extrema of `data`

    Returns
    -------
    extrema: ndarray
        Indices of the extrema, as an array
        of integers (same format as argmin, argmax

    See also
    --------
    argrelmin, argrelmax

    """
    results = _boolrelextrema(data, comparator,
                              axis, order, mode)
    if ~results.any():
        return (np.array([]),) * 2
    else:
        return np.where(results)

def mapping_reads_to_contig(contig,read):
	""" mapping raw chip-seq reads back to contigs"""
	os.system('bwa index -a is '+contig)
	os.system('bwa aln -t 8 '+contig+' '+read+' >'+contig+'.sample.sai')
	os.system('bwa samse '+contig+' '+contig+'.sample.sai '+read+' >'+contig+'.sample.sam')
	os.system('samtools view -bS '+contig+'.sample.sam >'+contig+'.sample.bam')
	os.system('samtools sort '+contig+'.sample.bam '+contig+'.sample.sorted')
	os.system('samtools index '+contig+'.sample.sorted.bam')
	os.system('bamToBed -i '+contig+'.sample.sorted.bam >'+contig+'.sample.bed')
	
def rename_fasta(contig):
	""" rename the header of each sequence in fasta file"""
	ifp = open(contig)
	ofp = open(contig+'.rename.fa','w')
	for line in ifp:
		if '>' in line:
			newline = '_'.join(line.rstrip().split())
			print >> ofp, newline
			continue
		print >> ofp, line.rstrip()
	ifp.close()
	ofp.close()

def merged_contigs(key, mydict):
	"""merged_contigs return the order of merged contig list"""
	if mydict[key][1] == 'null':
		return [key]
	return [key] + merged_contigs(mydict[key][1],mydict)

def usage():
	"""showing help information"""
	print 'Usage:'
	print '	python chipdenovo.py -i <filename.fasta> [opts]'
	print 'Opts:'
	print ' -K <int>	:kmer length (default:25, max=32)'
	print ' -L <int>	:min contig length to be reported (default:100)'
	print ' -r <int>	:min read coverage for contigs to be reported (default:2)'
	print ' -l <int>	:min merge suffix/prefix length (default:20)'
	print ' -h		:produce this menu'

if __name__ == "__main__":
	# parameters parsing.
	kmer = 25
	min_len = '100'
	read_cov = 2 
	min_merge_l = 20
	input = None
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:r:K:L:l:h')
	except getopt.GetoptError as err:
		print str(err)
		usage()
		sys.exit(2)
	for o,a in opts:
		if o == '-K':
			kmer = int(a)
		elif o == '-L':
			min_len = a
		elif o == '-i':
			input = a
		elif o == '-r':
			read_cov = int(a)
		elif o == '-l':
			min_merge_l = int(a)
		elif o == '-h':
			usage()
			sys.exit()
		else:
			assert False, "unhandled option"
	if not input:
		usage()
		sys.exit(2)

	#### STEP1: run inchworm for assembling the the raw reads

	os.system('inchworm --reads '+input+' --run_inchworm -K '+str(kmer)+' -L '+min_len+' --DS >'+input+'.contig')
	
	# rename contig fasta
	rename_fasta(input+'.contig')
	contig_dict = {}
	for read in SeqIO.parse(input+'.contig.rename.fa','fasta'):
		contig_dict[read.id] = read.seq
	print str(len(contig_dict))+' contigs assembled after inchworm step.'

	#### STEP2: merging contigs with common overlapped >= min_merge_l substring

	# look through generated contigs and generating suffix and prefix for each contig
	contig_num = 1e6
	for len_str in range(kmer-2,min_merge_l-1,-1):
		suffix = {}
		prefix = {}
		for key in contig_dict:
			try:
				prefix[contig_dict[key][0:len_str].tostring()].append(key+'#p')
			except KeyError:
				prefix[contig_dict[key][0:len_str].tostring()] = [key+'#p']
			try:
				prefix[contig_dict[key].reverse_complement()[0:len_str].tostring()].append(key+'#m')
			except KeyError:
				prefix[contig_dict[key].reverse_complement()[0:len_str].tostring()] = [key+'#m']
			try:
				suffix[contig_dict[key][-len_str:].tostring()].append(key+'#p')
			except KeyError:
				suffix[contig_dict[key][-len_str:].tostring()] = [key+'#p']
			try:
				suffix[contig_dict[key].reverse_complement()[-len_str:].tostring()].append(key+'#m')
			except KeyError:
				suffix[contig_dict[key].reverse_complement()[-len_str:].tostring()] = [key+'#m']
	
		# get overlapped substring
		common_str = set(suffix.keys()).intersection(set(prefix.keys()))
		linked = {}
		for each in common_str:
			name1 = suffix[each][0]
			name2 = prefix[each][0]
			if len(suffix[each]) == 1 and len(prefix[each]) == 1 and name1[0:-2] != name2[0:-2]:
				if name1 not in linked:
					linked[name1] = ['null',name2]
				else:
					linked[name1][1] = name2
				if name2 not in linked:
					linked[name2] = [name1,'null']
				else:
					linked[name2][0] = name1

		name_set = set()
		for key in linked:
			if linked[key][0] == 'null':
				merged_list = merged_contigs(key, linked)
				seq = ''
				name = ''
				for each in merged_list:
					name1 = each.split('#')
					if name1[1] == 'm':
						seq1 = contig_dict[name1[0]].reverse_complement()
					else:
						seq1 = contig_dict[name1[0]]
					seq += seq1[0:-len_str]
					name += name1[0]+'_+_'
				seq = seq + seq1[-len_str:]
				name = name[0:-3]
				names = ' '.join(sorted(name.split('_+_')))
				name = name + '_Kmer:_'+str(len_str)+'_len:_'+str(len(seq))
				if names in name_set:
					continue
				else:
					contig_dict[name] = seq
					name_set.add(names)
	
		for key in linked:
			name = key.split('#')
			try:
				del contig_dict[name[0]]
			except KeyError:
				pass
		if len(contig_dict) < contig_num:
			contig_num = len(contig_dict)
		else:
			break
	
	# generate merged contig file.
	ofp = open(input+'.contig.merged.fa','w')
	for key,value in contig_dict.iteritems():
		print >> ofp, '>'+key
		print >> ofp, value
	ofp.close()
	print str(len(contig_dict))+' contigs generated by merging step.'

	#### STEP3: filtering low coverage contigs by setting read coverage >= read_cov	

	# mapping raw reads to merged contigs
	mapping_reads_to_contig(input+'.contig.merged.fa',input)
	
	# calculating read coverage for merged contigs
	coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
	coverage2 = HTSeq.GenomicArray( "auto", stranded=True, typecode="i" )
	ifp = open(input+'.contig.merged.fa.sample.bed')
	for line in ifp:
		items = line.rstrip().split()
		if int(items[-2]) == 0:
			continue
		iv = HTSeq.GenomicInterval(items[0],int(items[1]),int(items[2]),'.')
		iv2 = HTSeq.GenomicInterval(items[0],int(items[1]),int(items[2]),items[-1])
		coverage[iv] += 1
		coverage2[iv2] += 1
	ifp.close()
	
	ofp = open(input+'.contig.merged.fa.coverage.fa','w')
	ofp1 = open(input+'.contig.merged.fa.coverage_p.fa','w')
	ofp2 = open(input+'.contig.merged.fa.coverage_m.fa','w')
	for item in coverage.chrom_vectors.keys():
		ctg_len = int(item.split('_')[-1])
		iv = HTSeq.GenomicInterval(item,0,ctg_len)
		cov = list(coverage[iv])
		avg_cov = sum(cov)/float(ctg_len)
		if avg_cov < read_cov:
			del contig_dict[item]
			continue
		iv_p = HTSeq.GenomicInterval(item,0,ctg_len,'+')
		iv_m = HTSeq.GenomicInterval(item,0,ctg_len,'-')
		cov1 = list(coverage2[iv_p])
		cov2 = list(coverage2[iv_m])
		print >> ofp, '>'+item
		print >> ofp, ' '.join(map(str,cov))
		print >> ofp1, '>'+item
		print >> ofp1, ' '.join(map(str,cov1))
		print >> ofp2, '>'+item
		print >> ofp2, ' '.join(map(str,cov2))
	ofp.close()
	ofp1.close()
	ofp2.close()
	print str(len(contig_dict))+' contigs left after filtering less than '+str(read_cov)+' read coverage.'

	#### STEP4: spliting possible worongly merged contigs based on the shape of read coverage curve	

	# smoothing read coverage for each generated contig using ksmooth in R
	fullpath = os.path.realpath(__file__)
	filename = sys.argv[0].split('/')[-1]
	path = fullpath.split(filename)[0]
	
	os.system('Rscript '+path+'smooth.r '+input+'.contig.merged.fa.coverage_p.fa '+input+'.contig.merged.fa.coverage_m.fa')
	
	# spliting merged contigs with double peak mode for each strand coverage curve
	in0 = iter(HTSeq.FastaReader(input+'.contig.merged.fa.coverage.fa'))
	in1 = iter(HTSeq.FastaReader(input+'.contig.merged.fa.coverage_p.fa.smoothed.fa'))
	in2 = iter(HTSeq.FastaReader(input+'.contig.merged.fa.coverage_m.fa.smoothed.fa'))
	for read0,read1,read2 in itertools.izip(in0,in1,in2):
		name = read1.name
		covp = map(float,read1.seq.split())
		covm = map(float,read2.seq.split())
		cov = map(int,read0.seq.split())
		loc_max_p = argrelextrema(np.array(covp),np.greater)
		loc_max_m = argrelextrema(np.array(covm),np.greater)
		if len(loc_max_p[0]) == 2 and len(loc_max_m[0]) == 2:
			new_array = np.hstack((loc_max_p[0],loc_max_m[0]))
			sort_index = np.argsort(new_array).tolist()
			if sort_index == [0,2,1,3] or sort_index == [2,0,3,1]:
				min_value = min(cov[new_array[sort_index[1]]:new_array[sort_index[2]]])
				min_index = cov[new_array[sort_index[1]]:new_array[sort_index[2]]].index(min_value)
				cut_index = min_index + new_array[sort_index[1]]
				contig_dict.update({'cut1_'+name:contig_dict[name][0:cut_index],'cut2_'+name:contig_dict[name][cut_index:-1]})
				del contig_dict[name]
				
	# generate final contig file after splitting.
	ofp = open(input+'.contig.final.fa','w')
	for key,value in contig_dict.iteritems():
		print >> ofp, '>'+key
		print >> ofp, value
	ofp.close()
	print str(len(contig_dict))+' contigs left after splitting step.'
	print 'contigs assembly is done!'
	os.mkdir('output')
	interfile = glob.glob('*.contig')
	interfile.extend(glob.glob('*.merged.fa*'))
	interfile.extend(glob.glob('*.rename.fa*'))
	for each in interfile:
		shutil.move(each,'output')
	print 'chipdenovo running is finished'
