Introduction
------------
ChIPdenovo is implemented by Python and it is designed to carry out de novo assembly of ChIP sequencing reads.

Pre-installtalation
-------------------
A. Programming environment and NGS tools
	1) Python version 2.5 or later (http://www.python.org/)
	2) R version 2.7 or later (http://www.r-project.org/)
	3) BWA aligner (http://bio-bwa.sourceforge.net/)
	4) SAMtools (http://samtools.sourceforge.net/)
	5) BEDtools (http://code.google.com/p/bedtools/)

B. Python packages:
	1) HTSeq (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
	2) Biopython (http://biopython.org/wiki/Main_Page)
	3) SciPy (http://www.scipy.org/)
	4) Numpy (http://www.numpy.org/)

Running ChIPdenovo
------------------
#### Command-line usage:
    python chipdenovo.py -i <filename.fasta> [opts]
#### Options:
 -K <int>    :kmer length (default:25, max=32)
 -L <int>    :min contig length to be reported (default:100)
 -r <int>    :min read coverage for contigs to be reported (default:2)
 -l <int>    :min merge suffix/prefix length (default:20)
 -h          :produce this menu
#### Example:
    python chipdenovo.py -i sample_data/lncap.simulated.fa -L 75 -l 10
Example data are provide in the directory *sample_data/*

Output
------
1. *.contig.final.fa is the final assembled contigs from ChIPdenovo
2.  directory *output/* stores all the intermediate resulst:
  2.1 *.contig is the output from Inchworm
  2.2 *.contig.merged.fa is generated by merging overlapped contigs from *.contig by at length n, where n is set by -l option in ChIPdenovo command line.
  2.3 *.contig.rename.fa is rename the headers in *.contig by using '_' instead of space

Downstream analysis
-------------------
1. mapping assembled contigs and classifying contigs based on their overlap with peak regions. 
Usage:    python classify_contig.py contig.fa peak.bed ref_path
*contig.fa is the assembled contigs output from ChIPdenovo
*peak.bed is the peaks identified by a peak calling tool in bed format.
*ref_path is the path for indexed reference geneome files
*classify_contig.py* will output classified contigs into four categories:.
*.unmapped.fa - the contigs can not be aligned to reference genome
*.repeat.fa/.bed - the contigs can not be uniquely aligned to reference genome
*.reliable_new.fa/.bed - the contigs reliably mapped to reference genome and excluded from identified peaks
*.reliable_peak.fa/.bed - the contigs reliably mapped to reference genome and included in identified peaks
*.mapped_new.bed/.bed - all the mapped novel peak regions identify by ChIPdenovo, the combined result form *repeat.bed and *reliable_new.bed

Contact us
----------
Questions, suggestions, comments, etc?
Author: Rendong Yang
Send email to cauyrd@gmail.com

Referencing ChIPdenovo
----------------------
updating soon!
