#inspired by https://url.au.m.mimecastprotect.com/s/md3sCq71grH8lE0DuEsKKiEY7ks?domain=biopython-tutorial.readthedocs.io
import gzip, sys

from Bio import SeqIO, bgzf

def trim_adaptors(records, adaptor, min_len):
	"""Trims perfect adaptor sequences, checks read length.

	This is a generator function, the records argument should
	be a list or iterator returning SeqRecord objects.
	"""

	len_adaptor = len(adaptor) #cache this for later
	for record in records:
		len_record = len(record) #cache this for later
		if len(record) < min_len:
			#Too short to keep
			continue
		index = record.seq.find(adaptor)
		if index == -1:
			#adaptor not found, so won't trim
			yield record
		elif index + len_adaptor - 1 >= min_len:
			#after trimming this will still be long enough
			yield record[:index + len_adaptor - 1]


infile = sys.argv[1]	#'B93449_AACAATG_psti-mspi.R1.fastq.gz'
infh = gzip.open(infile, 'rt')

outfile = sys.argv[2]	#'B93449_AACAATG_psti-mspi.R1.fastq.bgz'

original_reads = SeqIO.parse(infh, "fastq")
trimmed_reads = trim_adaptors(original_reads, "CTGCAG", 40)
trimmed_reads = trim_adaptors(trimmed_reads, "CCGG", 40)

with bgzf.BgzfWriter(outfile, "wb") as outgz:	#bgzf.BgzfWriter("test.fastq.bgz", "wb")
	count = SeqIO.write(trimmed_reads, handle=outgz, format="fastq")
	print("Saved %i reads" % count)

