# script to get and quantify junction reads around introns
# input is intron bed file (generated by splicingrates_intron.py), readlength, bam file
# output is one bed6 file per desired intron type

# Usage: python splicingrates_reads.py --help

import sys
import pdb
import os
import subprocess
import argparse
import collections
from os.path import isfile, join
import gzip
import getopt

from startsites import readBam, startsites_single, startsites_paired

def getRegions(intronbedname, readlength, overlap):
	intronbed = open(intronbedname)
	outname = intronbedname[:-4]
	outfile_exonup = open(outname +'_'+ str(readlength) +'nt_exonup_ee.bed', 'w')
	outfile_exondown = open(outname +'_'+ str(readlength) +'nt_exondown_ee.bed', 'w')
	outfile_intron = open(outname +'_'+ str(readlength) +'nt_intron_ie.bed', 'w')
	while('TRUE'):
		intronline = intronbed.readline().strip().split('\t')
		if not intronline or intronline[0] == '': break
		# get strand
		chrom = intronline[0]
		name = intronline[3]
		strand = intronline[5]		
		if strand == "+":
			# exon up
			exonup_start = int(intronline[1])-10
			exonup_end = int(intronline[1])
			# exon down
			exondown_start = int(intronline[2])
			exondown_end = int(intronline[2])+overlap
			# intron
			intron_start = int(intronline[2])-(readlength - overlap)
			intron_end = int(intronline[2])-overlap
		if strand == "-":
			# exon up
			exonup_start = int(intronline[2])
			exonup_end = int(intronline[2])+overlap
			# exon down
			exondown_start = int(intronline[1])-overlap
			exondown_end = int(intronline[1])
			# intron
			intron_start = int(intronline[1])+overlap
			intron_end = int(intronline[1])+(readlength-overlap)
		outfile_exonup.write(chrom +'\t'+ str(exonup_start) +'\t'+ str(exonup_end) +'\t'+ name +'\t.\t'+ strand +'\n') 
		outfile_exondown.write(chrom +'\t'+ str(exondown_start) +'\t'+ str(exondown_end) +'\t'+ name +'\t.\t'+ strand +'\n') 
		outfile_intron.write(chrom +'\t'+ str(intron_start) +'\t'+ str(intron_end) +'\t'+ name +'\t.\t'+ strand +'\n') 

def JunctionBam_SE(bam, bamname):
	#bamname = bam[:-4]
	outfile = bamname + '_junctions_read1.bam'   
	cmd_header = "samtools view -H " + bam + " > " + bamname + "_junctions_header.txt"
	cmd_junc = "samtools view -F 256 " + bam + " | awk '{if ($6 ~/N/) {print $0}}' | cat " + bamname + "_junctions_header.txt - | samtools view -bS - | samtools sort - -T " \
               + outfile + " -o " + outfile
	cmd_index = "samtools index " +outfile
	for command in (cmd_header, cmd_junc, cmd_index):
		 subprocess.call(command, shell=True)

def JunctionBam_PE(bam, bamname):
    #bamname = bam[:-4]
    cmd_header = "samtools view -H " + bam + " > " + bamname + "_junctions_header.txt"
    # read 1
    outfile1 = bamname +"_junctions_read1.bam"
    cmd_junc1 = "samtools view -f 64 -F 256 " + bam + " | awk '{if ($6 ~/N/) {print $0}}' | cat " + bamname + "_junctions_header.txt - | samtools view -bS - | samtools sort - -T " \
               + outfile1 + " -o " + outfile1
    cmd_index1 = "samtools index " +outfile1
    for command in (cmd_header, cmd_junc1, cmd_index1):
        subprocess.call(command, shell=True)
    # read 2
    outfile2 = bamname +"_junctions_read2.bam"
    cmd_junc2 = "samtools view -f 128 -F 256 " + bam + " | awk '{if ($6 ~/N/) {print $0}}' | cat " + bamname + "_junctions_header.txt - | samtools view -bS - | samtools sort - -T " \
               + outfile2 + " -o " + outfile2
    cmd_index2 = "samtools index " +outfile2
    for command in (cmd_junc2, cmd_index2):
        subprocess.call(command, shell=True)

def getCoverage_intron(bed, juncbam, readtype, readstrand):
# args.outbed +'_startsites_read1.bed.gz
	juncreadfile1 = juncbam + '_startsites_read1.bed.gz'
	juncreadfile2 = juncbam + '_startsites_read2.bed.gz'
	if readstrand == 'fr-unstrand':
		# single bam file, no strand
		print("NOTE: unstranded reads")
		command = "coverageBed -a " + bed + " -b " + juncreadfile1 + " | gzip > " + juncbam + "_iejunc.bed.gz"
		subprocess.call(command, shell=True)
	if readtype == 'single' and readstrand == 'fr-firststrand':
		# single bam file, read1=-S
		print("NOTE: fr-firststrand SE reads")
		command = "coverageBed -a " +  bed + " -b " + juncreadfile1 + " -S | gzip > " + juncbam + "_iejunc.bed.gz"
		subprocess.call(command, shell=True)
	if readtype == 'single' and readstrand == 'fr-secondstrand':
		# single bam file, read1=-s
		print("NOTE: fr-secondstrand SE reads")
		command = "coverageBed -a " + bed + " -b " + juncreadfile1 + " -s | gzip > " + juncbam + "_iejunc.bed.gz"
		subprocess.call(command, shell=True)
	if readtype == 'paired' and readstrand == 'fr-firststrand':
		# two bams, read1=-S, read2=-s-
		print("NOTE: fr-firststrand PE reads")
		read1command = "coverageBed -a " + bed + " -b " + juncreadfile1 + " -S > " + juncbam + "_iejunc_read1.bed"
		read2command = "coverageBed -a " + bed + " -b " + juncreadfile2 + " -s > " + juncbam + "_iejunc_read2.bed"
		joincommand = "cat " + juncbam + "_iejunc_read*.bed | gzip > " +juncbam +"_iejunc.bed.gz"
		rmcommand = "rm -f " + juncbam + "_iejunc_read*.bed"
		for command in (read1command, read2command, joincommand, rmcommand):
			subprocess.call(command, shell=True)
	if readtype == 'paired' and readstrand == 'fr-secondstrand':
		# two bams, read1=-s, read2=-S
		print("NOTE: fr-secondstrand PE reads")
		read1command = "coverageBed -a " + bed + " -b " + juncreadfile1 + " -s > " + juncbam + "_iejunc_read1.bed"
		read2command = "coverageBed -a " + bed + " -b " + juncreadfile2 + " -S > " + juncbam + "_iejunc_read2.bed"
		joincommand = "cat " + juncbam + "_iejunc_read*.bed | gzip > " +juncbam +"_iejunc.bed.gz"
		rmcommand = "rm -f " + juncbam + "_iejunc_read*.bed"
		for command in (read1command, read2command, joincommand, rmcommand):
			subprocess.call(command, shell=True)

def getCoverage_exon(bed, juncbam, readtype, readstrand, region):
	# outname + '_junctions_read1.bam'
	juncreadfile1 = juncbam + '_junctions_read1.bam'
	juncreadfile2 = juncbam + '_junctions_read2.bam'
	if readstrand == 'fr-unstrand':
		# single bam file, no strand
		print("NOTE: unstranded reads")
		command = "intersectBed -abam " + juncreadfile1 + " -b " + bed + " -bed -F 0.95 -wo  | gzip > " + juncbam + "_"+ region +"junc.bed.gz"
		subprocess.call(command, shell=True)
	if readtype == 'single' and readstrand == 'fr-firststrand':
		# single bam file, read1=-S
		print("NOTE: fr-firststrand SE reads")
		command = "intersectBed -abam " + juncreadfile1 + " -b " + bed + " -bed -wo -F 0.95 -S | gzip > " + juncbam + "_"+ region +"junc.bed.gz"
		subprocess.call(command, shell=True)
	if readtype == 'single' and readstrand == 'fr-secondstrand':
		# single bam file, read1=-s
		print("NOTE: fr-secondstrand SE reads")
		command = "intersectBed -abam " + juncreadfile1 + " -b " + bed + " -bed -wo -F 0.95 -s | gzip > " + juncbam + "_"+ region +"junc.bed.gz"
		subprocess.call(command, shell=True)
	if readtype == 'paired' and readstrand == 'fr-firststrand':
		# two bams, read1=-S, read2=-s-
		print("NOTE: fr-firststrand PE reads")
		read1command = "intersectBed -abam " + juncreadfile1 + " -b " + bed + " -bed -wo -F 0.95 -S > " + juncbam + "_"+ region +"junc_read1.bed"
		read2command = "intersectBed -abam " + juncreadfile2 + " -b " + bed + " -bed -wo -F 0.95 -s > " + juncbam + "_"+ region +"junc_read1.bed"
		joincommand = "cat " + juncbam + "_"+ region +"junc_read*.bed | gzip > " +juncbam +"_"+ region +"junc.bed.gz"
		rmcommand = "rm -f " + juncbam + "_"+ region +"junc_read*.bed"
		for command in (read1command, read2command, joincommand, rmcommand):
			subprocess.call(command, shell=True)
	if readtype == 'paired' and readstrand == 'fr-secondstrand':
		# two bams, read1=-s, read2=-S
		print("NOTE: fr-secondstrand PE reads")
		read1command = "intersectBed -abam " + juncreadfile1 + " -b " + bed + " -bed -wo -F 0.95 -s > " + juncbam + "_"+ region +"junc_read1.bed"
		read2command = "intersectBed -abam " + juncreadfile2 + " -b " + bed + " -bed -wo -F 0.95 -S > " + juncbam + "_"+ region +"junc_read1.bed"
		joincommand = "cat " + juncbam + "_"+ region +"junc_read*.bed | gzip > " +juncbam +"_"+ region +"junc.bed.gz"
		rmcommand = "rm -f " + juncbam + "_"+ region +"junc_read*.bed"
		for command in (read1command, read2command, joincommand, rmcommand):
			subprocess.call(command, shell=True)

def readExonFile(filename):
        # bedfile format:
        # chr, intron start, intron end, '.', intron name, strand
        if filename[-2:] == "gz": bedfile = gzip.open(filename, 'rt')
        else: bedfile = open(filename)
        beddict = collections.defaultdict(dict)
        while('TRUE'):
                bedline = bedfile.readline().strip().split('\t')
                if not bedline or bedline[0] == '': break
                readstart = bedline[1]
                if bedline[5] == "-": readstart = bedline[2]
                beddict[bedline[15]][bedline[3]] = readstart
        return(beddict)

def readIntronFile(filename):
        if filename[-2:] == "gz": bedfile = gzip.open(filename, 'rt')
        else: bedfile = open(filename)
        beddict = collections.defaultdict(dict)
        while('TRUE'):
                bedline = bedfile.readline().strip().split('\t')
                if not bedline or bedline[0] == '': break
                beddict[bedline[3]] = bedline[6]
        return(beddict)

def combineRegions(introndict, exonupdict, exondowndict, juncbam):
	outfilename = juncbam + "_junctionCombo.coverage.gz"
	outfile = gzip.open(outfilename, 'wt')
	outfile.write('intron\tie_count\tee_count\n')
#	for intron in exonupdict.keys():
	for intron in introndict.keys():
		overlapboth = list(set(exonupdict[intron].keys()) & set(exondowndict[intron].keys()))
		eecount = len(overlapboth)
		outfile.write(str(intron) +'\t'+ str(introndict[intron]) +'\t'+ str(eecount) +'\n')
	outfile.close()

def rmTempFiles(juncbam):
	cmd_startsites = "rm -f " + juncbam + "_startsites_read*.bed.gz"
	cmd_juncbam = "rm -f " + juncbam + "_junctions_read*.bam"
	cmd_bamindex = "rm -f " + juncbam + "_junctions_read*.bam.bai"
	cmd_header = "rm -f " + juncbam + "_junctions_header.txt"
	cmd_intronbed = "rm -f " + juncbam + '_iejunc.bed.gz'
	cmd_exonupbed = "rm -f " + juncbam + '_eeupjunc.bed.gz'
	cmd_exondownbed = "rm -f " + juncbam + '_eedownjunc.bed.gz'
	for command in (cmd_startsites, cmd_juncbam, cmd_bamindex, cmd_header):
		subprocess.call(command, shell=True)
	for command in (cmd_intronbed, cmd_exonupbed, cmd_exondownbed):
		subprocess.call(command, shell=True)



if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Get relevant intron-exon and exon-exon junction reads for splicing rate estimation.', add_help = True)
	parser.add_argument('--introns', metavar='x.introntype.bed', type=str, help='bed file with specific introns for which to estimate splicing rates', required=True)
	parser.add_argument('--readlength', type=int, help='length of reads', required=True)
	### region information
	group_region = parser.add_argument_group('region information')
	group_region.add_argument('--intronRegions', action='store_true', help='Calculate specific intronic regions to use for grabbing reads', required=False)
	group_region.add_argument('--overlap', type=int, metavar='', default=10, help='overlap of junction reads with region of interest (nt)', required=False)
	### read quantification
	group_reads = parser.add_argument_group('read quantification')
	group_reads.add_argument('--junctionReads', action='store_true', help='Quantify junction reads from specified introns', required=False)
	group_reads.add_argument('--bam', type = str, metavar='x.bam', help = 'bam from which to extract junction reads. required if --junctionReads', required=False, default="None")
	group_reads.add_argument('--outdir', type = str, metavar='', help = 'directory for temporary and final files. required if --junctionReads', required=False, default="None")
	group_reads.add_argument('--readtype', type=str, help = 'type of read', default='paired', required=False, choices=['single','paired'])
	group_reads.add_argument('--readstrand', type=str, help = 'directionality of RNA-seq data', default='fr-firststrand', required=False, choices=['fr-unstrand','fr-firststrand','fr-secondstrand'])
	args = parser.parse_args()

	if not args.intronRegions and not args.junctionReads:
		sys.exit("ERROR! Need to specify either --intronRegions or --junctionRegions to perform a function with this script!")

	if args.intronRegions:
		print("Getting Intron Regions...")
		getRegions(args.introns, args.readlength, args.overlap)

	if args.junctionReads:
		print("Getting Junction Reads...")
		if args.bam == "None" or args.outdir == "None":
			sys.exit("ERROR! Need to include --bam and --outdir when running --junctionReads")
		# make output filename
		name = args.bam
		bamname = args.outdir + name.split('/')[-1]
		# get non-junction read start sites
		print("... read starts for non-split reads.")
		bamhere = readBam(args.bam)
		if args.readtype == "single" or args.readstrand == 'fr-unstrand':
			startsites_single(bamhere, bamname, args.readtype)
		if args.readtype == "paired":
			startsites_paired(bamhere, bamname, args.readtype)

		# get junction reads
		print("... split junction reads.")
		if args.readtype == 'single' or args.readstrand == 'fr-unstrand':
			JunctionBam_SE(args.bam, bamname[:-4])
		if args.readtype == 'paired' and args.readstrand != 'fr-unstrand':
			JunctionBam_PE(args.bam, bamname[:-4])

		# names of region files
		intronbed = args.introns[:-4] +'_'+ str(args.readlength) +'nt_intron_ie.bed'
		exonupbed = args.introns[:-4] +'_'+ str(args.readlength) +'nt_exonup_ee.bed'
		exondownbed = args.introns[:-4] +'_'+ str(args.readlength) +'nt_exondown_ee.bed'

		# get junction coverage
		juncbam = bamname[:-4]
		print("... coverage of junction regions...")
		print("...... intron-exon junction coverage.")
		getCoverage_intron(intronbed, juncbam, args.readtype, args.readstrand)
		print("...... junction coverage of upstream exons.")
		getCoverage_exon(exonupbed, juncbam, args.readtype, args.readstrand, "eeup")
		print("...... junction coverage of downstream exons.")
		getCoverage_exon(exondownbed, juncbam, args.readtype, args.readstrand, "eedown")

		# get intron-specific coverage and make final file
		print("... making combined coverage file.")
		exonupdict = readExonFile(juncbam + "_eeupjunc.bed.gz")
		exondowndict = readExonFile(juncbam + "_eedownjunc.bed.gz")
		introndict = readIntronFile(juncbam + "_iejunc.bed.gz")
		combineRegions(introndict, exonupdict, exondowndict, juncbam)

		# cleanup files
		rmTempFiles(bamname)




