#!/bin/env python3


import sys

if sys.version_info[0] < 3:
    print("Package requires Python 3")
    sys.exit(1)


if "/vol/mauranolab/mauram01/hybridmice/dnase/src" not in sys.path:
       sys.path.insert(0, "/vol/mauranolab/mauram01/hybridmice/dnase/src")

import cigar
from sys import argv
import re
import argparse


import pysam


class SNV(object):
       def __init__(self, line):
              snpcoords = line.rstrip('\n').split('\t')
              self.mm10_chrom = snpcoords[0]
              self.mm10_chromStart = int(snpcoords[1])
              self.mm10_chromEnd = int(snpcoords[2])
              self.rsid = snpcoords[3]
              mm10_genotype = snpcoords[4].split("/")
              self.mm10_ref_allele = mm10_genotype[0]
              self.mm10_nonref_allele = mm10_genotype[1]
              self.denovo_chrom = snpcoords[6]
              self.denovo_chromStart = int(snpcoords[7])
              self.denovo_chromEnd = int(snpcoords[8])
              denovo_genotype = snpcoords[9].split("/")
              #NB denovo genotype is the same as mm10 genotype, so need to reverse alleles
              self.denovo_ref_allele = denovo_genotype[1]
              self.denovo_nonref_allele = denovo_genotype[0]
              self.distToNearestSNP = int(snpcoords[10])
       
       def __str__(self):
              return("Parsing SNP at mm10: " + self.mm10_chrom + ":" + str(self.mm10_chromStart) + "-" + str(self.mm10_chromEnd) + " " + self.rsid + " (" + self.mm10_ref_allele + "/" + self.mm10_nonref_allele + "; " + self.denovo_ref_allele + "/" + self.denovo_nonref_allele + ")")


def getSequenceAtSNV(pileupread):
       if pileupread is None:
              return None
       else:
              if pileupread.query_position is None:
                     #If pileup overlapped a deletion in read
                     return None
              else:
                     return pileupread.alignment.query_sequence[pileupread.query_position]

def printRead(prefix, pileupread, suffix, verbose):
       if verbose:
              try:
                     readname = pileupread.alignment.query_name
                     readname = re.sub(r"^(HISEQ\-|NS)[^:]+:", "", readname)
              
                     if pileupread.alignment.is_reverse:
                            read_strand = "-"
                     else:
                            read_strand = "+"
              
                     if pileupread.alignment.is_read1:
                            read_num = 1
                     else:
                            read_num = 2
              
                     read_mismatches = pileupread.alignment.get_tag("NM")
                     read_mapq = pileupread.alignment.mapping_quality
                     read_allele = getSequenceAtSNV(pileupread)
              
                     #Convert query_position to 1-based
                     print('\t%s%s R%s%s:%s:%s, %ibp, TLEN %i, flag %i, MAPQ %i, %s | Pos %i=%s; BaseQ %i; %imm' % (prefix, readname, str(read_num), suffix, pileupread.alignment.reference_start, read_strand, pileupread.alignment.query_length, pileupread.alignment.template_length, pileupread.alignment.flag, read_mapq, pileupread.alignment.cigarstring, getCycleNum(pileupread), read_allele, pileupread.alignment.query_qualities[pileupread.query_position], read_mismatches), file=sys.stderr)
              except:
                     #the "||" distinguishes failures like this
                     print('\t%s%s R%s%s:%s:%s, %ibp, TLEN %i, flag %i, MAPQ %i, %s || %imm' % (prefix, readname, str(read_num), suffix, pileupread.alignment.reference_start, read_strand, pileupread.alignment.query_length, pileupread.alignment.template_length, pileupread.alignment.flag, read_mapq, pileupread.alignment.cigarstring, read_mismatches), file=sys.stderr)


def isUnwantedChromosome(seqname):
       if re.search('hap', seqname) is not None or re.search('random', seqname) is not None or re.search('chrUn_', seqname) is not None or re.search('^scaffold', seqname) is not None or re.search('^C\d+', seqname) is not None:
              return True
       else:
              return False


def getCycleNum(pileupread):
       #Note that pileupread.query_position is 0-based
       #We return the cycle number starting from 1
       if pileupread.query_position is None:
              #If pileup overlapped a deletion in read
              return None
       elif pileupread.alignment.is_reverse:
              return pileupread.alignment.query_length - pileupread.query_position
       else:
              return pileupread.query_position + 1


def preprocessRead(pileupread, readHash, whichfile):
       #Add the read
       readHash[pileupread.alignment.query_name + str(pileupread.alignment.is_read1)] = pileupread


def getPreprocessedRead(reads, readname, whichfile):
       if readname in reads:
              pileupread = reads[readname]
              printRead("", pileupread, whichfile, args.verbose)
              return pileupread
       else:
              return None


def preValidateIndividualPileupRead(pileupread, args, curSNV):
       if pileupread is None:
              return None
       else:
              #Enforce filters that do not count as bias -- these reads will NOT be counted in totalReadsInclFailed
              if pileupread.query_position is not None and getCycleNum(pileupread) > args.cyclesToSkip and abs(pileupread.alignment.template_length) <= args.maxTemplateLength and pileupread.alignment.query_qualities[pileupread.query_position] >= args.minBaseQ and (not args.match_rsid or pileupread.alignment.query_name.split(';')[4]==curSNV.rsid):
                     return True
              else:
                     printRead("\tPrevalidate fail ", pileupread, "", args.verbose)
                     return False


def validateIndividualPileupRead(pileupread, read_allele, ref_allele, nonref_allele, numSNVs, args):
       #read_allele passed in purely for efficiency
       
       read_mismatches = pileupread.alignment.get_tag("NM")
       
       
       #Enforce filters that are indicative of bias -- these reads will be counted in totalReadsInclFailed
       #Used to set bitmask in pileup but has disappeared from docs and has no effect. Perhaps it was replaced by the new stepper param?
       
       pass_bitmask = (pileupread.alignment.flag & args.filter == 0) and (pileupread.alignment.flag & args.require == args.require)
       
       pass_mapq = pileupread.alignment.mapping_quality >= args.minMAPQ
       
       pass_mismatch = ((read_allele == ref_allele and read_mismatches <= args.permittedMismatches) or (read_allele == nonref_allele and read_mismatches <= args.permittedMismatches + numSNVs))
       
       pass_seq = re.search('[^ACGT]', pileupread.alignment.query_sequence) is None
       
       pass_cigar = re.search('[HSPDI]', pileupread.alignment.cigarstring) is None
       
       if pass_bitmask and pass_mapq and pass_seq and pass_mismatch and pass_cigar:
              return True
       else:
              if args.verbose:
                     #print("flag:", pileupread.alignment.flag, ", ", args.filter, ", ", args.require)
                     #print("seq:", pileupread.alignment.query_sequence)
                     print("\t\tValidate fail ", pileupread.alignment.query_name, ". bitmask=", pass_bitmask, ", mapq=", pass_mapq, ", mismatch=", pass_mismatch, ", seq=", pass_seq, ", cigar=", pass_cigar, sep="", file=sys.stderr)
              return False


class OtherSNVException(Exception):
       pass


def checkForOverlappingSNV(SNV, alignment, isRef, file1):
       if SNV is None:
              return 0
       
       if file1:
              chromStart = SNV.mm10_chromStart
              chromEnd = SNV.mm10_chromEnd
              ref_allele = SNV.mm10_ref_allele
              nonref_allele = SNV.mm10_nonref_allele
       else:
              chromStart = SNV.denovo_chromStart
              chromEnd = SNV.denovo_chromEnd
              ref_allele = SNV.denovo_ref_allele
              nonref_allele = SNV.denovo_nonref_allele
       
       #If there is a SNV overlapping by at least 1 bp
       if alignment.reference_start < chromEnd and alignment.reference_end > chromStart:
              #TODO Don't support multiple SNVs combined with indels/clipping yet
              if re.search('[HSPNDI]', alignment.cigarstring) is not None:
                     raise OtherSNVException
              
              #Assumes point variant
              offset = chromStart - alignment.reference_start
              readallele = alignment.query_sequence[offset]
              
              if args.verbose:
                     if file1:
                            whichfile = ", 1"
                     else:
                            whichfile = ", 2"

                     print("\t\tAnother SNV overlapping f", whichfile, ". isRef ", isRef, "; ", SNV.rsid, ": ", chromStart, " (", ref_allele, "/", nonref_allele, "). Offset ", offset, "=", readallele, ". ", sep="", file=sys.stderr)
              
              if isRef and readallele == ref_allele or not isRef and readallele == nonref_allele:
#                     if args.verbose:
#                            print("Overlap pass!!", file=sys.stderr)
                     return 1
              else:
                     raise OtherSNVException
       else:
              return 0


def validateDuallyAlignedPileupReads(file1pileupread, file2pileupread, prevSNV, curSNV, nextSNV, args):
       file1read_allele = getSequenceAtSNV(file1pileupread)
       file2read_allele = getSequenceAtSNV(file2pileupread)
       
       isRef = file1read_allele == curSNV.mm10_ref_allele
       if not isRef and file1read_allele != curSNV.mm10_nonref_allele:
              return False
       
       numSNVsfromFile1 = 1
       numSNVsfromFile2 = 1
       try:
              #NB we only examine the prev and next SNV, but anyway we don't permit >2 mm.
              numSNVsfromFile1 += checkForOverlappingSNV(prevSNV, file1pileupread.alignment, isRef, True)
              numSNVsfromFile1 += checkForOverlappingSNV(nextSNV, file1pileupread.alignment, isRef, True)
              
              numSNVsfromFile2 += checkForOverlappingSNV(prevSNV, file2pileupread.alignment, not isRef, False)
              numSNVsfromFile2 += checkForOverlappingSNV(nextSNV, file2pileupread.alignment, not isRef, False)
       except OtherSNVException as e:
              #print("Exceptionall!!", file=sys.stderr)
              return False
       if numSNVsfromFile1 != numSNVsfromFile2:
              return False
       
       #print("trying!!", numSNVsfromFile1, "/", numSNVsfromFile2, file=sys.stderr)
       
       file1read_mismatches = file1pileupread.alignment.get_tag("NM")
       file2read_mismatches = file2pileupread.alignment.get_tag("NM")
       
       return file1read_allele == file2read_allele and (isRef and file1read_mismatches + numSNVsfromFile1 == file2read_mismatches or not isRef and file1read_mismatches == file2read_mismatches + numSNVsfromFile2) and validateIndividualPileupRead(file1pileupread, file1read_allele, curSNV.mm10_ref_allele, curSNV.mm10_nonref_allele, numSNVsfromFile1, args) and validateIndividualPileupRead(file2pileupread, file2read_allele, curSNV.denovo_ref_allele, curSNV.denovo_nonref_allele, numSNVsfromFile2, args) and file1pileupread.query_position == file2pileupread.query_position and file1pileupread.alignment.is_reverse == file2pileupread.alignment.is_reverse


def computeDupKey(alignment):
       #Include reference sequence name in key for generality
       #NB SE read overlapping only one mate in PE pair will not be considered duplicates
       #TODO parse pileupread.cigartuples
       if alignment.is_reverse:
              fiveprime = alignment.reference_end - 1
       else:
              fiveprime = alignment.reference_start
       dupkey = str(alignment.reference_id) + ":" + str(fiveprime) + ":" + str(alignment.is_reverse)
       if alignment.is_paired:
              if not alignment.is_proper_pair:
                     raise Exception("Read must be properly paired " + alignment.query_name)
              #Reads are all properly paired so can infer strand of R1 as !R1
              if not alignment.is_reverse:
                     #Mate is reverse strand
                     next_cigar = alignment.get_tag("MC")
                     if cigar is not None:
                            next_cigarlength = cigar.aligned_length(cigar.parse(next_cigar))
                            fiveprime = alignment.next_reference_start + next_cigarlength - 1
                     else:
                            raise Exception("Require MC tag on alignment " + alignment.query_name)
              else:
                     fiveprime = alignment.next_reference_start
              
              dupkey = dupkey + str(alignment.next_reference_id) + ":" + str(fiveprime) + ":" + str(not alignment.is_reverse)
       return dupkey


def checkDupRead(pileupread1, pileupread2, dupReadHash):
       if pileupread1 is None and pileupread2 is None:
              return None
       
       if pileupread1 is not None:
              dupkey = computeDupKey(pileupread1.alignment)
       elif pileupread2 is not None:
              dupkey = computeDupKey(pileupread2.alignment)
       else:
              raise Exception("impossible")
       
       if not dupkey in dupReadHash:
              #Mark this location in case we see duplicate reads later on
              dupReadHash[dupkey] = True
              return False
       else:
              return True


def processSNV(prevSNV, curSNV, nextSNV):
       global lastChrom
       global totalChroms
       
       #Collect summary statistics
       global totalChroms
       global totalSNVsSkipped
       global totalSNVsCounted
       global grandTotalReadsSkipped
       global grandTotalReadsSkippedRef
       global grandTotalReadsSkippedNonref
       global grandTotalPassedReadsMarkedDup
       global grandTotalDupReadsRef
       global grandTotalDupReadsNonref
       global grandTotalReads
       global grandTotalReadsRef
       global grandTotalReadsNonref
       global grandTotalReadsInclFailed
       global grandTotalFailedReadsRef
       global grandTotalFailedReadsNonref
       global grandTotalFailedPairedReadsRef
       global grandTotalFailedPairedReadsNonref
       
       
       if curSNV.mm10_chrom != lastChrom:
              print(curSNV.mm10_chrom, file=sys.stderr)
              totalChroms += 1
              #Avoid comparing to SNV on wrong chromosome.
              #BUGBUG don't look ahead to check chromosome of nextSNV
              prevSNV = None
       lastChrom = curSNV.mm10_chrom
       
       if args.verbose: print(curSNV, file=sys.stderr)
       
       ###First, load and pre-process reads in pileup
       file1reads = dict()
       file2reads = dict()
       
       #NB only one column for a SNP (i think)
       for pileupcolumn in bamfile1.pileup(curSNV.mm10_chrom, curSNV.mm10_chromStart, curSNV.mm10_chromEnd, max_depth=100000, truncate=True, stepper="nofilter"):
              if pileupcolumn.nsegments >= args.minReads:
                     #Only examine pileup if there might be enough reads
                     for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                   preprocessRead(pileupread, file1reads, "f1")
       
       #Only do pileup of second bam if we enough reads have already passed
       if len(file1reads) >= args.minReads:
              for pileupcolumn in bamfile2.pileup(curSNV.denovo_chrom, curSNV.denovo_chromStart, curSNV.denovo_chromEnd, max_depth=100000, truncate=True, stepper="nofilter"):
                     if pileupcolumn.nsegments >= args.minReads:
                            #Only examine pileup if there might be enough reads
                            for pileupread in pileupcolumn.pileups:
                                   if not pileupread.is_del and not pileupread.is_refskip:
                                          preprocessRead(pileupread, file2reads, "f2")
              
              
              #Now examine preprocessed reads from each alignment
              # first get a list of all unique read names (not considering strand)
              file12readnames = set([])
              for read in file1reads.values():
                     file12readnames.add(read.alignment.query_name)
              for read in file2reads.values():
                     file12readnames.add(read.alignment.query_name)
              
              skippedReads = 0
              totalSkippedReadsRef = 0
              totalSkippedReadsNonref = 0
              passedReadsMarkedDup = 0
              totalDupReadsRef = 0
              totalDupReadsNonref = 0
              totalReads = 0
              totalReadsRef = 0
              totalReadsNonref = 0
              totalReadsInclFailed = 0 #All reads mapped to ref over SNP
              totalFailedReadsRef = 0
              totalFailedReadsNonref = 0
              totalFailedPairedReadsRef = 0
              totalFailedPairedReadsNonref = 0
              
              file1DupReadHash = dict()
              file2DupReadHash = dict()
              
              allelecounts = { 'A':0, 'C':0, 'G':0, 'T':0 }
              for readname in file12readnames:
                     file1pileupread1 = None
                     file2pileupread1 = None
                     file1pileupread2 = None
                     file2pileupread2 = None
                     
                     read_allele = None
                     
                     skipRead = False
                     
                     #Look for read 1 in a mate pair
                     file1pileupread1 = getPreprocessedRead(file1reads, readname + "True", "f1")
                     file2pileupread1 = getPreprocessedRead(file2reads, readname + "True", "f2")
                     read1validates = None
                     if file1pileupread1 != None or file2pileupread1 != None:
                            #Get the read allele from either read, later checked to make sure they match
                            read_allele = getSequenceAtSNV(file1pileupread1) or getSequenceAtSNV(file2pileupread1)
                            if preValidateIndividualPileupRead(file1pileupread1, args, curSNV)==False or preValidateIndividualPileupRead(file2pileupread1, args, curSNV)==False:
                                   skipRead = True
                                   
                                   #code to measure whether worth rescuing indel reads
#                                                        f1indel=None
#                                                        if file1pileupread1 != None:
#                                                               f1indel=re.search('[HSPDI]', file1pileupread1.alignment.cigarstring) is not None
#                                                        f2indel=None
#                                                        if file2pileupread1 != None:
#                                                               f2indel=re.search('[HSPDI]', file2pileupread1.alignment.cigarstring) is not None
#                                                        if f1indel==False and f2indel==True or f1indel==True and f2indel==False:
#                                                               print("could rescue indels: f1:", f1indel, "; f2:", f2indel, sep="", file=sys.stderr)

                            elif file1pileupread1 != None and file2pileupread1 != None:
                                   read1validates = validateDuallyAlignedPileupReads(file1pileupread1, file2pileupread1, prevSNV, curSNV, nextSNV, args)
                            else:
                                   read1validates = False
                     
                     #Look for read 2 in a mate pair
                     file1pileupread2 = getPreprocessedRead(file1reads, readname + "False", "f1")
                     file2pileupread2 = getPreprocessedRead(file2reads, readname + "False", "f2")
                     read2validates = None
                     if readname + "False" in file1reads or readname + "False" in file2reads:
                            #Get the read allele from either read, later checked to make sure they match
                            read_allele = getSequenceAtSNV(file1pileupread2) or getSequenceAtSNV(file2pileupread2)
                            if preValidateIndividualPileupRead(file1pileupread2, args, curSNV)==False or preValidateIndividualPileupRead(file2pileupread2, args, curSNV)==False:
                                   skipRead = True
                            elif readname + "False" in file1reads and readname + "False" in file2reads:
                                   read2validates = validateDuallyAlignedPileupReads(file1pileupread2, file2pileupread2, prevSNV, curSNV, nextSNV, args)
                            else:
                                   read2validates = False
                     
                     #Check if both reads in mate pair overlap SNP
                     read12validates = None
                     if read1validates != None and read2validates != None:
                            #f1 == f2 is checked above, so just check R1 == R2 for f1
                            if read1validates == True and read2validates == True and getSequenceAtSNV(file1pileupread1) == getSequenceAtSNV(file1pileupread2):
                                   read12validates = True
                            else:
                                   read12validates = False
                     
                     
                     file1isdup = None
                     file2isdup = None
                     #At this point all reads pass primary filters, so here we make sure they also pass the AI-specific ones
                     if not skipRead and (read1validates == True or read1validates == None) and (read2validates == True or read2validates == None) and (read12validates == True or read12validates == None):
                            #Check if these reads are a duplicate
                            file1isdup = checkDupRead(file1pileupread1, file1pileupread2, file1DupReadHash)
                            file2isdup = checkDupRead(file2pileupread1, file2pileupread2, file2DupReadHash)
                     
                     
                     if skipRead:
                            passRead = None
                            countRead = False
                            skippedReads += 1
                            
                            if read_allele == curSNV.mm10_ref_allele:
                                   totalSkippedReadsRef += 1
                            else:
                                   totalSkippedReadsNonref += 1
                            
                     elif (read1validates == True and read2validates == None or read1validates == None and read2validates == True or read1validates == True and read2validates == True and read12validates == True):
                            passRead = True
                            if (file1isdup == False and (file2isdup is None or file2isdup == False)):
                                   countRead = True
                                   totalReads += 1
                                   totalReadsInclFailed += 1
                                   allelecounts[read_allele] += 1
                                   
                                   #Count each PE read since we are interested mainly in sequencing/mapping bias
                                   #can use either f1 or f2 for cycle number
                                   if file1pileupread1 is not None:
                                          if read_allele == curSNV.mm10_ref_allele:
                                                 numRefReadsByCycle[getCycleNum(file1pileupread1) - 1] += 1
                                                 numRefReadsByReadLength[file1pileupread1.alignment.query_length - 1] += 1
                                          numReadsByCycle[getCycleNum(file1pileupread1) - 1] += 1
                                          numReadsByReadLength[file1pileupread1.alignment.query_length - 1] += 1
                                   if file1pileupread2 is not None:
                                          if read_allele == curSNV.mm10_ref_allele:
                                                 numRefReadsByCycle[getCycleNum(file1pileupread2) - 1] += 1
                                                 numRefReadsByReadLength[file1pileupread2.alignment.query_length - 1] += 1
                                          numReadsByCycle[getCycleNum(file1pileupread2) - 1] += 1
                                          numReadsByReadLength[file1pileupread2.alignment.query_length - 1] += 1
                                   
                                   if read_allele == curSNV.mm10_ref_allele:
                                          totalReadsRef += 1
                                   else:
                                          totalReadsNonref += 1

                            else:
                                   countRead = False
                                   #If the read failed only because it is a duplicate, don't include it in the total counts
                                   passedReadsMarkedDup += 1
                                   
                                   if read_allele == curSNV.mm10_ref_allele:
                                          totalDupReadsRef += 1
                                   else:
                                          totalDupReadsNonref += 1
                     else:
                            #Failed read
                            passRead = False
                            countRead = False
                            totalReadsInclFailed += 1
                            
                            if read_allele == curSNV.mm10_ref_allele:
                                   totalFailedReadsRef += 1
                            else:
                                   totalFailedReadsNonref += 1
                            
                            #Did both mates overlap SNV?
                            if read12validates is False:
                                   if read_allele == curSNV.mm10_ref_allele:
                                          totalFailedPairedReadsRef += 1
                                   else:
                                          totalFailedPairedReadsNonref += 1

                     
                     if args.verbose: print("\t\tSkip=", skipRead, ". Validation: R1=", read1validates, ", R2=", read2validates, ", R1R2=", read12validates, ", Pass all=", passRead, ". Duplicate: F1=", file1isdup, ", F2=", file2isdup, ". Count read=", countRead, sep="", file=sys.stderr)
              
              if args.verbose: print("\tCounts: ", allelecounts[curSNV.mm10_ref_allele], " ref/", allelecounts[curSNV.mm10_nonref_allele], " nonref. Skipped=", skippedReads, ", failed=", totalReadsInclFailed-totalReads, ", dups=", passedReadsMarkedDup, ", TotalReads=", totalReads, sep="", file=sys.stderr)
              
              if totalReads >= args.minReads:
                     #Print SNP
                     print(curSNV.mm10_chrom, curSNV.mm10_chromStart, curSNV.mm10_chromEnd, curSNV.rsid, totalReads, curSNV.mm10_ref_allele, allelecounts[curSNV.mm10_ref_allele], curSNV.mm10_nonref_allele, allelecounts[curSNV.mm10_nonref_allele], totalReadsInclFailed, skippedReads, curSNV.distToNearestSNP, sep="\t")
                     
                     #Only count reads if the SNV has enough reads to include
                     totalSNVsCounted += 1
                     grandTotalReadsSkipped += skippedReads
                     grandTotalReadsSkippedRef += totalSkippedReadsRef
                     grandTotalReadsSkippedNonref += totalSkippedReadsNonref
                     grandTotalPassedReadsMarkedDup += passedReadsMarkedDup
                     grandTotalDupReadsRef += totalDupReadsRef
                     grandTotalDupReadsNonref += totalDupReadsNonref
                     grandTotalReads += totalReads
                     grandTotalReadsRef += totalReadsRef
                     grandTotalReadsNonref += totalReadsNonref
                     grandTotalReadsInclFailed += totalReadsInclFailed
                     grandTotalFailedReadsRef += totalFailedReadsRef
                     grandTotalFailedReadsNonref += totalFailedReadsNonref
                     grandTotalFailedPairedReadsRef += totalFailedPairedReadsRef
                     grandTotalFailedPairedReadsNonref += totalFailedPairedReadsNonref
              else:
                     totalSNVsSkipped += 1


if __name__=='__main__':
       #Can add allow_abbrev=False in python 3.5+
       parser = argparse.ArgumentParser(prog = "countReads", description = "Counts reads for each allele", add_help=True)
       parser.add_argument('snpfilename', action='store', help='snpfilename')
       parser.add_argument('bamfile1name', action='store', help='bamfile1name')
       parser.add_argument('bamfile2name', action='store', help='bamfile2name')
       parser.add_argument('--maxReadLength', action='store', type=int, default=36, help='maxReadLength used to allocate space for recordkeeping [%(default)s]')
       parser.add_argument('--sample', action='store', default='', help='Optional sample label used for some output')
       parser.add_argument('--strain', action='store', default='', help='Optional strain label used for some output')
       parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
       parser.add_argument("--match-rsid", action='store_true', default=False, help = "Enforce that rsid in read name matches SNP (for mappability tracks)")
       parser.add_argument('--chrom', action='store', help='Jump to and process data for given <chromosome> only. NB applies to reference only (not non-reference, if different) [%(default)s]')
       parser.add_argument('--minReads', action='store', type=int, default=8, help='Minimum required reads passing all filters at SNP [%(default)s]')
       parser.add_argument('--minBaseQ', action='store', type=int, default=20, help='minBaseQ used to skip reads [%(default)s]')
       parser.add_argument('--cyclesToSkip', action='store', type=int, default=3, help='cyclesToSkip at 5\' end [%(default)s]')
       parser.add_argument('--maxTemplateLength', action='store', type=int, default=1000, help='Skip reads over maxTemplateLength bp [%(default)s]') #tried 150 before
       parser.add_argument('--filter', action='store', type=int, default=524, help='Bitmask of required flag [%(default)s]') #1548
       parser.add_argument('--require', action='store', type=int, default=524, help='Bitmask of flags used to fail reads [%(default)s]') #1548
       parser.add_argument('--permittedMismatches', action='store', type=int, default=1, help='In addition to one mismatch permitted when mapping to the genome for the other (ie ref vs nonref) allele [%(default)s]')
       parser.add_argument('--minMAPQ', action='store', type=int, default=30, help='minMAPQ to fail reads [%(default)s]')
       parser.add_argument('--version', action='version', version='%(prog)s 1.0')

       #argparse does not set an exit code upon this error
       #https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
       try:
              args = parser.parse_args()
       except argparse.ArgumentError as exc:
              print(exc.message, '\n', exc.argument)
              sys.exit(2)



       print("Counting alleles.\nParameters:", file=sys.stderr)
       print(args, file=sys.stderr)


       if(args.snpfilename =="-"):
              snpfile = sys.stdin
       else:
              snpfile = open(args.snpfilename, 'r') 
       bamfile1 = pysam.AlignmentFile(args.bamfile1name, "rb")
       bamfile2 = pysam.AlignmentFile(args.bamfile2name, "rb")

       lastChrom=""


       #Collect summary statistics
       totalChroms = 0
       totalSNVsSkipped = 0
       totalSNVsCounted = 0
       grandTotalReadsSkipped = 0
       grandTotalReadsSkippedRef = 0
       grandTotalReadsSkippedNonref = 0
       grandTotalPassedReadsMarkedDup = 0
       grandTotalDupReadsRef = 0
       grandTotalDupReadsNonref = 0
       grandTotalReads = 0
       grandTotalReadsRef = 0
       grandTotalReadsNonref = 0
       grandTotalReadsInclFailed = 0
       grandTotalFailedReadsRef = 0
       grandTotalFailedReadsNonref = 0
       grandTotalFailedPairedReadsRef = 0
       grandTotalFailedPairedReadsNonref = 0


       numRefReadsByCycle = [0] * args.maxReadLength
       numReadsByCycle = [0] * args.maxReadLength

       numRefReadsByReadLength = [0] * args.maxReadLength
       numReadsByReadLength = [0] * args.maxReadLength


       curSNV = None
       nextSNV = None

       try:
              with snpfile as source:
                     for line in source:
                            #print("read line:", line, file=sys.stderr)
                     
                            #Maintain pipeline of next/last SNV
                            #BUGBUG the first/last 2 SNVs on a chromosome won't be processed right
                            prevSNV= curSNV
                            curSNV = nextSNV
                     
                            needNextSNV = True
                            while(needNextSNV):
                                   nextSNV = SNV(line)
                                   needNextSNV = isUnwantedChromosome(nextSNV.mm10_chrom) or isUnwantedChromosome(nextSNV.denovo_chrom) or (args.chrom is not None and args.chrom != nextSNV.mm10_chrom)
                     
                            if curSNV is not None:
                                   processSNV(prevSNV, curSNV, nextSNV)
       
              #Do last SNV
              processSNV(curSNV, nextSNV, None)


       except KeyboardInterrupt:
              print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)


       finally:
              bamfile1.close()
              bamfile2.close()
              snpfile.close()
       
              print("\n\ncountReads.py statistics:", file=sys.stderr)
              print("#Sample", "strain", "totalChroms", "totalSNVsSkipped", "totalSNVsCounted", "grandTotalReadsSkipped", "grandTotalReadsSkippedRef", "grandTotalReadsSkippedNonref", "grandTotalPassedReadsMarkedDup", "grandTotalDupReadsRef", "grandTotalDupReadsNonref", "grandTotalReads", "grandTotalReadsRef", "grandTotalReadsNonref", "grandTotalReadsInclFailed", "grandTotalFailedReadsRef", "grandTotalFailedReadsNonref", "grandTotalFailedPairedReadsRef", "grandTotalFailedPairedReadsNonref", sep="\t", file=sys.stderr)
              print(args.sample, args.strain, totalChroms, totalSNVsSkipped, totalSNVsCounted, grandTotalReadsSkipped, grandTotalReadsSkippedRef, grandTotalReadsSkippedNonref, grandTotalPassedReadsMarkedDup, grandTotalDupReadsRef, grandTotalDupReadsNonref, grandTotalReads, grandTotalReadsRef, grandTotalReadsNonref, grandTotalReadsInclFailed, grandTotalFailedReadsRef, grandTotalFailedReadsNonref, grandTotalFailedPairedReadsRef, grandTotalFailedPairedReadsNonref, sep="\t", file=sys.stderr)
              print("pctRef by cycle: ", [None if numReadsByCycle[i]==0 else float(numRefReadsByCycle[i]) / float(numReadsByCycle[i]) for i in range(args.maxReadLength)], sep="", file=sys.stderr)
              print("pctRef by readLength: ", [None if numReadsByReadLength[i]==0 else float(numRefReadsByReadLength[i]) / float(numReadsByReadLength[i]) for i in range(args.maxReadLength)], sep="", file=sys.stderr)


       print("\nDone!", file=sys.stderr)
