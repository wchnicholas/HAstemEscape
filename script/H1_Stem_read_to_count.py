#!/usr/bin/python
import os
import sys
import glob
import string
from Bio import SeqIO
from collections import Counter

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def readingFile2SampleID(file2SampleID_file):
  file2SampleID_dict = {}
  infile = open(file2SampleID_file,'r')
  for line in infile.xreadlines():
    if 'R1File' in line: continue
    line = line.rstrip().rsplit("\t")
    file2SampleID_dict[line[0]] = line[1]
  infile.close()
  return file2SampleID_dict

def readingBarcode2Resi(Barcode2Resi_file):
  Barcode2Resi_dict = {}
  infile = open(Barcode2Resi_file,'r')
  for line in infile.xreadlines():
    if 'Barcode' in line: continue
    line = line.rstrip().rsplit("\t")
    Barcode2Resi_dict[line[0]] = line[1]
  infile.close()
  return Barcode2Resi_dict

def readingWTcodon(WTcodon_file):
  WTcodon_dict = {}
  infile = open(WTcodon_file,'r')
  for line in infile.xreadlines():
    if 'Codon' in line: continue
    line = line.rstrip().rsplit("\t")
    WTcodon_dict[line[0]] = line[1]
  infile.close()
  return WTcodon_dict

def Processlib(R1file, WTpep, residues):
  R2file = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file, "fastq")
  R2records = SeqIO.parse(R2file,"fastq") 
  muts = []
  count_record = 0
  BC_count_mismatch = 0
  BC_count_good = 0
  BC_count_bad  = 0
  for R1record in R1records:
    count_record += 1
    #if count_record == 500: break
    R2record = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    Mutcodon_dict = {
      'R1resi42': str(R1seq[25:28]),   'R2resi42': str(rc(R2seq[235:238])),
      'R1resi45': str(R1seq[34:37]),   'R2resi45': str(rc(R2seq[226:229])),
      'R1resi46': str(R1seq[37:40]),   'R2resi46': str(rc(R2seq[223:226])),
      'R1resi47': str(R1seq[40:43]),   'R2resi47': str(rc(R2seq[220:223])),
      'R1resi48': str(R1seq[43:46]),   'R2resi48': str(rc(R2seq[217:220])),
      'R1resi49': str(R1seq[46:49]),   'R2resi49': str(rc(R2seq[214:217])),
      'R1resi52': str(R1seq[55:58]),   'R2resi52': str(rc(R2seq[205:208])),
      'R1resi111':str(R1seq[232:235]), 'R2resi111':str(rc(R2seq[28:31]))
      }
    R1roi = Mutcodon_dict['R1resi42']+Mutcodon_dict['R1resi45']+Mutcodon_dict['R1resi46']+Mutcodon_dict['R1resi47']+ \
            Mutcodon_dict['R1resi48']+Mutcodon_dict['R1resi49']+Mutcodon_dict['R1resi52']+Mutcodon_dict['R1resi111']
    R2roi = Mutcodon_dict['R2resi42']+Mutcodon_dict['R2resi45']+Mutcodon_dict['R2resi46']+Mutcodon_dict['R2resi47']+ \
            Mutcodon_dict['R2resi48']+Mutcodon_dict['R2resi49']+Mutcodon_dict['R2resi52']+Mutcodon_dict['R2resi111']
    if 'N' in R1roi or 'N' in R2roi: continue
    R1pep = translation(R1roi)
    R2pep = translation(R2roi)
    if R1pep == R2pep:
      mut = []
      for wtaa, mutaa, resi in zip(WTpep, R1pep, residues):
        if wtaa!=mutaa:
          mut.append(wtaa+resi+mutaa)
      if len(mut)==0:   mut = 'WT'
      elif len(mut)==1: mut = mut[0]
      else: continue
      muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def Output(mutcount_dict, outfile):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','count'])+"\n")
  for mut in mutcount_dict.keys():
    outfile.write("\t".join(map(str,[mut,mutcount_dict[mut]]))+"\n")
  outfile.close()

def main():
  WTseq_dict = {'Mich15': 'QIDKITVH', 'SI06': 'QINGITVH'}
  residues   = ['42','45','46','47','48','49','52','111']
  R1filenames = glob.glob('fastq/*R1*.fastq')
  file2SampleID_dict = readingFile2SampleID('doc/SampleID.tsv')
  WTcodon_dict       = readingWTcodon('doc/WTcodon.tsv')
  print "Analyzing %i samples" % len(R1filenames)
  for R1filename in sorted(R1filenames,key=lambda x:int(x.rsplit('_S')[1].rsplit('_')[0])):
    SampleID = file2SampleID_dict[R1filename]
    strain   = SampleID.rsplit('_')[0]
    WTpep    = WTseq_dict[strain]
    outfile  = 'count/count_'+SampleID+'.tsv'
    print "Analyzing %s (%s)" % (R1filename, SampleID)
    mutcount_dict = Processlib(R1filename, WTpep, residues)
    Output(mutcount_dict, outfile)

if __name__ == "__main__":
  main()
