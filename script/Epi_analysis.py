#!/usr/bin/python
import os
import sys

def reading_Dresist_file(filename):
  infile = open(filename,'r')
  Dresist_dict = {}
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    Dmut = line[0]
    fit  = line[4]
    ipt_count     = line[3]
    resist_CR9114 = line[10]
    resist_FI6v3  = line[12]
    if float(ipt_count) < 20: continue
    if float(fit) < 0.5: continue
    Dresist_dict[Dmut] = {'CR9114': resist_CR9114, 'FI6v3': resist_FI6v3}
  infile.close()
  return Dresist_dict

def reading_Sresist_file(filename):
  infile = open(filename, 'r')
  Sresist_dict = {}
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = line[4]
    ipt_count     = line[3]
    resist_CR9114 = line[10]
    resist_FI6v3  = line[12]
    if float(ipt_count) < 20: continue
    if float(fit) < 0.5: continue
    Sresist_dict[mut] = {'CR9114': resist_CR9114, 'FI6v3': resist_FI6v3}
  infile.close()
  return Sresist_dict

def compute_expect_fit(outfile, Dresist_dict, Sresist_dict):
  print "writing: %s" % outfile
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['Dmut', 'resist_CR9114', 'expect_resist_CR9114', 'resist_FI6v3', 'expect_resist_FI6v3'])+"\n")
  for Dmut in Dresist_dict.keys():
    mut1 = Dmut.rsplit('/')[0]
    mut2 = Dmut.rsplit('/')[1]
    if mut1 in Sresist_dict.keys() and mut2 in Sresist_dict.keys():
      resist_CR9114 = Dresist_dict[Dmut]['CR9114']
      resist_FI6v3  = Dresist_dict[Dmut]['FI6v3']
      expect_resist_CR9114 = float(Sresist_dict[mut1]['CR9114'])*float(Sresist_dict[mut2]['CR9114'])
      expect_resist_FI6v3 = float(Sresist_dict[mut1]['FI6v3'])*float(Sresist_dict[mut2]['FI6v3'])
      resist_CR9114 = Dresist_dict[Dmut]['CR9114']
      outfile.write("\t".join(map(str,[Dmut, resist_CR9114, expect_resist_CR9114, resist_FI6v3, expect_resist_FI6v3]))+"\n")
  outfile.close()
    
def main():
  outfile      = "result/Resist_epi.tsv"
  Dresist_file = "result/Resist_D.tsv"
  Sresist_file = "result/Resist_S.tsv"
  Dresist_dict = reading_Dresist_file(Dresist_file)
  Sresist_dict = reading_Sresist_file(Sresist_file)
  compute_expect_fit(outfile, Dresist_dict, Sresist_dict)

if __name__ == "__main__":
  main()
