#!/usr/bin/python
import os
import sys
import glob

def compile_fitness(Input_dict, R1A0_dict, R2A0_dict, R1F1_dict, R2F1_dict, mut_count_cutoff, outfile):
  print "writing %s" % outfile
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','Input_count',
                           'rep1_no_antibody_fit','rep2_no_antibody_fit',
                           'rep1_300ng_FI6v3_fit','rep2_300ng_FI6v3_fit'])+"\n")
  muts_dict = {}
  muts = []
  [muts.extend(mut_dict.keys()) for mut_dict in [Input_dict, R1A0_dict, R2A0_dict, R1F1_dict, R2F1_dict]]
  for mut in muts:
    if mut.count('/') >= mut_count_cutoff: continue
    Input_count = float(Input_dict[mut]) + float(1) if Input_dict.has_key(mut) else float(1)
    R1A0_count  = float(R1A0_dict[mut]) + float(1)  if R1A0_dict.has_key(mut)  else float(1)
    R2A0_count  = float(R2A0_dict[mut]) + float(1)  if R2A0_dict.has_key(mut)  else float(1)
    R1F1_count  = float(R1F1_dict[mut]) + float(1)  if R1F1_dict.has_key(mut)  else float(1)
    R2F1_count  = float(R2F1_dict[mut]) + float(1)  if R2F1_dict.has_key(mut)  else float(1)
    muts_dict[mut]  = {'Input_count': Input_count,
                       'R1A0_ratio': R1A0_count/Input_count,
                       'R2A0_ratio': R2A0_count/Input_count,
                       'R1F1_ratio': R1F1_count/Input_count,
                       'R2F1_ratio': R2F1_count/Input_count}
  for mut in muts_dict.keys():
    Input_count = muts_dict[mut]['Input_count']-1
    R1A0_fit    = muts_dict[mut]['R1A0_ratio']/muts_dict['WT']['R1A0_ratio']
    R2A0_fit    = muts_dict[mut]['R2A0_ratio']/muts_dict['WT']['R2A0_ratio']
    R1F1_fit    = muts_dict[mut]['R1F1_ratio']/muts_dict['WT']['R1F1_ratio']
    R2F1_fit    = muts_dict[mut]['R2F1_ratio']/muts_dict['WT']['R2F1_ratio']
    outfile.write("\t".join(map(str,[mut, int(Input_count), R1A0_fit, R2A0_fit, R1F1_fit, R2F1_fit]))+"\n")
  outfile.close()

def reading_countfile(countfile):
  mut_dict = {}
  infile = open(countfile,'r')
  for line in infile.xreadlines():
    if 'mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut   = line[0]
    count = line[1]
    mut_dict[mut] = float(count)
  infile.close()
  return mut_dict

def main():
  SI06_outfile   = 'result/Fitness_SI06.tsv'
  Mich15_outfile   = 'result/Fitness_Mich15.tsv'
  SI06_input_dict       =  reading_countfile('count/count_SI06_input.tsv')
  SI06_rep1_no_antibody_dict   = reading_countfile('count/count_SI06_rep1_no_antibody.tsv')
  SI06_rep2_no_antibody_dict   = reading_countfile('count/count_SI06_rep2_no_antibody.tsv')
  SI06_rep1_300ng_FI6v3_dict   = reading_countfile('count/count_SI06_rep1_300ng_FI6v3.tsv')
  SI06_rep2_300ng_FI6v3_dict   = reading_countfile('count/count_SI06_rep2_300ng_FI6v3.tsv')
  Mich15_input_dict       = reading_countfile('count/count_Mich15_input.tsv')
  Mich15_rep1_no_antibody_dict   = reading_countfile('count/count_Mich15_rep1_no_antibody.tsv')
  Mich15_rep2_no_antibody_dict   = reading_countfile('count/count_Mich15_rep2_no_antibody.tsv')
  Mich15_rep1_300ng_FI6v3_dict   = reading_countfile('count/count_Mich15_rep1_300ng_FI6v3.tsv')
  Mich15_rep2_300ng_FI6v3_dict   = reading_countfile('count/count_Mich15_rep2_300ng_FI6v3.tsv')
  compile_fitness(SI06_input_dict, SI06_rep1_no_antibody_dict, SI06_rep2_no_antibody_dict,
                  SI06_rep1_300ng_FI6v3_dict, SI06_rep2_300ng_FI6v3_dict, 1, SI06_outfile)
  compile_fitness(Mich15_input_dict, Mich15_rep1_no_antibody_dict, Mich15_rep2_no_antibody_dict,
                  Mich15_rep1_300ng_FI6v3_dict, Mich15_rep2_300ng_FI6v3_dict, 1, Mich15_outfile)

if __name__ == "__main__":
  main()
