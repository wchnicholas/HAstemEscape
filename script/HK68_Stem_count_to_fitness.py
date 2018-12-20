#!/usr/bin/python
import os
import sys
import glob

def compile_fitness(Input_dict, R1A0_dict, R2A0_dict, R1A1_dict, R2A1_dict, R1A5_dict, R2A5_dict, 
                                                      R1F1_dict, R2F1_dict, R1F5_dict, R2F5_dict, mut_count_cutoff, outfile):
  print "writing %s" % outfile
  outfile = open(outfile,'w')
  outfile.write("\t".join(['mut','Input_count','rep1_no_antibody_fit','rep2_no_antibody_fit',
                           'rep1_2ug_CR9114_fit','rep2_2ug_CR9114_fit','rep1_10ug_CR9114_fit','rep2_10ug_CR9114_fit',
                           'rep1_300ng_FI6v3_fit','rep2_300ng_FI6v3_fit','rep1_2500ng_FI6v3_fit','rep2_2500ng_FI6v3_fit'])+"\n")
  muts_dict = {}
  muts = []
  [muts.extend(mut_dict.keys()) for mut_dict in [Input_dict, R1A0_dict, R2A0_dict, R1A1_dict, R2A1_dict, R1A5_dict, R2A5_dict,
                                                                                   R1F1_dict, R2F1_dict, R1F5_dict, R2F5_dict]]
  for mut in muts:
    if mut.count('/') >= mut_count_cutoff: continue
    Input_count = float(Input_dict[mut]) + float(1) if Input_dict.has_key(mut) else float(1)
    R1A0_count  = float(R1A0_dict[mut]) + float(1)  if R1A0_dict.has_key(mut)  else float(1)
    R2A0_count  = float(R2A0_dict[mut]) + float(1)  if R2A0_dict.has_key(mut)  else float(1)
    R1A1_count  = float(R1A1_dict[mut]) + float(1)  if R1A1_dict.has_key(mut)  else float(1)
    R2A1_count  = float(R2A1_dict[mut]) + float(1)  if R2A1_dict.has_key(mut)  else float(1)
    R1A5_count  = float(R1A5_dict[mut]) + float(1)  if R1A5_dict.has_key(mut)  else float(1)
    R2A5_count  = float(R2A5_dict[mut]) + float(1)  if R2A5_dict.has_key(mut)  else float(1)
    R1F1_count  = float(R1F1_dict[mut]) + float(1)  if R1F1_dict.has_key(mut)  else float(1)
    R2F1_count  = float(R2F1_dict[mut]) + float(1)  if R2F1_dict.has_key(mut)  else float(1)
    R1F5_count  = float(R1F5_dict[mut]) + float(1)  if R1F5_dict.has_key(mut)  else float(1)
    R2F5_count  = float(R2F5_dict[mut]) + float(1)  if R2F5_dict.has_key(mut)  else float(1)
    muts_dict[mut]  = {'Input_count': Input_count,
                       'R1A0_ratio': R1A0_count/Input_count,
                       'R2A0_ratio': R2A0_count/Input_count,
                       'R1A1_ratio': R1A1_count/Input_count,
                       'R2A1_ratio': R2A1_count/Input_count,
                       'R1A5_ratio': R1A5_count/Input_count,
                       'R2A5_ratio': R2A5_count/Input_count,
                       'R1F1_ratio': R1F1_count/Input_count,
                       'R2F1_ratio': R2F1_count/Input_count,
                       'R1F5_ratio': R1F5_count/Input_count,
                       'R2F5_ratio': R2F5_count/Input_count}
  for mut in muts_dict.keys():
    Input_count = muts_dict[mut]['Input_count']-1
    R1A0_fit    = muts_dict[mut]['R1A0_ratio']/muts_dict['WT']['R1A0_ratio']
    R2A0_fit    = muts_dict[mut]['R2A0_ratio']/muts_dict['WT']['R2A0_ratio']
    R1A1_fit    = muts_dict[mut]['R1A1_ratio']/muts_dict['WT']['R1A1_ratio']
    R2A1_fit    = muts_dict[mut]['R2A1_ratio']/muts_dict['WT']['R2A1_ratio']
    R1A5_fit    = muts_dict[mut]['R1A5_ratio']/muts_dict['WT']['R1A5_ratio']
    R2A5_fit    = muts_dict[mut]['R2A5_ratio']/muts_dict['WT']['R2A5_ratio']
    R1F1_fit    = muts_dict[mut]['R1F1_ratio']/muts_dict['WT']['R1F1_ratio']
    R2F1_fit    = muts_dict[mut]['R2F1_ratio']/muts_dict['WT']['R2F1_ratio']
    R1F5_fit    = muts_dict[mut]['R1F5_ratio']/muts_dict['WT']['R1F5_ratio']
    R2F5_fit    = muts_dict[mut]['R2F5_ratio']/muts_dict['WT']['R2F5_ratio']
    outfile.write("\t".join(map(str,[mut, int(Input_count), R1A0_fit, R2A0_fit, R1A1_fit, R2A1_fit, R1A5_fit, R2A5_fit,
                                                                                R1F1_fit, R2F1_fit, R1F5_fit, R2F5_fit]))+"\n")
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
  S_outfile   = 'result/Fitness_S.tsv'
  D_outfile   = 'result/Fitness_D.tsv'
  single_input_dict       =  reading_countfile('count/count_single_input.tsv')
  single_no_antibody_dict = reading_countfile('count/count_single_rep1_no_antibody.tsv')
  single_rep1_2ug_CR9114_dict    = reading_countfile('count/count_single_rep1_2ug_CR9114.tsv')
  single_rep1_10ug_CR9114_dict   = reading_countfile('count/count_single_rep1_10ug_CR9114.tsv')
  single_rep2_no_antibody_dict   = reading_countfile('count/count_single_rep2_no_antibody.tsv')
  single_rep2_2ug_CR9114_dict    = reading_countfile('count/count_single_rep2_2ug_CR9114.tsv')
  single_rep2_10ug_CR9114_dict   = reading_countfile('count/count_single_rep2_10ug_CR9114.tsv')
  single_rep1_300ng_FI6v3_dict   = reading_countfile('count/count_single_rep1_300ng_FI6v3.tsv')
  single_rep1_2500ng_FI6v3_dict  = reading_countfile('count/count_single_rep1_2500ng_FI6v3.tsv')
  single_rep2_300ng_FI6v3_dict   = reading_countfile('count/count_single_rep2_300ng_FI6v3.tsv')
  single_rep2_2500ng_FI6v3_dict  = reading_countfile('count/count_single_rep2_2500ng_FI6v3.tsv')
  double_input_dict       = reading_countfile('count/count_double_input.tsv')
  double_rep1_no_antibody_dict   = reading_countfile('count/count_double_rep1_no_antibody.tsv')
  double_rep1_2ug_CR9114_dict    = reading_countfile('count/count_double_rep1_2ug_CR9114.tsv')
  double_rep1_10ug_CR9114_dict   = reading_countfile('count/count_double_rep1_10ug_CR9114.tsv')
  double_rep2_no_antibody_dict   = reading_countfile('count/count_double_rep2_no_antibody.tsv')
  double_rep2_2ug_CR9114_dict    = reading_countfile('count/count_double_rep2_2ug_CR9114.tsv')
  double_rep2_10ug_CR9114_dict   = reading_countfile('count/count_double_rep2_10ug_CR9114.tsv')
  double_rep1_300ng_FI6v3_dict   = reading_countfile('count/count_double_rep1_300ng_FI6v3.tsv')
  double_rep1_2500ng_FI6v3_dict  = reading_countfile('count/count_double_rep1_2500ng_FI6v3.tsv')
  double_rep2_300ng_FI6v3_dict   = reading_countfile('count/count_double_rep2_300ng_FI6v3.tsv')
  double_rep2_2500ng_FI6v3_dict  = reading_countfile('count/count_double_rep2_2500ng_FI6v3.tsv')
  compile_fitness(single_input_dict, single_no_antibody_dict, single_rep2_no_antibody_dict,
                  single_rep1_2ug_CR9114_dict, single_rep2_2ug_CR9114_dict,
                  single_rep1_10ug_CR9114_dict, single_rep2_10ug_CR9114_dict,
                  single_rep1_300ng_FI6v3_dict, single_rep2_300ng_FI6v3_dict,
                  single_rep1_2500ng_FI6v3_dict, single_rep2_2500ng_FI6v3_dict, 1, S_outfile)
  compile_fitness(double_input_dict, double_rep1_no_antibody_dict, double_rep2_no_antibody_dict,
                  double_rep1_2ug_CR9114_dict, double_rep2_2ug_CR9114_dict,
                  double_rep1_10ug_CR9114_dict, double_rep2_10ug_CR9114_dict,
                  double_rep1_300ng_FI6v3_dict, double_rep2_300ng_FI6v3_dict,
                  double_rep1_2500ng_FI6v3_dict, double_rep2_2500ng_FI6v3_dict, 2, D_outfile)

if __name__ == "__main__":
  main()
