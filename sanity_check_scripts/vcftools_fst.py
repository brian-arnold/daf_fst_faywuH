#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt


def main():

  d = '/scratch/gpfs/bjarnold/chernobylWolves/analyses/Fst'
  # fst_file = "CEZ_v_BLR.windowed.weir.fst"
  fst_file = "YLS_v_BLR.windowed.weir.fst"
  positions = []
  fsts = []
  with open(f'{d}/{fst_file}', 'r') as f:
    for line in f:
      if line.startswith("chrX"):
        line = line.strip().split()
        start = int(line[1])
        end = int(line[2])
        mid = ((start + end) / 2)
        positions.append(mid)
        fsts.append(float(line[5]))
  scat = sns.scatterplot(x=positions, y=fsts, linewidth=0, alpha=0.1)
  plt.title("YLS_v_BLR")  
  sns.despine()
  scat.get_figure().savefig("vcftools_fst_vs_position.png", dpi=300)


if __name__ == '__main__':
  main()
