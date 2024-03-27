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
  fst_file = 'combined_windowed_fst.txt'
  positions = []
  # fsts = []
  fsts = defaultdict(list)
  with open(f'{d}/{fst_file}', 'r') as f:
    for line in f:
      if line.startswith("chrX"):
        line = line.strip().split()
        start = int(line[1])
        end = int(line[2])
        mid = ((start + end) / 2)
        positions.append(mid)
        # fsts.append(float(line[5]))
        fsts['CEZ_YLS'].append(float(line[3]))
        fsts['CEZ_BLR'].append(float(line[4]))
        fsts['YLS_BLR'].append(float(line[5]))

  fig, axs = plt.subplots(3,1, figsize=(10, 10))
  scat = sns.scatterplot(x=positions, y=fsts['CEZ_YLS'], linewidth=0, alpha=0.1, ax=axs[0])
  axs[0].set_title("CEZ_YLS")
  sns.scatterplot(x=positions, y=fsts['CEZ_BLR'], linewidth=0, alpha=0.1, ax=axs[1])
  axs[1].set_title("CEZ_BLR")
  sns.scatterplot(x=positions, y=fsts['YLS_BLR'], linewidth=0, alpha=0.1, ax=axs[2])
  axs[2].set_title("YLS_BLR")

  sns.despine()
  # save figure
  scat.get_figure().savefig("vcftools_fst_vs_position_3pop.png", dpi=300)

if __name__ == '__main__':
  main()
