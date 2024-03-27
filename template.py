#!/usr/bin/python -tt

"""
"""

import re
import sys
import os
from collections import defaultdict
import pandas as pd
import numpy as np
import seaborn as sns

def main():

  df_21 = pd.read_csv('CEZ_BLR_YLS_WIN21_OFF1_MINAF1_FayWuH.txt', sep='\t')
  df_41 = pd.read_csv('CEZ_BLR_YLS_WIN41_OFF1_MINAF1_FayWuH.txt', sep='\t')

  print(np.var(df_21.FayWuH))
  print(np.var(df_41.FayWuH))

  #print(np.mean(abs(df.midSNPpos - df.midPos)))
  #print(np.median(abs(df.midSNPpos - df.midPos)))
  #g = sns.histplot(abs(df.midSNPpos - df.midPos))
  #g = sns.histplot(df.maxPos - df.minPos)
  # save seaborn figure to file
  #g.figure.savefig('test.png')
  

if __name__ == '__main__':
  main()
