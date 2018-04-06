#!/usr/bin/env python

#This code was initially writen by Fidel and Vivek in a python notebook where they were preparing the figures for their Nature paper. I just adapted the code and made it generalize to be able to use it on any upcoming HiC data.

import os
import sys
import argparse 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def parse_args():
  """

  """
  parser=argparse.ArgumentParser()
  #required argumnets:
  parser.add_argument("--Matrix",
                      "-m",
                      dest="TADmatrix",
                      type=str,
                      metavar="STR",
                      help="A matrix obtained from deeptools computeMatrix scale-regions with padded regions on both side. The input for the matrix is a TADboundery.bed from hicexplorer and an active_histonmark.bw",
                      required=True)

  parser.add_argument("--output",
                      "-o",
                      dest="output",
                      type=str,
                      metavar="STR",
                      help="output figure",
                      required=True)

  parser.add_argument("--figsize",
                      "-fs",
                      dest="FigureSize",
                      type= int,
                      metavar="INT",
                      help="figure size",
                      default = 3.5)

  parser.add_argument("--regionLength",
                      "-rl",
                      dest="regionLength",
                      type= str,
                      metavar="STR",
                      help="length of the region",
                      default = "1kb")

  return parser


def main():
   """
   Main function to visualise TAD correlation
   """
   parser = parse_args()
   args = parser.parse_args()
   table = pd.read_csv(args.TADmatrix,  sep='\t', skiprows=1, header=None)   
   values = table.iloc[:, 6:]
   values = np.array(values)   
   values[np.isnan(values)]=0

   matrix_size = values.shape[1]
   correlation_table = np.zeros((matrix_size, matrix_size))
   for col1 in range(matrix_size):
       for col2 in range(matrix_size):
           if col1 >= col2:
              correlation_table[col1,col2] = spearmanr(values[:,col1],values[:,col2])[0]
   final_correlation_table = correlation_table + correlation_table.T - np.diag(correlation_table.diagonal())
   print(final_correlation_table)   
   heatmap = plt.figure(figsize=(args.FigureSize,args.FigureSize))
   ax = heatmap.add_subplot(111)
   cax = ax.imshow(final_correlation_table, interpolation='nearest', vmin=0, vmax=1)
   cbar = plt.colorbar(cax, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1], fraction=0.046, pad=0.04)
   plt.xticks([0,14,29, 44], ['-'+args.regionLength,'TAD start','TAD end', args.regionLength])
   plt.yticks([0,14,29, 44], ['-'+args.regionLength,'TAD start','TAD end', args.regionLength])
   plt.tight_layout()
   plt.savefig(args.output)


if __name__ == "__main__":
    main()

