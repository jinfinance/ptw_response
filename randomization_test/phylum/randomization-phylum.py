'''
2017/2018

This code performs a permutation test to evaluate the significance of the measured correlation
between virus and host abundances in the phylum datasets. For each host dataset, we generated 10^4
randomized samples of size n_i (the number of data points for virus-host pair i).
In each sample, the virus abundances are permuted without replacement while holding the host
abundances constant. Then the distributions of measured correlations are compared to the original
samples. The results are displayed in a violin diagram form.

'''


import numpy as np
import pandas as pd
import seaborn as sns
import csv
import matplotlib.pyplot as plt
from scipy import stats
from random import randint
from mpmath import mp
import matplotlib
import matplotlib.gridspec as gridspec
from collections import OrderedDict
import copy
import math
from numpy import median

df = pd.read_csv('phFile.csv', header=0,names=["Virome_Sample","Metagenome_Sample","Phylum","Microbe","Phage_Abundance","Microbial_Abundance","Phage_Prevalence","Microbial_Prevalence","Phage_Host_Co_Detections"])

microbes = df.iloc[:,3]
microbe_list = list(OrderedDict.fromkeys(microbes))
allMicrobs = len(microbe_list)
microbe_list2 = copy.copy(microbe_list)

#filter the microbes
for i in microbe_list2:
   for index, row in df.iterrows():
     if (row.iloc[3] == i) and (row.iloc[7] < 5):
       microbe_list.remove(i)
       break


f = open('phylumSlope.txt','w')
f2 = open('phylumSpearsmanRank.txt','w')
gs = gridspec.GridSpec(len(microbe_list),1)
fig = plt.figure()
fig.set_size_inches(20, 7, forward=True)
num=-1 #To change the position of each violin plot on the y-axis

for i in microbe_list:

  print 'Processing host: ', i
  host=list()
  virus=list()
  stat_slope=list()
  stat_rho=list()
  num+=1
 
  sample_size=0
  for index, row in df.iterrows():
       if (row.iloc[3] == i):
         if (row.iloc[5] > 0) and (row.iloc[4] > 0):
              host.append(row.iloc[5])
              virus.append(row.iloc[4])
              sample_size+=1


  loghost  = np.log10(host)
  logvirus = np.log10(virus)
  slopeO, interceptO, r_valueO, p_valueO, std_errO = stats.linregress(loghost,logvirus) 
  rhoO, pvalO = stats.spearmanr(loghost,logvirus)
     
  #---randomization---
  for j in range(0,10000):

      bs_loghost  = loghost
      bs_logvirus = np.random.choice(logvirus,len(logvirus),replace=False)

      slope, intercept, r_value, p_value, std_err = stats.linregress(bs_loghost,bs_logvirus) 
      stat_slope.append(slope)

      rho, pval = stats.spearmanr(bs_loghost,bs_logvirus)
      stat_rho.append(rho)
      
  #plot  
  st=i+ '\n' + ' (' + str(sample_size) + ')'
 
  if(num==0):
        ax = fig.add_subplot(gs[num])
  else:
        ax = fig.add_subplot(gs[num],sharex=ax)
 
  sns.boxplot(stat_slope, showfliers=False, showbox=False, whis=[0.5,99.5], color='gray', linewidth=1.75)

  if( (slopeO < np.percentile(stat_slope, [0.5])) or (slopeO > np.percentile(stat_slope, [99.5]))):
      ax.set_ylabel(st, size =6, rotation=0, color='gray', horizontalalignment='right')
      sns.violinplot(stat_slope, color='mistyrose', linewidth=0.85, inner="box")

  else:
      ax.set_ylabel(st, size =6, rotation=0, horizontalalignment='right')
      sns.violinplot(stat_slope, color='lightgray', linewidth=0.85, inner="box")

  plt.axvline(x=slopeO, color='royalblue', linestyle='--', linewidth=2)
  
  ax.get_yaxis().set_label_coords(-0.01,-0.03)
  ax.set_xlim([-1.25,1.25])
  ax.yaxis.set_ticks_position('none')
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  if(num!=len(microbe_list)-1):
       plt.tick_params(axis='x',labelbottom='off', bottom='off')
       ax.spines['bottom'].set_visible(False)
  if(num==len(microbe_list)-1):
       ax.set_xlabel('Slope', size =7)
 
    
  print 'sample size for this host = ', sample_size
  print '----------------------------------------------------------'

  #write to txt files
  line = i + ' & ' + str(sample_size) + ' & ' + str("%.4f" % slopeO) + ' & ' + str("%.4f" % np.median(stat_slope)) + ' & ' + str("%.4f" % np.percentile(stat_slope, [2.5])) + ' & ' + str("%.4f" % np.percentile(stat_slope, [97.5])) + ' & ' + str("%.4f" % np.percentile(stat_slope, [0.5])) + ' & ' + str("%.4f" % np.percentile(stat_slope, [99.5])) + '\\' + '\\' +  '\n'
  f.write(line)
  line = i + ' & ' + str(sample_size) + ' & ' + str("%.4f" % rhoO) + ' & ' + str("%.4f" % np.median(stat_rho)) + ' & ' + str("%.4f" % np.percentile(stat_rho, [2.5])) + ' & ' + str("%.4f" % np.percentile(stat_rho, [97.5])) + ' & ' + str("%.4f" % np.percentile(stat_rho, [0.5])) + ' & ' + str("%.4f" % np.percentile(stat_rho, [99.5])) + '\\' + '\\' + '\n'
  f2.write(line)
  
f.close()
f2.close()
plt.show()
#fig.savefig('phylum-violin-plot-n=1000.pdf')
print 'Done..'
