import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# parse the phylum file
df = pd.read_csv('phFile.csv', header=0,names=["Virome_Sample","Metagenome_Sample","Phylum","Microbe","Phage_Abundance","Microbial_Abundance","Phage_Prevalence","Microbial_Prevalence","Phage_Host_Co_Detections"])

#get the complete microbe list
microbe_list = df.iloc[:,3]
microbe_list = set(microbe_list) #to remove duplicates
allMicrobs = len(microbe_list)

print microbe_list 
print len(microbe_list)


#filter the microbes
#if 'Microbial_Prevalence' is below 5, drop that microbe from the list
for i in microbe_list.copy():
   for index, row in df.iterrows():
     if (row.iloc[3] == i) and (row.iloc[7] < 5):
       microbe_list.remove(i)
       break
      

#set the figure size
fig = plt.figure(figsize=(6,6), dpi=200)
ax1 = fig.add_subplot(111)

number = len(microbe_list)
colors = plt.cm.tab20b(np.linspace(0, 1, number))
num=0   #to change the color of each scatter plot



for i in microbe_list:
  host=list()
  virus=list()
  vhr=list()
  for index, row in df.iterrows():
       if (row.iloc[3] == i):
         if (row.iloc[5] > 0) and (row.iloc[4] > 0): #both phage and host abundance > 0
            host.append(row.iloc[5])
            virus.append(row.iloc[4])
            vhr.append(row.iloc[4]/row.iloc[5]) #vhr = virus abundance / host abundance

  #log transformation
  loghost  = np.log10(host)
  logvirus = np.log10(virus)
  logvhr   = np.log10(vhr) 
  
  c=colors[num]
  #curve fitting with polyfit
  ax1.plot(loghost, np.poly1d(np.polyfit(loghost, logvhr, 1))(loghost), color=c) #best fit line
  ax1.scatter(loghost, logvhr, s=5, label=i, color=c)  
  num=num+1 
  
plt.legend(bbox_to_anchor=(-0.01, 1.01, 1., .102),labelspacing=0.5, handlelength=0.5, handletextpad=0.5,
           frameon=True, ncol=7, columnspacing=0.5, prop={'size': 5})

#add a -1 slope base line, and annotation
plt.plot([-2,0], [-5,-7], 'k--')
plt.annotate('-1 slope', xy=(-1.2, -5.5), xytext=(-1.5, -5.9),rotation=-32,
            )

ax1.grid(False)
plt.xlim(-2,2)
plt.ylim(-7,-1)
plt.xlabel("log10(host abundance)")
plt.ylabel("log10(VHR)")
plt.show() 
fig.savefig('phylum-host-vhr-log-plot.png')


print 'Done..'
print 'no. Of all microbes = ', allMicrobs
print 'no. Of included microbes = ', len(microbe_list)

