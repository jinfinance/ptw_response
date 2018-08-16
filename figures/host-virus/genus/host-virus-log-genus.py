import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# parse the genus file
df = pd.read_csv('gFile.csv', header=0,names=["Virome_Sample","Metagenome_Sample","Phylum","Microbe","Phage_Abundance","Microbial_Abundance","Phage_Prevalence","Microbial_Prevalence","Phage_Host_Co_Detections"])

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
  sample_size=0
  for index, row in df.iterrows():
       if (row.iloc[3] == i):
         if (row.iloc[5] > 0) and (row.iloc[4] > 0): #both phage and host abundance > 0
            host.append(row.iloc[5])
            virus.append(row.iloc[4])
            sample_size+=1
     
  #log transformation
  loghost  = np.log10(host)
  logvirus = np.log10(virus)
  
  c=colors[num]
  #curve fitting with polyfit
  ax1.plot(loghost, np.poly1d(np.polyfit(loghost, logvirus, 1))(loghost), color=c) #best fit line
  #scatter plotting
  ax1.scatter(loghost, logvirus, s=5, label=i, color=c)

  num=num+1 
  print sample_size
    
plt.legend(bbox_to_anchor=(0.06, 1.05, 1., .102),labelspacing=0.5, handlelength=0.5, handletextpad=0.5,
           frameon=True, ncol=8, columnspacing=0.5, prop={'size': 4.2})


ax1.grid(False)
plt.xlim(-2,2)
plt.xlabel("log10(host abundance)")
plt.ylabel("log10(virus abundance)")
plt.show() 
fig.savefig('genus-host-virus-log-plot.png')


print 'Done..'
print 'no. Of all microbes = ', allMicrobs
print 'no. Of included microbes = ', len(microbe_list)

