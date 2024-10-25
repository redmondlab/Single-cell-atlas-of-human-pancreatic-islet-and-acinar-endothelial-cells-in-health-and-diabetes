import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import numpy as np

#PATH TO FILTERED DATAFRAME 
PATH_DF = ""

#Which folder is diabetic either 'Folder1' or 'Folder2'
DIABETIC_FOLDER = ""
NONDIABETIC_FOLDER = ""

#MARKER BEING TESTED
MARKER = ""

df = pd.read_csv(PATH_DF)


df['Condition'] = ''
df.loc[df['Folder'] == DIABETIC_FOLDER,'Condition'] = 'Diabetic'
df.loc[df['Folder'] == NONDIABETIC_FOLDER,'Condition'] = 'Non-Diabetic'


y = 'bg_total_MARKER' 
one = df.groupby(['Condition'])[y].mean().tolist()[0]
two = df.groupby(['Condition'])[y].mean().tolist()[1]

sns.set_style("whitegrid")
plt.figure(figsize= (8,11))

ax = sns.violinplot(data=df, x="Condition", y=y, hue="Condition", order= [ 'Non-Diabetic','Diabetic',],hue_order=['Non-Diabetic','Diabetic'],log_scale=True,inner=None)#,inner_kws=dict(box_width=12, whis_width=2, color=".75"))
yposlist = df.groupby(['Condition'])[y].median().tolist()
xx = yposlist[1]
yposlist[1] = yposlist[0]
yposlist[0] = xx
xposlist = range(len(yposlist))
stringlist = ['Mdn = ' + str(int(np.round(two,0))),'Mdn = ' + str(int(np.round(one,0)))]


sns.boxplot(x='Condition', y=y, data=df, width=0.1,
            boxprops={'zorder': 2},color="gray", ax=ax,log_scale=True,showfliers=False,whis=(5,95))

for i in range(len(stringlist)):
    ax.text(xposlist[i] + .06, yposlist[i]-1000, stringlist[i],color='black', fontsize=24)
ax.set_xlabel('Condition',fontsize=18)
ax.set_ylabel('Total Normalized ' + MARKER + ' (log10 scale)',fontsize=18)

ax.tick_params(axis='x', labelsize=18)
plt.savefig("tr.pdf",format='pdf', bbox_inches="tight")
plt.show()


mean1 = np.mean(df[df['Folder'] == "Folder1"][y])
mean2 = np.mean(df[df['Folder'] == "Folder2"][y])
print("Mean Folder1:", np.round(mean1,2), "  Mean Folder2:",np.round(mean2,2))

res = stats.ranksums(df[df['Folder'] == "Folder1"][y], df[df['Folder'] == "Folder2"][y],nan_policy='omit')
print("T Statistic:", np.round(res.statistic,5))
#print("P Value:",np.round(res.pvalue,12))
print("P Value:",np.format_float_scientific(res.pvalue, exp_digits=2,precision=12))
