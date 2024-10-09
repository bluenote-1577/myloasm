import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#tsv = 'x_hu-99'
tsv = 'x_hu-98'
df = pd.read_csv(tsv, sep='\t')
#first column is type
#second column is dist

def overlap(r1, r2):
    spl1 = r1.split('-')
    spl2 = r2.split('-')
    r1x = int(spl1[0])
    r1y = int(spl1[1])
    r2x = int(spl2[0])
    r2y = int(spl2[1])
    
    ol1 = r1x >= r2x and r1x <= r2y
    ol2 = r1y >= r2x and r1y <= r2y
    contained1 = r1x <= r2x and r1y >= r2y
    contained2 = r2x <= r1x and r2y >= r1y
    return ol1 or ol2 or contained1 or contained2

df.columns = ['type', 'dist', 'r1', 'r2']

# r1 is like 500-1000, r2 is also a range like 1000-2000
# create a new column checking if they overlap

df['overlap'] = df.apply(lambda x: overlap(x['r1'], x['r2']), axis=1)
df['final_type'] = df.apply(lambda x: x['type'] + str(x['overlap']), axis=1)

print(df)
df.to_csv('x_hu-99.8.df', sep='\t', index=False)
sns.histplot(df, x='dist', hue='final_type', bins=20, kde=True)
plt.show()

