import sys
import pandas as pd

path = sys.argv[1]
meg_file = sys.argv[2]

df = pd.read_csv(path , index_col=0)
all_files = list(df.columns)
s = 0
df = df.applymap(str)
with open(meg_file , 'w') as f:
    f.write("#mega\n")
    f.write("!TITLE: " + meg_file + ";" + '\n')
    f.write("!Format DataType=distance NTaxa=" + str(len(all_files)) + ';' + '\n')
    f.write(' \n')

    for i in all_files:
        i = i[11:]
        str = i.replace(',_complete_genom)', '')
        str = str.replace(',_complete_cd)' , '')
        str = str.replace(',_complete_CD)' , '')
        str = str.replace(',_complete_sequenc)' , '')
        str = str.replace('#' , '_')
        str = str.replace("'" , '')
        str = str.replace('((' , '(')
        f.write("#" + str + '\n')
    f.write('\n')

    for a in all_files:
        for b in all_files:
            if df[b][a] != 'nan':
                f.write(df[b][a] + '\t')
            else:
                f.write('\t')
        f.write('\n')
print(s)