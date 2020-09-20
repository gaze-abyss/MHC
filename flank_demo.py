import pandas as pd
ensemble = pd.read_csv("C:/Users/hjwang/Desktop/homo_pfams_202003101726.csv")
def regularaa(seq):
    '''
    filter for natural amino acid
    '''
    aalist = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    try:
        for j in seq:
            if j in aalist:
                continue
            else:
                return False
        return True
    except:
        return False

def judgement(seq_file,flank_file):
    '''
    filter according to 1)natural peptide 2)sequence length 3)not blank
    '''
    judge = []
    for i in range(len(seq_file)):
        if len(seq_file[i]) in [8,9,10,11] and regularaa(seq_file[i]) and regularaa(flank_file[i]) and flank_file[i] !='':
            judge.append(True)
        else:
            judge.append(False)
    return judge

df = pd.read_csv("")
mix_protein = []

def match0flank(seq,seq_all,length = 30):
    seq_all = '-'*length+seq_all+'-'*length
    num = seq_all.find(seq)
    return seq_all[num-length:num+len(seq)+length]

mix_peptide = df["sequence"].values.tolist()

mix_flank = []
for i in range(len(df)):
    mix_flank.append(match0flank(mix_peptide[i],mix_protein[i],length=30))
df['protein'] = mix_protein
df['flanking_seq30'] = mix_flank

judge = judgement(mix_peptide,mix_flank)
df_filtered = df[judge]

df.to_csv('')