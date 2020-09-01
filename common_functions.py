from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


pow4 = [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536,
    262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824, 4294967296,
    17179869184, 68719476736, 274877906944, 1099511627776, 4398046511104, 17592186044416,
    70368744177664, 281474976710656, 1125899906842624, 4503599627370496, 18014398509481984,
    72057594037927936, 288230376151711744, 1152921504606846976, 4611686018427387904]

pow4_ACGT = {'A' : [i*0 for i in pow4],
             'C' : [i*1 for i in pow4],
             'G' : [i*2 for i in pow4],
             'T' : [i*3 for i in pow4]}

V = {'A':0 , 'C':1 , 'G':2 , 'T':3}
RC = {'A':'T' , 'C':'G' , 'G':'C' , 'T':'A'}

def get_name(description):
    name = description[description.find('.')+3:description.find(",")]
    name = name[name.find("(")+1:name.find(")")]
    if name.find('(') != -1:
        name += ')'
    name = name.replace(' ', '_')
    name = name.replace('/', '_')
    print(name)
    return name

def get_info_fasta(path):
    fasta = SeqIO.parse(path, 'fasta')
    name = None
    s = 0
    for record in fasta:
        if not name:
            name = get_name(record.description)
        s += len(record.seq)
    return s, name


def reversed_comp(seq):
    rc = ''
    for i in range(len(seq)-1, -1, -1):
        rc += RC[seq[i]]
    return rc


def hashit(kmer, z, key_len):
    k = len(kmer)
    hashed = 0
    for i in z:
        if i >= key_len:
            hashed += pow4_ACGT[kmer[i]][k - (i+1)]
    return hashed


def draw_dot_plot(dot_plot, path, db_name, query_name, db_len, query_len, dim, label):
    
    plt.figure(figsize=(10,10))
    plt.imshow(-dot_plot, origin="upper", cmap='gray', extent=[0,query_len,0,db_len])
    plt.xlabel(query_name)
    plt.ylabel(db_name)
    plt.title("Dotplot for "+label)
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.savefig(path + db_name + '_' + query_name + '_' + label + '_dot_plot.png')
    plt.close()


def z_to_label(z, k):
    label = ''
    for i in range(k):
        if i in z:
            label += '1'
        else:
            label += '*'
    return label