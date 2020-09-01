from database import *
from common_functions import *
import os.path
from tqdm import tqdm
from pathlib import Path


def chromeister(path, query, ref, z, kmer_len=32, key_len=12, dim=1000, log=True):
    db = database(z=z, kmer_len=kmer_len, key_len=key_len)

    ref_path = path+ref
    ref_len, ref_name = get_info_fasta(ref_path)
    db.update_ref_len(ref_len)

    query_path = path+query
    query_len, query_name = get_info_fasta(query_path)
    db.update_query_len(query_len)

    label = z_to_label(z, kmer_len)
    res_path = path + 'results/db/'
    image_res_path = path + 'results/'

    ref_file = res_path + ref_name + '_' + label + '_unique_inexact_hits.csv'
    ref_query_file =res_path + ref_name + '_' + query_name + '_' + label + '_hits.csv'

    Path(path + '/results/db/').mkdir(parents=True, exist_ok=True)


    if os.path.isfile(ref_query_file):
        db.load(ref_query_file)
        print("[INFO] Found "+str(db.n_hits)+" hits between reference and query.")        

    else:
        if os.path.isfile(ref_file):
            db.load(ref_file)

        else:
            print('Reading reference fasta and making inexact matching database ...')
            ref = SeqIO.parse(ref_path, 'fasta')
            kmer = ''
            curr_pos = 0
            if log:
                t1 = tqdm(total=ref_len)

            for record in ref:
                sequence = record.seq
                for i in range(len(sequence)):
                    if len(kmer) == kmer_len:
                        db.add(kmer, curr_pos)
                        kmer = ''

                    if sequence[i] in V:
                        kmer += sequence[i]
                    
                    else:
                        kmer = ''
            
                    curr_pos += 1
                    if log:
                        t1.update()
            if log:
                t1.close()
            db.remove_repeats()
            db.save(ref_file)
        print('[INFO] Found '+str(db.n_uniques)+' uniques in reference.')

        print('Reading query fasta and findig hits between reference and query ...')
        query = SeqIO.parse(query_path, 'fasta')
        kmer = ''
        curr_pos = 0
        if log:
            t2 = tqdm(total=query_len)

        for record in query:
            sequence = record.seq
            for i in range(len(sequence)):
                if len(kmer) == kmer_len:
                    db.add_query(kmer, curr_pos)
                    db.add_query(reversed_comp(kmer), curr_pos)
                    kmer = kmer[5:]

                if sequence[i] in V:
                    kmer += sequence[i]
                
                else:
                    kmer = ''
        
                curr_pos += 1
                if log:
                    t2.update()
        if log:
            t2.close()
        db.remove_not_hits()
        db.save(ref_query_file)
        print("[INFO] Found "+str(db.n_hits)+" hits between reference and query.")

    db.create_hitmatrix(dim=dim)

    print("[INFO] Found "+str(np.sum(db.hit_matrix))+" hits for z ="+str(db.z))
    
    db.create_dotplot()

    print("[INFO] Found "+str(np.sum(db.dot_plot))+" hits after filtering.")

    score = db.score()
    
    draw_dot_plot(db.dot_plot, image_res_path , ref_name, query_name, ref_len, query_len, dim, label)

    return score, query_name, ref_name
