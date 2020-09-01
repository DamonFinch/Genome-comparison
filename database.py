import csv
from common_functions import *
import numpy as np
import pandas as pd


class Hit:
    def __init__(self, kmer_hash, kmer_pos):
        self.kmer_hash = kmer_hash
        self.kmer_pos = kmer_pos
        self.repetition = False
        self.pos_in_query = -1
        self.hit_counts = 0


class database:
    def __init__ (self, z, kmer_len=12, key_len=32):
        self.hit_table ={}
        self.key_len = key_len
        self.kmer_len = kmer_len
        self.z = z
        self.n_uniques = 0
        self.n_hits = 0
        self.ref_len = 0
        self.query_len = 0
        self.dot_plot = None

    def update_ref_len(self, ref_len):
        self.ref_len = ref_len

    def update_query_len(self, query_len):
        self.query_len = query_len

    def add (self, kmer, index):
        key = kmer[:self.key_len]
        if key in self.hit_table:
            if self.hit_table[key].repetition == False:
                self.hit_table[key].repetition = True
                self.n_uniques -= 1
            
        else:
            self.n_uniques += 1
            self.hit_table[key] = Hit(hashit(kmer, self.z, self.key_len), index)

    def add_query(self, kmer, index):
        key = kmer[:self.key_len]
        hash_mer = hashit(kmer, self.z, self.key_len)
        if key in self.hit_table:
            if (self.hit_table[key].repetition == False) and (hash_mer == self.hit_table[key].kmer_hash) :
                self.hit_table[key].hit_counts += 1
                self.hit_table[key].pos_in_query = index
                self.n_hits += 1


    def remove_repeats(self):
        for hit_key in list(self.hit_table):
            if self.hit_table[hit_key].repetition == True:
                del self.hit_table[hit_key]

    def remove_not_hits(self):
        for hit_key in list(self.hit_table):
            if self.hit_table[hit_key].hit_counts == 0:
                del self.hit_table[hit_key]

    def save(self, file):
        print("Saving data ...")
        try:
            with open(file, 'w') as csvfile:
                writer = csv.writer(csvfile)
                for hit_key in self.hit_table:
                    hit = self.hit_table[hit_key]
                    writer.writerow([hit.kmer_hash, hit.kmer_pos, hit.repetition, hit.pos_in_query, hit.hit_counts, hit_key])              
        except IOError:
            print("I/O error" + file)



    def load(self, csv_file):
        print("Loading data ...")
        try:
            with open(csv_file, 'r') as csvfile:  
                reader = csv.reader(csvfile)
                for row in reader:
                    kmer_hash = int(row[0])
                    index = int(row[1])
                    repetition = (row[2] == 'True')
                    pos_in_query = int(row[3])
                    hit_counts = int(row[4])
                    key = str(row[5])
                    if repetition == False: self.n_uniques += 1
                    if hit_counts == 1: self.n_hits += 1
                    self.hit_table[key] = Hit(kmer_hash, index)
                    self.hit_table[key].repetition = repetition
                    self.hit_table[key].pos_in_query = pos_in_query
                    self.hit_table[key].hit_counts = hit_counts
        except IOError:
            print("I/O error")


    def create_hitmatrix(self, dim):
        self.hit_matrix = np.zeros((dim, dim), dtype=np.int)
        for hit in self.hit_table.values():
            if hit.hit_counts  == 2 :
                redir_db = min(dim - 1, int(hit.kmer_pos * dim / self.ref_len * 1.0))
                redir_query = min(dim - 1, int(hit.pos_in_query * dim / self.query_len * 1.0))
                self.hit_matrix[redir_db][redir_query] += 1



    def create_dotplot(self, diag_len=4):
        dim = len(self.hit_matrix)

        score_density = np.zeros((dim, dim), dtype=np.int)
        aux_density = self.hit_matrix


        for i in range(dim):
            cmax_pos = np.argmax(aux_density[i, :])
            if aux_density[i, cmax_pos] > 0:
                score_density[i, :] = 0
                score_density[i, cmax_pos] = 1


        for i in range(dim):
            cmax_pos = np.argmax(aux_density[:, i])
            if aux_density[cmax_pos, i] > 0:
                score_density[:, i] = 0
                score_density[cmax_pos, i] = 1

        score_copy = np.copy(score_density)

        for i in range(5, dim - 5):
            for j in range(5, dim - 5):
                value = 0
                for w in range(int(-diag_len / 2), int(diag_len / 2)):
                    if score_density[i + w, j + w] > 0:
                        value += 1

                if value >= diag_len:
                    for k in range(1, 5):
                        score_copy[i + k, j + k] = 1
                        score_copy[i - k, j - k] = 1

        for i in range(5, dim - 5):
            for j in range(5, dim - 5):
                value = 0
                for w in range(int(-diag_len / 2), int(diag_len / 2)):
                    if score_density[i - w, j + w] > 0:
                        value += 1

                if value >= diag_len:
                    for k in range(1, 5):
                        score_copy[i - k, j + k] = 1
                        score_copy[i + k, j - k] = 1

        for i in range(dim): 
            for j in range(dim):
                min_i = max(1, i - 2)
                max_i = min(dim, i + 2)
                min_j = max(1, j - 2)
                max_j = min(dim, j + 2)

                value = np.sum(score_copy[min_i:max_i, min_j:max_j])

                if value < (int(diag_len / 2)):
                    score_copy[i, j] = 0

        self.dot_plot = score_copy
        return 

    def score(self, dist_th=1.5):
        dim = len(self.hit_matrix)
        score = 0

        dvec1 = abs(np.argmax(self.dot_plot[:, 1]) - np.argmax(self.dot_plot[:, 0]))
        dvec2 = abs(np.argmax(self.dot_plot[:, 2]) - np.argmax(self.dot_plot[:, 1]))
        dvec3 = abs(np.argmax(self.dot_plot[:, 3]) - np.argmax(self.dot_plot[:, 2]))

        for i in range(4, dim):
            distance = np.mean([dvec1, dvec2, dvec3])

            dvec1 = dvec2
            dvec2 = dvec3
            dvec3 = abs(np.argmax(self.dot_plot[:, i]) - np.argmax(self.dot_plot[:, i - 1]))

            if distance > dist_th or distance == 0:
                score += dim

        score /= (dim ** 2)
        return score
