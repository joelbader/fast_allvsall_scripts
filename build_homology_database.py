#!/usr/bin/python3

import glob
import pickle
import pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess
import time
import parasail
from queue import Queue
import threading
import os

#Parameters and filenames:
num_threads=8
req_match_fraction = 0.50 #This fraction of query length must match to any similar sequences in order to be included in the pickle.
num_entries = -1 #Number of proteins to query from each genome. Set this to -1 to run through all genes. Useful for doing small testing runs

#output_pickle_filename = glob.glob("*.pickle")[0]
output_pickle_filename = 'mhdb_'+str(round(100*req_match_fraction))+'.pickle'

def main():
    #Get names of all *.fa files in current directory 
    genomes = glob.glob("genomes/*.fa")
    #If no database (pickle) is provided then create it, otherwise read it in
    try:
        #Read it in
        fp = open(output_pickle_filename, 'rb')
        db = pickle.load(fp)
        fp.close()
    except FileNotFoundError:
        #Create the dataframe
        db = pd.DataFrame()

#    if pick:
#        #Read it in
#        with open(pick[0], 'rb') as fp:
#            db = pickle.load(fp)
#        fp.close()
#    else:
#        #Create the dataframe
#        db = pd.DataFrame()
    
    #Compare current database with available *.fa - if any rows/columns don't exist then expand the matrix to include them. If any are extraneous then delete them.
    if not (set(db.columns) == set(db.index)):
        exit("Pickle in directory does not have the same rows as columns. This should not happen - please correct")
    
    if not bool(set(db.columns).issubset(genomes)):
        #Remove db.columns not in genomes
        to_drop = set(db.columns).difference(genomes)
        db.drop(to_drop,axis = 0)
        db.drop(to_drop,axis = 1)
    if not bool(set(genomes).issubset(db.columns)):
        #Add genomes (subtract db.columns) to db
        new = set(genomes).difference(db.columns)
        for c in new:
            db[c] = np.nan
            db.loc[c] = np.nan
    
    db = db.astype('object')
    #Iterate over missing elements (call add_genome()) - export to pickle after each genome_vs_genome
    q = Queue()
    #Locks necessary to protect writing out
    db_write_lock = threading.Lock()
    for c in db.columns:
        for r in db.index:
            if pd.isnull(db.at[r,c]):
                q.put((r,c))

    start = time.time()
    for x in range(num_threads):
        t=threading.Thread(target=threader, args=(q,db,db_write_lock), daemon=True)
        t.start()

    q.join()
    print('Entire job took (sec):',time.time()-start)

def threader(q,db,db_write_lock):
    while True:
        gen1,gen2 = q.get()
        temp_dict = genome_vs_genome(gen1,gen2)
        with db_write_lock:
            db.at[gen1,gen2] = temp_dict
            with open(output_pickle_filename, 'wb') as fp:
                pickle.dump(db, fp)
        q.task_done()

def genome_vs_genome(gen1, gen2):
    homolog_db = dict()

    #Define names of files to avoid thread collisions
    fname = gen1[8:-3] + '_' + gen2[8:-3]
    query_gene = 'query_gene_'+fname+'.fa'
    opal_output = 'OPAL_OUTPUT_'+fname+'.txt'
    parse_out = 'PARSE_OUT_'+fname+'.fa'

    iden_matrix = parasail.Matrix("WM_IDENTITY_MATRIX_parasail.txt")
    i = 0
    for gene in SeqIO.parse(gen1, "fasta"):
        gene.id = gene.id.split("|")[0] #Necessary for H37Rv genome
        #print(threading.current_thread().name, gen2)
        SeqIO.write(gene, query_gene, "fasta")
        #Call query_gene against gen2 with opal
        #t0 = time.time()
        subprocess.run('./opal_aligner -o 0 -e 0 -a NW -f Identity_score_matrix.txt -x 1 '+query_gene+' '+gen2+' > '+opal_output, shell=True)
        #t1 = time.time()
        #Read in opal output
        subprocess.run('./parse_opal.py '+opal_output+' '+gen2+' '+parse_out+' -t '+str(req_match_fraction), shell=True)
        #t2 = time.time()
        #Follow-up opal hits with parasail to count matches
        profile = parasail.profile_create_stats_16(str(gene.seq), iden_matrix)
        homolog_ls = list()
        for s2 in SeqIO.parse(parse_out, "fasta"):
            s2.id = s2.id.split("|")[0] #Necessary for H37Rv genome
            #print("1"+gene.id,s2.id)
            result = parasail.nw_stats_scan_profile_16(profile,str(s2.seq),20,1)
            #print("2"+gene.id,s2.id)
            match_fraction = result.matches/len(gene.seq)
            if match_fraction > req_match_fraction:
                homolog_ls.append((s2.id,match_fraction))
        #t3 = time.time()
        #print("t1-t0: ",t1-t0,"t2-t1: ",t2-t1,"t3-t2: ",t3-t2)

        homolog_db[gene.id]=homolog_ls
        i += 1
        if i == num_entries:
            break
    os.remove(query_gene)
    os.remove(opal_output)
    os.remove(parse_out)

    return(homolog_db)

if __name__ == '__main__':
    main()
