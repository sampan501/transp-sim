import numpy as np

def get_ave_distance(chrom):
    dist = 0
    for i in range(chrom.shape[0]-1):
        dist += abs(chrom[i]-chrom[i+1])
    return dist/(chrom.shape[0]-1)
def check_replication(list_chroms,v):
    dist = np.zeros(list_chroms.shape[0])
    for i,chrom in enumerate(list_chroms):
        dist[i] = get_ave_distance(chrom)
    # PLAY WITH THRESHOLD
    threshold = max(dist)*0.1
    if v:
        print('subgenome_1 average distance between transposons',dist[0])
        print('subgenome_2 average distance between transposons',dist[1])
        print('subgenome_3 average distance between transposons',dist[2])
    # check if any have a similar average distance between transposons
    sub1_2 = (print('Subgenome 1 and 2 originate from a single genome') if abs(dist[0]-dist[1]) < threshold else False)
    sub1_3 = (print('Subgenome 1 and 3 originate from a single genome') if abs(dist[0]-dist[2]) < threshold else False)
    sub2_3 = (print('Subgenome 2 and 3 originate from a single genome') if abs(dist[1]-dist[2]) < threshold else False)

def get_transposons_in_common(chrom1,chrom2,r):
    chrom2 = [round(t,r) for t in chrom2]
    transposons_in_common = []
    for transposon in chrom1: # look at each transposon
        t = round(transposon,r)
        if t in chrom2: # count transposons in common
            transposons_in_common.append(transposon)
    return len(transposons_in_common)

def find_most_ancestral(list_chroms,v,r):
    common_t1_2 = get_transposons_in_common(np.array(list_chroms[0]),np.array(list_chroms[1]),r)
    common_t1_3 = get_transposons_in_common(np.array(list_chroms[0]),np.array(list_chroms[2]),r)
    common_t2_3 = get_transposons_in_common(np.array(list_chroms[1]),np.array(list_chroms[2]),r)
    
    ave_1 = (common_t1_2+common_t1_3)/2
    ave_2 = (common_t1_2+common_t2_3)/2
    ave_3 = (common_t2_3+common_t1_3)/2
    
    max_c = max(ave_1,ave_2,ave_3)
    if max_c == ave_1:
        print('Subgenome 1 is the most ancestral')
        most_ancestral = list_chroms[0]
    elif max_c == ave_2:
        print('Subgenome 2 is the most ancestral')
        most_ancestral = list_chroms[1]
    else:
        print('Subgenome 3 is the most ancestral')
        most_ancestral = list_chroms[2]
#     return most_ancestral, max(common_t1_2,common_t1_3,common_t2_3)