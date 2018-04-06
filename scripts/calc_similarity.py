import numpy as np
import glob
import sys
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm # a progress bar

# Usage: python calc_similarity.py <data_dir> <experiment>
# Examples:
# python calc_similarity.py arm_gecco_variationdgx arm
# python calc_similarity.py hexa_gecco_variationdgx hexa
# python calc_similarity.py schwefel_variationdgx schwefel

archives_dir = sys.argv[1]
experiment_name = sys.argv[2] # arm, hexa, schwefel

parallel = True
num_cores = multiprocessing.cpu_count()

if experiment_name == 'arm':
    behavior_dim = 2
    genotype_dim = 12
    phenotype_range = [-np.pi, np.pi]
elif experiment_name == 'hexa':
    behavior_dim = 6
    genotype_dim = 36
    phenotype_range = [0.0, 1.0]
elif experiment_name == 'schwefel':
    behavior_dim = 2
    genotype_dim = 100
    phenotype_range = [-5.0, 5.0]
else:
    print "Invalid experiment name"
    exit()

genotype_start_index = 1+behavior_dim
genotype_end_index = genotype_start_index+genotype_dim
max_dist = np.linalg.norm(np.array([phenotype_range[1]]*genotype_dim) - np.array([phenotype_range[0]]*genotype_dim))

def load_archives(archives_dir, generation):
    f_list = glob.glob(archives_dir + '/*/*/archive_' + str(generation) + '.dat')
    archives = []
    for f in f_list:
        archives.append(np.loadtxt(f))
    
    return archives

def get_similarity(archive):
    genotypes = archive[:,genotype_start_index:genotype_end_index]

    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=len(genotypes), algorithm='ball_tree').fit(genotypes)
    distances, indices = nbrs.kneighbors(genotypes)

    # calculate the similarity
    mean_dist = np.mean(distances)
    similarity = 1.0 - mean_dist/max_dist

    return similarity

def load_process_output(archives_dir, generation):
    print "Generation",str(generation)

    archives = load_archives(archives_dir, generation)

    if parallel:
        similarities = Parallel(n_jobs=num_cores)(delayed(get_similarity)(archive) for archive in tqdm(archives))

        filename = archives_dir + '/' + experiment_name + '_similarity_' + str(gen) + '.dat'
        with open(filename, 'w') as f:
            for similarity in similarities:
                f.write(str(similarity) + "\n")
    else:
        filename = archives_dir + '/' + experiment_name + '_similarity_' + str(gen) + '.dat'
        with open(filename, 'a') as f:
            for i,archive in enumerate(archives):
                print 'writing',i
                similarity = get_similarity(archive)
                f.write(str(similarity) + "\n")
        

if experiment_name == 'arm' or experiment_name == 'schwefel':
    gen = 0
    load_process_output(archives_dir, gen)
    gen = 500
    load_process_output(archives_dir, gen)
    gen = 998
    load_process_output(archives_dir, gen)
elif experiment_name == 'hexa':
    gen = 0
    load_process_output(archives_dir, gen)
    gen = 2500
    load_process_output(archives_dir, gen)
    gen = 5000
    load_process_output(archives_dir, gen)
