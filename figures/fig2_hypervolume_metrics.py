import numpy as np
import matplotlib.pyplot as plt

dim = 2

max_dist = np.linalg.norm(np.array([1]*dim))

def calc_spread_similarity(points):
    from sklearn.neighbors import NearestNeighbors
    nbrs = NearestNeighbors(n_neighbors=len(points), algorithm='ball_tree', n_jobs=8).fit(points)
    distances, indices = nbrs.kneighbors(points)

    # get distances to nearest neighbors
    distance_to_closest_neighbor = distances[:,1]
    mean_dist_to_closest_neighbor = np.mean(distance_to_closest_neighbor)
    mean_distances = np.mean(distances)

    spread = mean_dist_to_closest_neighbor / max_dist
    similarity = 1.0 - mean_distances/max_dist

    return spread, similarity

height = 5
num_plots = 3
fig, axes = plt.subplots(1, num_plots, figsize=(num_plots*height,height), facecolor='white', edgecolor='white')
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.05, hspace=0.05)

data_high_spread_low_similarity = np.random.rand(200,2)
data_low_spread_high_similarity = 0.5 + np.random.randn(200,2) * 0.02

data_low_spread_low_similarity = np.zeros((200,2))
data_low_spread_low_similarity[:50,:] = np.array([0.05,0.05]) + np.random.randn(50,2) * 0.01
data_low_spread_low_similarity[50:100,:] = np.array([0.05,0.95]) + np.random.randn(50,2) * 0.01
data_low_spread_low_similarity[100:150,:] = np.array([0.95,0.05]) + np.random.randn(50,2) * 0.01
data_low_spread_low_similarity[150:,:] = np.array([0.95,0.95]) + np.random.randn(50,2) * 0.01


spread0, similarity0 = calc_spread_similarity(data_high_spread_low_similarity)
print spread0, similarity0
spread1, similarity1 = calc_spread_similarity(data_low_spread_high_similarity)
print spread1, similarity1
spread2, similarity2 = calc_spread_similarity(data_low_spread_low_similarity)
print spread2, similarity2


for ax in axes:
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.tick_params(axis='both', length=0, labelsize=15)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
        ax.spines[axis].set_linewidth(4)

axes[0].scatter(data_high_spread_low_similarity[:,0], data_high_spread_low_similarity[:,1])
axes[1].scatter(data_low_spread_high_similarity[:,0], data_low_spread_high_similarity[:,1])
axes[2].scatter(data_low_spread_low_similarity[:,0], data_low_spread_low_similarity[:,1])

# plt.show()
fig.savefig('figure2.pdf')
