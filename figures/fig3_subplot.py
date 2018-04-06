import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.spatial import Voronoi, voronoi_plot_2d
import sys
from matplotlib.ticker import FuncFormatter

cdict = {'red': [(0.0,  0.0, 0.0),
                 (0.33, 0.0, 0.0),
                 (0.66,  1.0, 1.0),
                 (1.0,  1.0, 1.0)],
         'blue': [(0.0,  0.0, 0.0),
                  (0.33, 1.0, 1.0),
                  (0.66,  0.0, 0.0),
                  (1.0,  0.0, 0.0)],
         'green': [(0.0,  0.0, 0.0),
                   (0.33, 0.0, 0.0),
                   (0.66,  0.0, 0.0),
                   (1.0,  1.0, 1.0)]}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 256)

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def load_data(filename):
    print "Loading ",filename
    data = np.loadtxt(filename)
    num_joints = 2

    indices = data[:,:1]
    descriptors = data[:,1:3]
    genotypes = data[:,3:num_joints+3]
    fitnesses = data[:,num_joints+3:num_joints+4]

    return indices,descriptors,genotypes,fitnesses

def plot_task_space(ax, indices, descriptors, fitnesses, regions, vertices):
    norm = matplotlib.colors.Normalize(vmin=min(fitnesses), vmax=max(fitnesses))

    # colorize
    for i, region in enumerate(regions):
        polygon = vertices[region]
        ax.fill(*zip(*polygon), alpha=0.05, edgecolor='black', facecolor='white', lw=1)

    # print indices
    for i in range(len(indices)):
        region = regions[int(indices[i][0])]
        polygon = vertices[region]
        ax.fill(*zip(*polygon), alpha=0.9, color=my_cmap(norm(fitnesses[i]))[0])

    sc = ax.scatter(descriptors[:,0], descriptors[:,1], c=fitnesses, cmap=my_cmap, s=10, zorder=0)
    return sc


if len(sys.argv) < 3:
    sys.exit('Usage: %s centroids_file data_dir' % sys.argv[0])

points = np.loadtxt(sys.argv[1])
rowsPoints, colsPoints = np.shape(points)
print 'Centroids: rows=',rowsPoints,' columns=',colsPoints
# compute Voronoi tesselation
vor = Voronoi(points)
regions, vertices = voronoi_finite_polygons_2d(vor)

data_dir = sys.argv[2]
filenames = ["archive_0.dat", "archive_50.dat", "archive_500.dat", "archive_5000.dat"]

# Loading data
indices = []
descriptors = []
genotypes = []
fitnesses = []
for f in filenames:
    ind,desc,gen,fit = load_data(data_dir + '/' + f)
    indices.append(ind)
    descriptors.append(desc)
    genotypes.append(gen)
    fitnesses.append(fit)


# cmap=my_cmap
titles = ["Generation 0", "Generation 50", "Generation 500", "Generation 5000"]
title = "2-DOF Arm"
tick_size = 18

# Plot
fig, axes = plt.subplots(2, 4, figsize=(22, 10), facecolor='white', edgecolor='white')
axes[0,0].set_ylabel("Task Space", rotation="vertical", fontsize=22)
axes[0,0].get_yaxis().set_label_coords(-0.15,0.5)

axes[1,0].set_ylabel("Parameter Space", rotation="vertical", fontsize=22)
axes[1,0].get_yaxis().set_label_coords(-0.15,0.5)

fig.subplots_adjust(left=0.05, right=0.92, top=0.95, bottom=0.05, wspace=0.2, hspace=0.2)

# Task Space
for i, ax in enumerate(axes[0,:]):
    print "Plotting task space",i
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.set_ylim(-1,1)
    ax.set_xlim(-1,1)
    # ax.set_ylim(-np.pi,np.pi)
    # ax.set_xlim(-np.pi,np.pi)
    ax.set_title(titles[i], fontsize=22)
    sc = plot_task_space(ax, indices[i], descriptors[i], fitnesses[i], regions, vertices)

# Parameter Space
for i, ax in enumerate(axes[1,:]):
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    ax.set_ylim(-np.pi,np.pi)
    ax.set_xlim(-np.pi,np.pi)
    print "Plotting parameter space",i
    sc = ax.scatter(genotypes[i][:,0], genotypes[i][:,1], c=fitnesses[i], cmap=my_cmap, s=35, zorder=0, marker='.')

cbar_ax = fig.add_axes([0.95, 0.1, 0.02, 0.8])
cbar_ax.tick_params(labelsize=tick_size)
# fig.colorbar(sc, ax=axes.ravel().tolist())
fig.colorbar(sc, cax=cbar_ax)

fig.savefig('figure3.png')
