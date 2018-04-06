import sys

if len(sys.argv) < 3:
    sys.exit('Usage: %s input_dir [VARIANT: spread/similarity]' % sys.argv[0])

input_dir = sys.argv[1]
variant = sys.argv[2]

if variant != 'spread' and variant != 'similarity':
    sys.exit('Variant should either be "spread" or "similarity"')

from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import brewer2mpl

bmap_schwefel = brewer2mpl.get_map('Blues', 'sequential', 3)
colors_schwefel = bmap_schwefel.mpl_colors
bmap_arm = brewer2mpl.get_map('Reds', 'sequential', 3)
colors_arm = bmap_arm.mpl_colors
bmap_hexa = brewer2mpl.get_map('Greens', 'sequential', 3)
colors_hexa = bmap_hexa.mpl_colors

experiments_titles = {
    'schwefel' : 'Schwefel Function',
    'arm' : 'Arm Repertoire',
    'hexa' : 'Hexapod Locomotion'
}

experiments_ticks = {
    'schwefel'  : ['0', '50k', '100k'],
    'arm'       : ['0', '50k', '100k'],
    'hexa'      : ['0', '250k', '500k']
}
experiments_filenames_spread = {
    'schwefel'  : ['schwefel_spread_0.dat', 'schwefel_spread_500.dat', 'schwefel_spread_998.dat'],
    'arm'       : ['arm_spread_0.dat', 'arm_spread_500.dat', 'arm_spread_998.dat'],
    'hexa'      : ['hexa_spread_0.dat', 'hexa_spread_2500.dat', 'hexa_spread_5000.dat']
}

experiments_filenames_similarity = {
    'schwefel'  : ['schwefel_similarity_0.dat', 'schwefel_similarity_500.dat', 'schwefel_similarity_998.dat'],
    'arm'       : ['arm_similarity_0.dat', 'arm_similarity_500.dat', 'arm_similarity_998.dat'],
    'hexa'      : ['hexa_similarity_0.dat', 'hexa_similarity_2500.dat', 'hexa_similarity_5000.dat']
}

# variant = 'spread'
# variant = 'similarity'

if variant == 'spread':
    experiments_filenames = experiments_filenames_spread
elif variant == 'similarity':
    experiments_filenames = experiments_filenames_similarity

# The data that I need to load has the format:
# Gen Median Perc25 Perc75
# The filename determines the task and metric (e.g., hexa_spread, arm_volume)

task_experiments = \
[
    [
        experiments_filenames['schwefel'][0],
        experiments_filenames['schwefel'][1],
        experiments_filenames['schwefel'][2]
    ],
    [
        experiments_filenames['arm'][0],
        experiments_filenames['arm'][1],
        experiments_filenames['arm'][2]
    ],
    [
        experiments_filenames['hexa'][0],
        experiments_filenames['hexa'][1],
        experiments_filenames['hexa'][2]
    ]
]

experiment_color = {}
experiment_color[experiments_filenames['schwefel'][0]] = colors_schwefel[0]
experiment_color[experiments_filenames['schwefel'][1]] = colors_schwefel[1]
experiment_color[experiments_filenames['schwefel'][2]] = colors_schwefel[2]
experiment_color[experiments_filenames['arm'][0]] = colors_arm[0]
experiment_color[experiments_filenames['arm'][1]] = colors_arm[1]
experiment_color[experiments_filenames['arm'][2]] = colors_arm[2]
experiment_color[experiments_filenames['hexa'][0]] = colors_hexa[0]
experiment_color[experiments_filenames['hexa'][1]] = colors_hexa[1]
experiment_color[experiments_filenames['hexa'][2]] = colors_hexa[2]


def do_boxplot(sp, data, experiments):
    global next_color
    bp = sp.boxplot(data, notch=0, sym='b+', vert=1, whis=1.5,
                    positions=None, widths=0.6)

    for i in range(len(bp['boxes'])):
        color = experiment_color[experiments[i]]
        box = bp['boxes'][i]
        box.set_linewidth(0)
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
            boxCoords = zip(boxX, boxY)
            boxPolygon = Polygon(boxCoords, facecolor=color, linewidth=0)
            sp.add_patch(boxPolygon)

    for i in range(0, len(bp['boxes'])):
        color = experiment_color[experiments[i]]
        bp['boxes'][i].set_color(color)
        # we have two whiskers!
        bp['whiskers'][i*2].set_color(color)
        bp['whiskers'][i*2 + 1].set_color(color)
        bp['whiskers'][i*2].set_linewidth(2)
        bp['whiskers'][i*2 + 1].set_linewidth(2)
        # top and bottom fliers
        bp['fliers'][i].set(markerfacecolor=color,
                                marker='o', alpha=0.75, markersize=6,
                                markeredgecolor='none')
        bp['medians'][i].set_color('black')
        bp['medians'][i].set_linewidth(3)
        # and 4 caps to remove
        for c in bp['caps']:
            c.set_linewidth(0)


def arrange_subplot(sp, title, labels):
    sp.set_title(title, fontsize=20)
    sp.spines['top'].set_visible(False)
    sp.spines['right'].set_visible(False)
    sp.spines['left'].set_visible(False)
    sp.get_xaxis().tick_bottom()
    sp.get_yaxis().tick_left()
    sp.tick_params(axis='x', direction='in', labelsize=20)
    sp.tick_params(axis='y', length=0, labelsize=20)
    sp.set_xticklabels(labels)
    sp.grid(axis='y', color="0.85", linestyle='--', linewidth=1.25)
    sp.set_axisbelow(True)


def load_and_plot(task, ax):
    task_data = []
    for i, exp in enumerate(experiments_filenames[task]):
        # Load Data
        print exp
        d = np.loadtxt(input_dir + '/' + exp)
        task_data.append(d)
        
    do_boxplot(ax, task_data, experiments_filenames[task])
    arrange_subplot(ax, experiments_titles[task], experiments_ticks[task])
    print np.shape(task_data)
    median = np.median(task_data, axis=1)
    print np.shape(median)
    # exit()

    ax.plot([1,2,3], median, color='black')

    print "Statistical Test:",task
    num_exp = len(task_data)
    for i in range(0,num_exp-1):
        for j in range(i+1,num_exp):
            z,p = scipy.stats.mannwhitneyu(task_data[i], task_data[j])
            p_value = p * 2.0
            print experiments_filenames[task][i],experiments_filenames[task][j],"p=",p_value


if variant == 'spread':
    ylabel = 'Genotypic Spread'
    filename = 'figure4'
elif variant == 'similarity':
    ylabel = 'Genotypic Similarity'
    filename = 'figure5'

fig, axes = plt.subplots(1, len(experiments_filenames), sharey=True, figsize=(10,4), facecolor='white', edgecolor='white')
axes[0].set_ylabel(ylabel, rotation="vertical", fontsize=22)
from matplotlib.ticker import MaxNLocator
axes[0].yaxis.set_major_locator(MaxNLocator(5))
fig.subplots_adjust(left=0.11, right=0.97, top=0.9, bottom=0.08, wspace=0.1, hspace=0.25)

load_and_plot('schwefel',axes[0])
load_and_plot('arm',axes[1])
load_and_plot('hexa',axes[2])


fig.savefig(filename + '.pdf')
# plt.show()
