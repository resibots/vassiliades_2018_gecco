import glob
from pylab import *
import brewer2mpl
import sys
from joblib import Parallel, delayed
import multiprocessing
import matplotlib.ticker as ticker
import scipy.stats

params = {'text.usetex': True,
          'text.latex.preamble': [r'\usepackage{cmbright}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)

experiment_color = {}

tasks_experiments = [
[
"schwefel_variationisolinedd",
"schwefel_variationlinedd",
"schwefel_variationline",
"schwefel_variationiso",
"schwefel_variationisodd",
"schwefel_variationisosa",
"schwefel_variationgc",
"schwefel_variationsbx"
],
[
"arm_variationisolinedd",
"arm_variationlinedd",
"arm_variationline",
"arm_variationiso",
"arm_variationisodd",
"arm_variationisosa",
"arm_variationgc",
"arm_variationsbx"
],
[
"hexa_variationisolinedd",
"hexa_variationlinedd",
"hexa_variationline",
"hexa_variationiso",
"hexa_variationisodd",
"hexa_variationisosa",
"hexa_variationgc",
"hexa_variationsbx"
]
]

titles = (
    "Iso+LineDD (proposed approach)",
    "LineDD",
    "Line",
    "Iso",
    "IsoDD",
    "IsoSA",
    "Global Correlation",
    "SBX"
    )

# Assign colors
bmap = brewer2mpl.get_map('Dark2', 'Qualitative', 8)
colors = bmap.mpl_colors

for i, task in enumerate(tasks_experiments):
    for j, exp in enumerate(task):
        experiment_color[exp] = colors[j]

def processLoadData(f):
    # <archive size> <mean fitness> <max fitness>
    data = np.loadtxt(f)
    return [data[:,2], data[:,4], data[:,5]]

def loadCVT(dir):
    data = []
    f_list = glob.glob(dir + '/*/*/progress_archive.dat')

    num_cores = multiprocessing.cpu_count()
    data = Parallel(n_jobs=num_cores)(delayed(processLoadData)(f) for f in f_list)

    min_nb_generations = 99999999999
    for d in data:
        if len(d) < min_nb_generations:
            min_nb_generations = len(d)
    for i,d in enumerate(data):
        data[i] = data[i][:min_nb_generations]

    return np.asarray(data)

def perc(data):
    num_archives, vals, num_gens = np.shape(data)
    assert vals == 3
    print "num_archives=",num_archives, " num_gens",num_gens
    median = np.zeros((num_gens, vals))
    perc_25 = np.zeros((num_gens, vals))
    perc_75 = np.zeros((num_gens, vals))
    for i in range(0, len(median)):
        median[i,:] = np.median(data[:,:,i], axis=0)
        perc_25[i,:] = np.percentile(data[:,:,i], 25, axis=0)
        perc_75[i,:] = np.percentile(data[:,:,i], 75, axis=0)
    return median, perc_25, perc_75

def arrange_subplot(ax,ylabel,xlabel='Evaluations (x$10^{3}$)', grid=True):
    ax.set_xlabel(xlabel, fontsize=20)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x', direction='out', labelsize=20)
    if grid:
        ax.tick_params(axis='y', length=0, labelsize=20)
        ax.grid(axis='y', color="0.85", linestyle='--', linewidth=1.25)
        ax.grid(axis='x', color="0.85", linestyle='--', linewidth=1.25)
    else:
        ax.tick_params(axis='y', direction='in', labelsize=20)

    ax.set_axisbelow(True)

def non_dominated(i, exp):
    if i == 0 and (exp == "schwefel_variationisolinedd" or exp == "schwefel_variationline" or exp == "schwefel_variationsbx"):
        return True

    if i == 1 and (exp == "arm_variationisolinedd" or exp == "arm_variationlinedd"):
        return True

    if i == 2 and (exp == "hexa_variationisolinedd" or exp == "hexa_variationline"):
        return True

    return False

def stars(p):
   if p < 0.0001:
       return "****"
   elif (p < 0.001):
       return "***"
   elif (p < 0.01):
       return "**"
   elif (p < 0.05):
       return "*"
   else:
       return "-"

# performs a statistical test between all pairs of solutions and all metrics
# at the specified generation and outputs the results in filename
def statistical_test(gen, filename):
    # <exp> <exp> <metric> <pvalue> <significance (number of stars)>
    with open(filename,'w') as f:
        # for each task
        for i in range(len(tasks_experiments)):

            # for each metric
            for j in range(3):

                # for each pair of experiments
                for k1, exp1 in enumerate(tasks_experiments[i]):
                    if len(data[exp1][0,j,:]) < gen:
                        continue
                    for k2 in range(k1+1, len(tasks_experiments[i])):
                        exp2 = tasks_experiments[i][k2]

                        # print exp1, exp2
                        f.write("<" + exp1 + ">\t<" + exp2 + ">\t")
                        if j == 0:
                            f.write("<archivesize>\t")
                        elif j == 1:
                            f.write("<meanfitness>\t")
                        elif j == 2:
                            f.write("<maxfitness>\t")

                        # Do a statistical test
                        data1 = data[exp1][:,j,gen]
                        data2 = data[exp2][:,j,gen]

                        z,p = scipy.stats.mannwhitneyu(data1, data2)
                        p_value = p * 2.0

                        f.write(" " + str(p_value) + "\t" + stars(p_value) + "\n")


# Create the plots
fig, axes = plt.subplots(4, 3, figsize=(14,13), facecolor='white', edgecolor='white')
fig.subplots_adjust(left=0.1, right=0.97, top=0.97, bottom=0.17, wspace=0.25, hspace=0.34)

axes[0,0].set_title('Schwefel Function', fontsize=22)
axes[0,1].set_title('Arm Repertoire', fontsize=22)
axes[0,2].set_title('Hexapod Locomotion', fontsize=22)

from matplotlib.ticker import MaxNLocator

labels_arm_schwefel = ['0','20','40','60','80','100']
labels_hexa = ['0', '100', '200', '300', '400', '500']

axes[0,0].set_ylim(2000,10000)
axes[0,0].yaxis.set_major_locator(MaxNLocator(4))
axes[0,0].set_xticklabels(labels_arm_schwefel)

axes[0,1].set_ylim(2000,7500)
axes[0,1].set_xticklabels(labels_arm_schwefel)

axes[0,2].set_ylim(2000,7500)
axes[0,2].set_xticklabels(labels_hexa)

axes[1,0].set_ylim(-40000,-1000)
axes[1,0].yaxis.set_major_locator(MaxNLocator(4))
axes[1,0].set_xticklabels(labels_arm_schwefel)

axes[1,1].set_ylim(-30,0)
axes[1,1].set_xticklabels(labels_arm_schwefel)

axes[1,2].set_ylim(0.1,0.43)
axes[1,2].yaxis.set_major_locator(MaxNLocator(5))
axes[1,2].set_xticklabels(labels_hexa)

axes[2,0].set_ylim(-2000,0)
axes[2,0].set_xticklabels(labels_arm_schwefel)
axes[2,1].set_ylim(-10,0)
axes[2,1].set_xticklabels(labels_arm_schwefel)
axes[2,2].set_ylim(0.5,1.35)
axes[2,2].yaxis.set_major_locator(MaxNLocator(5))
axes[2,2].set_xticklabels(labels_hexa)


axes[3,0].set_xlim(9690,9960)
axes[3,0].set_ylim(-10500,-1000)

axes[3,1].set_xlim(5900,7700)

axes[3,2].set_ylim(0.25,0.43)
axes[3,2].set_xlim(3500,7500)

axes[0,0].annotate('Archive Size',  # Your string
            # The point that we'll place the text in relation to
            xy=(0.015, 0.5),
            # Interpret the x as figure coords, and the y as axes coords
            xycoords=('figure fraction', 'axes fraction'),
            # The distance from the point that the text will be at
            xytext=(0, 0),
            # Interpret `xytext` as an offset in points...
            textcoords='offset points',
            # Any other text parameters we'd like
            size=20, rotation="vertical", ha="center", va="center")

axes[1,0].annotate('Mean Archive Fitness',  # Your string
            # The point that we'll place the text in relation to
            xy=(0.015, 0.5),
            # Interpret the x as figure coords, and the y as axes coords
            xycoords=('figure fraction', 'axes fraction'),
            # The distance from the point that the text will be at
            xytext=(0, 0),
            # Interpret `xytext` as an offset in points...
            textcoords='offset points',
            # Any other text parameters we'd like
            size=20, rotation="vertical", ha="center", va="center")

axes[2,0].annotate('Max Archive Fitness',  # Your string
            # The point that we'll place the text in relation to
            xy=(0.015, 0.5),
            # Interpret the x as figure coords, and the y as axes coords
            xycoords=('figure fraction', 'axes fraction'),
            # The distance from the point that the text will be at
            xytext=(0, 0),
            # Interpret `xytext` as an offset in points...
            textcoords='offset points',
            # Any other text parameters we'd like
            size=20, rotation="vertical", ha="center", va="center")

axes[3,0].annotate('Mean Archive Fitness',  # Your string
            # The point that we'll place the text in relation to
            xy=(0.015, 0.5),
            # Interpret the x as figure coords, and the y as axes coords
            xycoords=('figure fraction', 'axes fraction'),
            # The distance from the point that the text will be at
            xytext=(0, 0),
            # Interpret `xytext` as an offset in points...
            textcoords='offset points',
            # Any other text parameters we'd like
            size=20, rotation="vertical", ha="center", va="center")

# Get the data
dataDir = sys.argv[1]

# Load data
data = {}
for task in tasks_experiments:
    min_nb_generations = 99999999999
    for exp in task:
        print 'Loading ',exp
        data[exp] = loadCVT(dataDir + '/' + exp + '/')
        print np.shape(data[exp]) # shape should be: nb_runs 3 nb_gens

        if data[exp].shape[2] < min_nb_generations:
            min_nb_generations = data[exp].shape[2]

# Calculate statistics
medians = {}
percentiles_25 = {}
percentiles_75 = {}

gen_for_statistics = 600

if gen_for_statistics < 0:
    gen_for_statistics = -1
    filename_statistics = dataDir + '/statistics_last_gen.dat'
    filename_statistical_test = dataDir + '/statistical_test_last_gen.dat'
else:
    filename_statistics = dataDir + '/statistics_gen_' + str(gen_for_statistics) + '.dat'
    filename_statistical_test = dataDir + '/statistical_test_gen_' + str(gen_for_statistics) + '.dat'


# Perform the statistical test between all experiments
print 'Performing statistical tests'
statistical_test(gen_for_statistics, filename_statistical_test)

with open(filename_statistics,'w') as f:
    for task in tasks_experiments:
        for exp in task:
            print 'Calculating statistics for ',exp
            data[exp] = data[exp][:,:,:min_nb_generations]
            print np.shape(data[exp])
            medians[exp], percentiles_25[exp], percentiles_75[exp] = perc(data[exp])

            if len(medians[exp][:,0]) >= gen_for_statistics:
                f.write(exp + " : median = " + str(medians[exp][gen_for_statistics,:]) + "\t perc25 = " + str(percentiles_25[exp][gen_for_statistics,:]) + "\t perc75 = " + str(percentiles_75[exp][gen_for_statistics,:]) + "\n")
                # f.write(exp + " : median = " + str(medians[exp][200,:]) + "\t perc25 = " + str(percentiles_25[exp][200,:]) + "\t perc75 = " + str(percentiles_75[exp][200,:]) + "\n")

ylabels = ['Archive Size', 'Mean Archive Fitness', 'Max Archive Fitness']

lines = []

# for each task
for i in range(len(tasks_experiments)):
    # for each metric (archivesize, meanfitness, maxfitness)
    for j in range(3):
        ax = axes[j,i]

        arrange_subplot(ax,ylabels[j])

        # Plot
        for k, exp in enumerate(tasks_experiments[i]):
            x = np.arange(0, len(medians[exp]))
            print 'i=',i,' j=',j,': plotting ',exp
            ax.fill_between(x, percentiles_25[exp][:,j], percentiles_75[exp][:,j], alpha=0.25, linewidth=0, color=experiment_color[exp])

            if exp == "schwefel_variationisolinedd" or exp == "arm_variationisolinedd" or exp == "hexa_variationisolinedd":
                l, = ax.plot(x, medians[exp][:,j], linewidth=2, color=experiment_color[exp], marker='o', markevery=int(len(x)/10), zorder=10)
            # elif exp == "schwefel_variationlinedd" or exp == "arm_variationlinedd" or exp == "hexa_variationlinedd":
            #     l, = ax.plot(x, medians[exp][:,j], linewidth=2, color=experiment_color[exp], marker='s', markevery=int(len(x)/10))
            # elif exp == "schwefel_variationline" or exp == "arm_variationline" or exp == "hexa_variationline":
            #     l, = ax.plot(x, medians[exp][:,j], linewidth=2, color=experiment_color[exp], marker='h', markevery=int(len(x)/10))
            else:
                l, = ax.plot(x, medians[exp][:,j], linewidth=2, color=experiment_color[exp])

            if len(lines) <= k:
                lines.append(l)

    # Plot the Pareto fronts
    ax = axes[3,i]
    for k, exp in enumerate(tasks_experiments[i]):
        # x is median of archive size, y is median of mean archive fitness
        point = [medians[exp][-1,0] , medians[exp][-1,1]]

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        if non_dominated(i, exp):
            ax.scatter(point[0], point[1], color=experiment_color[exp], s=70, marker='s', zorder=20)
        else:
            ax.scatter(point[0], point[1], color=experiment_color[exp], s=30)
        arrange_subplot(ax, 'Mean Archive Fitness', 'Archive Size', grid=False)
        ax.xaxis.set_tick_params(pad=-10)
        ax.axhline(point[1], color=experiment_color[exp], linestyle='--', alpha=0.5)
        ax.axvline(point[0], color=experiment_color[exp], linestyle='--', alpha=0.5)

fig.legend(tuple(lines), titles, loc='lower center', fontsize=20, ncol=3)

fig.savefig('figure7.pdf')
# plt.show()
