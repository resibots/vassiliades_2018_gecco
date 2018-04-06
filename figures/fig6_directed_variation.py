import numpy as np
import matplotlib.pyplot as plt
params = {'text.usetex': True,
          'text.latex.preamble': [r'\usepackage{cmbright}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)

dimensionality = 2

p1 = [0.4,0.4]
p2 = [0.8,0.8]
p3 = [0.45,0.45]
p4 = [0.8,0.4]


sigma1 = 0.02   # for width: 0 creates a line 
sigma2 = 0.0    # for direction: 0 creates isotropic

def sample(p1, p2, n=2000):
    # Find a direction vector
    direction1 = np.array(p2) - np.array(p1)
    # print direction1

    points = []
    for k in range(n):
        point = p1 + sigma1 * np.random.randn(len(p1)) + sigma2 * direction1 * np.random.randn()
        points.append(point)

    return np.array(points)


fig_number_size = 22
font_size = 18
sigma1_pos = (0.55, 0.15)
sigma2_pos = (0.55, 0.05)

fig, axes = plt.subplots(1, 3, figsize=(8,2.6), facecolor='white', edgecolor='white')
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.05, hspace=0.05)

axes[0].annotate('A (Iso)', xy=(0.03,0.9), xycoords=('axes fraction', 'axes fraction'), size=fig_number_size, fontweight='bold')
axes[1].annotate('B (LineDD)', xy=(0.03,0.9), xycoords=('axes fraction', 'axes fraction'), size=fig_number_size, fontweight='bold')
axes[2].annotate('C (Iso+LineDD)', xy=(0.03,0.9), xycoords=('axes fraction', 'axes fraction'), size=fig_number_size, fontweight='bold')

for ax in axes:
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.tick_params(axis='both', length=0, labelsize=15)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
        ax.spines[axis].set_linewidth(4)

p = p2
sigma1 = 0.02   # for width: 0 creates a line 
sigma2 = 0.0    # for direction: 0 creates isotropic
sigma3 = 0.0
offspring1 = sample(p1,p)
ax = axes[0]
ax.annotate('p1',xy=tuple(p1),xycoords='data',xytext=(0.1, 0.5),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), size=font_size)
ax.annotate('p2',xy=tuple(p),xycoords='data',xytext=(0.5, 0.8),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), size=font_size)
ax.scatter(p1[0], p1[1], s=100, zorder=10, marker='.', color='black')
ax.scatter(p[0], p[1], s=100, zorder=10, marker='.', color='black')
ax.scatter(offspring1[:,0], offspring1[:,1], s=20, marker='.', color='red', alpha=0.5)
ax.annotate('$\sigma_1=$ '+str(sigma1),xy=sigma1_pos,xycoords='data',xytext=sigma1_pos,textcoords='data', size=font_size)
ax.annotate('$\sigma_2=$ '+str(sigma2),xy=sigma2_pos,xycoords='data',xytext=sigma2_pos,textcoords='data', size=font_size)


sigma1 = 0.0   # for width: 0 creates a line 
sigma2 = 0.2    # for direction: 0 creates isotropic
sigma3 = 0.2    # for direction: 0 creates isotropic
offspring1 = sample(p1,p)
ax = axes[1]
ax.annotate('p1',xy=tuple(p1),xycoords='data',xytext=(0.1, 0.5),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), size=font_size)
ax.annotate('p2',xy=tuple(p),xycoords='data',xytext=(0.5, 0.8),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), size=font_size)
ax.scatter(p1[0], p1[1], s=100, zorder=10, marker='.', color='black')
ax.scatter(p[0], p[1], s=100, zorder=10, marker='.', color='black')
ax.scatter(offspring1[:,0], offspring1[:,1], s=20, marker='.', color='red', alpha=0.5)
ax.annotate('$\sigma_1=$ '+str(sigma1),xy=sigma1_pos,xycoords='data',xytext=sigma1_pos,textcoords='data', size=font_size)
ax.annotate('$\sigma_2=$ '+str(sigma2),xy=sigma2_pos,xycoords='data',xytext=sigma2_pos,textcoords='data', size=font_size)


sigma1 = 0.02   # for width: 0 creates a line 
sigma2 = 0.2    # for direction: 0 creates isotropic
offspring1 = sample(p1,p)
ax = axes[2]
ax.annotate('p1',xy=tuple(p1),xycoords='data',xytext=(0.1, 0.5),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), size=font_size)
ax.annotate('p2',xy=tuple(p),xycoords='data',xytext=(0.5, 0.8),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle="arc3"), size=font_size)
ax.scatter(p1[0], p1[1], s=100, zorder=10, marker='.', color='black')
ax.scatter(p[0], p[1], s=100, zorder=10, marker='.', color='black')
ax.scatter(offspring1[:,0], offspring1[:,1], s=20, marker='.', color='red', alpha=0.5)
ax.annotate('$\sigma_1=$ '+str(sigma1),xy=sigma1_pos,xycoords='data',xytext=sigma1_pos,textcoords='data', size=font_size)
ax.annotate('$\sigma_2=$ '+str(sigma2),xy=sigma2_pos,xycoords='data',xytext=sigma2_pos,textcoords='data', size=font_size)

fig.savefig('figure6.pdf')
# plt.show()
