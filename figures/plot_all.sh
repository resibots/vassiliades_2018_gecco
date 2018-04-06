# Figure 1
echo '\nFigure 1\n'
inkscape fig1_conceptual.svg -D --export-pdf=figure1.pdf

# Figure 2
echo '\nFigure 2\n'
python fig2_hypervolume_metrics.py

# Figure 3
echo '\nFigure 3\n'
python fig3_subplot.py ../centroids/centroids_10000_2_normalized_-1.0_1.0.dat ../results/2_joints/.

# Figure 4
echo '\nFigure 4\n'
python fig4_5_spread_similarity_boxplots.py ../results/spread_similarity/. spread

# Figure 5
echo '\nFigure 5\n'
python fig4_5_spread_similarity_boxplots.py ../results/spread_similarity/. similarity

# Figure 6
echo '\nFigure 6\n'
python fig6_directed_variation.py

# Figure 7
echo '\nFigure 7\n'
python fig7_results.py ../results/.
