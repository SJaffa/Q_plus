# Q_plus
Q+: an algorithm for statistically quantifying sub-structure in 
star clusters

# INPUTS:
x,y: two-dimensional position of stars in parsecs

# OUTPUTS
filename
(Possbile Q warning if centrally concentrated)
PC1 and PC2 values
(Possibe edge of parameter space warning)
[D,D_err],[L,L_err],[C,C_err]

Where D,L,C areestimated fractal dimension, 
number of levels, density scaling exponent
and D_err,L_err,C_err are errors on these.

# DOCUMENTATION
This algorithm:
	-builds the complete graph and minimum spanning tree of the points,
	-derives various statistics of structure from these,
	-compares these statistics to a set with known structure, and
	-estimates the most likely structure parameters and their error.
See paper for a full explanation of the algorithm.

Usage:

In a terminal run:
python gridrecon.py <filename> <columns>

Where <filename> is the path to the text file you want to run 
and <columns> is a list of the indices if the x and y positions.
Example:

python gridrecon.py /home/sjaffa/Dropbox/Data/IC348.txt [0,1]

