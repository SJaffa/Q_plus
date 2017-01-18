# Q_plus
Q+: an algorithm for statistically quantifying sub-structure in star clusters

# INPUTS:
x,y: two-dimensional position of stars in parsecs

# OUTPUTS
[D,L,C]: estimated fractal dimension, number of levels, density scaling exponent
[D_err,L_err,C_err]: errors on  estimated fractal dimension, number of levels, density scaling exponent

# DOCUMENTATION
This algorithm:
	-builds the complete graph and minimum spanning tree of the points,
	-derives various statistics of structure from these,
	-compares these statistics to a set with known structure, and
	-estimates the most likely structure parameters and their error.
See paper for a full explanation of the algorithm.


