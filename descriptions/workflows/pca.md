Principal Component Analysis
---------------

Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert
a set of observations of possibly correlated variables (entities each of which takes on various numerical values)
into a set of values of linearly uncorrelated variables called principal components.

The calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix,
not by using eigen on the covariance matrix. This is generally the preferred method for numerical accuracy.