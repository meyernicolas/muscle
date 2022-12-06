# Multivariate Sparse CLustering for Extremes - MUSCLE

The code is based on the study of extremal dependence through L^1 minimization developed in the article
Meyer and Wintenberger (2021). Multivariate Sparse CLustering for Extremes.

The first file consists in the algorithm MUSCLE which provides the extremal directions (clusters) of a sample X_1,...,X_n. We then illustrate our algorithm on several examples.

In the first example we deal with a Gaussian copula with a common correlation parameter ρ<1 and marginal distributions satisfying P(X > x) = 1/x, ie Pareto(1) distribution. We compare our approach with the DAMEX algorithm introduced by Goix et al.(2017) (see 'damex').

The second example is taken from Simpson et al. (2020) and consists of a mixture of Gaussian and extreme-value logistic distributions (in 'swt'). We compare our approach with the one introduced by Goix et al.(2017) and by Simpson et al.(2020).

Finally we apply MUSCLE to wind speed data and financial return data.
