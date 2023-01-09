# Multivariate Sparse CLustering for Extremes - MUSCLE

The code is based on the study of extremal dependence through L1 minimization developed in the article Meyer and Wintenberger (2021) entitled "Multivariate Sparse CLustering for Extremes".

The file 'muscle' contains the algorithm MUSCLE developed in this article which provides the extremal directions (clusters) of a sample X1,...,Xn.

This algorithm is compared with the DAMEX algorithm introduced by Goix et al.(2017) (see 'damex') and a procedure based on hodden regular variation proposed by Simpson et al. (2020) (see 'swt').
The comparison is achieved on two types of examples :
- in the first example we deal with a Gaussian copula with a common correlation parameter ρ<1 and marginal distributions satisfying P(X > x) = 1/x, ie Pareto(1) distribution.
- the second example is taken from Simpson et al. (2020) and consists of a mixture of Gaussian and extreme-value logistic distributions. We compare our approach with the one introduced by Goix et al.(2017) and by Simpson et al.(2020). Since the true distribution of Z is not known, we compute it via MOnte-Carlo simulations. (see 'monte_carlo_max_mixture').

Finally we apply the algorithm MUSCLE to windspeed data and financial return data.
