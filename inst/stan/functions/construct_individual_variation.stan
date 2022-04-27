// Construct individual variation terms
// 
// Combines standard normal variables per individual and term with standard
// deviation effects per term. Optionally supports uncorrelated random effects,
// correlated random effects, and no random effects
// 
// @param ieta A matrix of standard normal effects with P rows and K columns.
// 
// @param beta A vector of standard deviations indicating the population
// distribution for each individual.
// 
// @param P Integer indicating the number of individuals.
// 
// @param K Integer indicating the number of terms
//
// @param corr Binary indicating if variation should be correlated.
//
// @param var_exists Binary indicating if variation should be present
//
// @return A matrix of individual level effects based on population-level
// variation and optionally correlation.
// 
// @author Sam Abbott
matrix construct_individual_variation(matrix ieta, vector beta, int P, int K,
  int corr, matrix L_Omega, int var_exists
) {
    matrix[P, K] eta;
    if (corr) {
      // Cholesky factor of the covariance matrix
      matrix[K, K] L;
      L = diag_pre_multiply(beta, L_Omega);

      // Calculate per-infection correlated effects
      eta = (L * ieta)';
    }else{
      if (var_exists) {
        // All effects are independent
        for (i in 1:K) {
          eta[1:P, i] = to_vector(ieta[i, 1:P]) * beta[i];
        }
      }else{
        // No infection level differences
        eta = rep_matrix(0, P, K);
      }
    }
    return(eta);
  }