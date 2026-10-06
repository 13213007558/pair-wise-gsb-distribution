require 'distribution/kolmogorov_smirnov/ruby'

module Distribution
  # Exact finite-sample distribution of the one-sample two-sided
  # Kolmogorov-Smirnov statistic
  #   D_n = sup_x |F_n(x) - F(x)|
  # under a continuous null hypothesis (the distribution does not
  # depend on F).
  #
  # cdf, sf and p_value are exact for finite n (Pomeranz's recursion);
  # they are appropriate for small samples, where the asymptotic
  # Kolmogorov distribution is a poor approximation. The asymptotic
  # approximation is available separately as cdf_asymptotic and
  # sf_asymptotic and is clearly marked as approximate.
  #
  # == Usage:
  #    Distribution::KolmogorovSmirnov.cdf(0.3, 20)
  #    => 0.95761...
  #    Distribution::KolmogorovSmirnov.sf(0.3, 20)
  #    => 0.04238...
  #    Distribution::KolmogorovSmirnov.p_value(0.95, 20)
  #    => 0.29407...
  module KolmogorovSmirnov
    SHORTHAND = 'ks'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: cdf(d, n)
    # Exact P(D_n <= d), the finite-n cdf of the two-sided
    # Kolmogorov-Smirnov statistic under a continuous null.

    ##
    # :singleton-method: sf(d, n)
    # Exact P(D_n > d), the finite-n survival function. Computed with
    # a dedicated high-precision upper-tail algorithm, so very small
    # right tails stay accurate (it is not 1 - cdf throughout).

    ##
    # :singleton-method: p_value(pr, n)
    # Exact quantile of D_n: the smallest d with cdf(d) >= pr.

    ##
    # :singleton-method: cdf_asymptotic(d, n)
    # APPROXIMATE cdf from Kolmogorov's limiting distribution
    # (n -> infinity). Not exact for finite n.

    ##
    # :singleton-method: sf_asymptotic(d, n)
    # APPROXIMATE survival function from Kolmogorov's limiting
    # distribution (n -> infinity). Not exact for finite n.
  end
end
