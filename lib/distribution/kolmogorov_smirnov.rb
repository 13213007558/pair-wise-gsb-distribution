require 'distribution/kolmogorov_smirnov/ruby'

module Distribution
  # Finite-sample distribution of the one-sample two-sided
  # Kolmogorov-Smirnov statistic under a fully specified continuous null.
  module KolmogorovSmirnov
    SHORTHAND = 'ks'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: cdf(d, n)
    # Returns the exact finite-n CDF Pr(D_n <= d).

    ##
    # :singleton-method: sf(d, n)
    # Returns the exact finite-n survival function Pr(D_n > d).

    ##
    # :singleton-method: quantile(p, n)
    # Returns a finite-n quantile satisfying the probability enclosure.
  end

  KS = KolmogorovSmirnov
end
