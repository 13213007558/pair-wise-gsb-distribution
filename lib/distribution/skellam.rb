require 'distribution/skellam/ruby'

module Distribution
  # Difference of two independent Poisson random variables.
  module Skellam
    SHORTHAND = 'skel'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pmf(k, mu1, mu2)
    # PMF for the Skellam distribution. +k+ is an integer and +mu1+
    # and +mu2+ are the non-negative means of the two Poisson variables.

    ##
    # :singleton-method: log_pmf(k, mu1, mu2)
    # Log-PMF for the Skellam distribution.

    ##
    # :singleton-method: cdf(k, mu1, mu2)
    # CDF for the Skellam distribution.

    ##
    # :singleton-method: survival(k, mu1, mu2)
    # Directly computed upper tail P(X > k).

    ##
    # :singleton-method: p_value(probability, mu1, mu2)
    # Smallest integer +k+ for which the CDF is at least +probability+.
  end
end
