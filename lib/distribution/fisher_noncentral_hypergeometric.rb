require 'distribution/fisher_noncentral_hypergeometric/ruby'

module Distribution
  # Fisher's noncentral hypergeometric distribution. It describes the number
  # of items from the first of two fixed groups in a fixed-size sample, where
  # every item in the first group is +odds+ times as likely to be selected as
  # every item in the second group. The margins m1, m2 and n are fixed.
  module FisherNoncentralHypergeometric
    SHORTHAND = 'fnchg'
    extend Distributable

    create_distribution_methods

    ##
    # :singleton-method: pmf(k,m1,m2,n,odds)
    # Probability mass function for Fisher's noncentral hypergeometric
    # distribution.

    ##
    # :singleton-method: cdf(k,m1,m2,n,odds)
    # Cumulative distribution function.

    ##
    # :singleton-method: sf(k,m1,m2,n,odds)
    # Survival function P(X > k).

    ##
    # :singleton-method: quantile(p,m1,m2,n,odds)
    # Inverse cumulative distribution function.

    ##
    # :singleton-method: mean(m1,m2,n,odds)
    # Distribution mean.

    ##
    # :singleton-method: variance(m1,m2,n,odds)
    # Distribution variance.
  end
end
