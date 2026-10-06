require 'distribution/fisher_noncentral_hypergeometric/ruby'

module Distribution
  # Fisher's noncentral hypergeometric distribution.
  #
  # Models a sample of +n+ items drawn without replacement from a
  # population of +m1+ type-1 and +m2+ type-2 items, where the two types
  # have different chances of being sampled: each type-1 item is +odds+
  # times as likely to be sampled as a type-2 item. The number of drawn
  # type-1 items follows
  #
  #   P(X = k) = C(m1,k) * C(m2,n-k) * odds^k / Z,
  #
  # with Z the normalizing constant and k bounded by the fixed margins,
  # max(0, n-m2) <= k <= min(n, m1).
  #
  # This is Fisher's definition (a fixed-size weighted sample, equivalent
  # to a conditional binomial model), not Wallenius' definition based on
  # sequential weighted draws. For odds == 1 it reduces to the
  # hypergeometric distribution.
  module FisherNoncentralHypergeometric
    SHORTHAND = 'fnchypg'
    extend Distributable

    create_distribution_methods

    ##
    # :singleton-method: pdf(k,m1,m2,n,odds)
    # Probability of drawing +k+ type-1 items in a sample of +n+ from a
    # population of +m1+ type-1 and +m2+ type-2 items with sampling
    # +odds+ for type-1 relative to type-2.

    ##
    # :singleton-method: cdf(k,m1,m2,n,odds)
    # Cumulative distribution function P(X <= k).

    ##
    # :singleton-method: sf(k,m1,m2,n,odds)
    # Survival function P(X > k).

    ##
    # :singleton-method: quantile(pr,m1,m2,n,odds)
    # Smallest k such that cdf(k) >= pr.

    ##
    # :singleton-method: mean(m1,m2,n,odds)
    # Expected number of drawn type-1 items.

    ##
    # :singleton-method: variance(m1,m2,n,odds)
    # Variance of the number of drawn type-1 items.
  end
end
