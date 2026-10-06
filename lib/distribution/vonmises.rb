require 'distribution/vonmises/ruby'

module Distribution
  # From Wikipedia:
  #   The von Mises distribution (also known as the circular normal
  #   distribution) is a continuous probability distribution on the
  #   circle. It is a close approximation to the wrapped normal
  #   distribution and is the circular analogue of the normal
  #   distribution, used for directional data: angles close to each
  #   other across the 0/2pi boundary (like 359 deg and 1 deg) are
  #   treated as neighbours, not as distant values.
  #
  # All angles are in radians. The mean direction +mu+ accepts any
  # real value (it may exceed one period); the concentration +kappa+
  # ranges from 0 (circular uniform) to arbitrarily large values.
  module Vonmises
    SHORTHAND = 'vonmises'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pdf(x, mu, kappa)
    # Returns the von Mises PDF at angle +x+ (radians)

    ##
    # :singleton-method: cdf(x, mu, kappa, cut)
    # Returns the von Mises CDF anchored at +cut+, monotone on
    # [cut, cut + 2*PI]

    ##
    # :singleton-method: quantile(qn, mu, kappa, cut)
    # Returns the inverse CDF on the period starting at +cut+

    ##
    # :singleton-method: arc_probability(start, length, mu, kappa)
    # Returns the probability of the directed arc starting at +start+
    # with signed +length+ (radians)

    ##
    # :singleton-method: rng(mu, kappa, seed)
    # Returns a lambda which returns a random angle from the
    # von Mises distribution
  end
end
