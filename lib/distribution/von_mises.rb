require 'distribution/von_mises/ruby'

module Distribution
  # Circular normal distribution on angles measured in radians. Angles are
  # treated modulo 2*pi, so angles separated by a full turn are equivalent.
  module VonMises
    SHORTHAND = 'vmis'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pdf(theta, mu = 0, kappa = 1)
    # Circular probability density.

    ##
    # :singleton-method: cdf(x, mu = 0, kappa = 1, cut = 0)
    # CDF along [cut, cut + 2*pi].

    ##
    # :singleton-method: quantile(probability, mu = 0, kappa = 1, cut = 0)
    # Inverse CDF in the one-turn interval beginning at +cut+.

    ##
    # :singleton-method: arc_probability(start_angle, arc_length, mu = 0, kappa = 1)
    # Probability of a directed circular arc.
  end
end
