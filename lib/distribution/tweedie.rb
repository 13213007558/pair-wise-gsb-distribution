require 'distribution/tweedie/ruby'

module Distribution
  # Tweedie distribution with variance power between one and two.
  #
  # This member of the Tweedie family is a compound Poisson sum of Gamma
  # claim sizes. It has a probability mass at zero and a continuous density
  # on the positive real line.
  module Tweedie
    extend Distributable
    SHORTHAND = 'tweedie'
    create_distribution_methods

    ##
    # :singleton-method: pdf(x, mean, dispersion, power, tolerance = 1e-12)
    # Positive-part probability density. Negative +x+ and +x+ equal to zero
    # return zero; use zero_probability for the separate point mass.

    ##
    # :singleton-method: log_pdf(x, mean, dispersion, power, tolerance = 1e-12)
    # Logarithm of the positive-part probability density.

    ##
    # :singleton-method: zero_probability(mean, dispersion, power)
    # Probability mass at zero.
    #
    # :singleton-method: rng(mean, dispersion, power, seed = nil, random: nil)
    # Return a reproducible sampler from the compound Poisson-Gamma law.
  end
end
