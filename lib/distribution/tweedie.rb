require 'distribution/tweedie/ruby'

module Distribution
  # Tweedie distribution with variance power 1 < p < 2.
  #
  # This member of the Tweedie family is a compound Poisson-Gamma
  # distribution with an atom at zero and a continuous component for x > 0.
  module Tweedie
    SHORTHAND = 'tweedie'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pdf(x, mean, dispersion, power, tolerance: 1e-12)
    # Probability density on x > 0. Zero is handled by zero_mass/pmf.

    ##
    # :singleton-method: log_pdf(x, mean, dispersion, power, tolerance: 1e-12)
    # Log probability density on x > 0. Zero is handled by zero_mass/pmf.

    ##
    # :singleton-method: zero_mass(mean, dispersion, power)
    # Probability mass P(X = 0).

    ##
    # :singleton-method: rng(mean, dispersion, power, seed = nil, random: nil)
    # Returns a lambda which generates Tweedie variates.
  end
end
