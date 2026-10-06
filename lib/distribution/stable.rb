require 'distribution/stable/ruby'

module Distribution
  # Nolan's S0 four-parameter stable distribution.
  #
  # Parameters are +alpha+, +beta+, +delta+ and +gamma+. Public PDF, CDF and
  # RNG entry points always use the S0 parameterization.
  module Stable
    SHORTHAND = 'stab'
    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pdf(x, alpha:, beta: 0, delta: 0, gamma: 1, diagnostics: false, tolerance: 1e-10)
    # Returns the S0 stable probability density at +x+.

    ##
    # :singleton-method: cdf(x, alpha:, beta: 0, delta: 0, gamma: 1, diagnostics: false, tolerance: 1e-10)
    # Returns the S0 stable cumulative probability at +x+.

    ##
    # :singleton-method: rng(alpha:, beta: 0, delta: 0, gamma: 1, source: nil)
    # Returns a proc that emits S0 stable variates. +source+ may be a Random
    # instance or any zero-argument callable returning Uniform(0, 1) values.
  end
end
