require 'distribution/poisson_binomial/ruby'

module Distribution
  # Probabilities for a sum of independent Bernoulli variables with
  # potentially different success probabilities.
  module PoissonBinomial
    SHORTHAND = 'pbino'

    extend Distributable
    create_distribution_methods
  end
end
