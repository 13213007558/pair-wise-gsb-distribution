require 'distribution/compound_poisson/ruby'

module Distribution
  # Compound Poisson distribution: a Poisson number of independent requests,
  # each adding a non-negative integer cost from a caller-supplied PMF.
  module CompoundPoisson
    SHORTHAND = 'cmpois'
    extend Distributable
    create_distribution_methods
  end
end
