require 'distribution/noncentral_f/ruby'
require 'distribution/f'

module Distribution
  # Noncentral F distribution with numerator degrees of freedom +d1+,
  # denominator degrees of freedom +d2+, and noncentrality parameter +ncp+.
  module NoncentralF
    SHORTHAND = 'ncfdist'
    extend Distributable
    create_distribution_methods
  end
end
