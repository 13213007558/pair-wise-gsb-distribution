require 'distribution/skellam/ruby'

module Distribution
  # Difference of two independent Poisson random variables.
  module Skellam
    SHORTHAND = 'skel'

    extend Distributable
    create_distribution_methods
  end
end
