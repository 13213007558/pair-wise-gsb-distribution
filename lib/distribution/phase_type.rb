require 'distribution/phase_type/ruby'

module Distribution
  # Continuous phase-type distribution represented by an initial transient
  # probability vector and a continuous-time sub-generator.
  module PhaseType
    SHORTHAND = 'phtype'
    extend Distributable

    class Error < StandardError; end
    class ConvergenceError < RuntimeError; end

    def self.new(alpha, t, **options)
      Ruby_::Distribution.new(alpha, t, **options)
    end
  end
end
