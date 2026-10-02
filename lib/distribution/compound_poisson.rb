require 'distribution/compound_poisson/ruby'

module Distribution
  # Compound Poisson distribution over a non-negative integer total cost.
  #
  # The number of requests N follows a Poisson distribution with mean
  # +lambda+, and each request carries an independent, discrete,
  # non-negative integer cost whose probability table is supplied by the
  # caller. The distribution of the total cost
  # S = X_1 + ... + X_N is computed exactly with the Panjer recursion,
  # so the Poisson count is never truncated at an arbitrary bound.
  #
  # == Usage
  #   table = {0 => 0.2, 1 => 0.5, 3 => 0.3}
  #   res = Distribution::CompoundPoisson.pmf_table(2.5, table, 0..10)
  #   res.pmf          # => {0 => ..., 1 => ..., ..., 10 => ...}
  #   res.cdf          # => cumulative P(S <= k) for k in 0..10
  #   res.lower_mass   # => P(S < 0)
  #   res.upper_mass   # => P(S > 10)
  #   res.omitted_mass # => mass outside the requested range
  module CompoundPoisson
    SHORTHAND = 'cpoi'

    # Result of CompoundPoisson.pmf_table: PMF and CDF over the
    # requested range plus the probability mass omitted on each side.
    class Result
      include Enumerable

      # Range of totals covered by +pmf+ and +cdf+.
      attr_reader :range
      # Hash of {total => P(S = total)} for totals in +range+.
      attr_reader :pmf
      # Hash of {total => P(S <= total)} for totals in +range+.
      attr_reader :cdf
      # P(S < range.begin).
      attr_reader :lower_mass
      # P(S > range.end).
      attr_reader :upper_mass

      def initialize(range, pmf, cdf, lower_mass, upper_mass) # :nodoc:
        @range = range
        @pmf = pmf
        @cdf = cdf
        @lower_mass = lower_mass
        @upper_mass = upper_mass
      end

      # Total probability mass outside +range+.
      def omitted_mass
        @lower_mass + @upper_mass
      end

      # Total probability mass inside +range+.
      def covered_mass
        @cdf[@range.end] - @lower_mass
      end

      def [](k)
        @pmf[k]
      end

      def each(&block)
        @pmf.each(&block)
      end
    end

    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pdf(k, lambda_value, cost_table)
    # P(S = k) for the compound Poisson total cost S.

    ##
    # :singleton-method: cdf(k, lambda_value, cost_table)
    # P(S <= k) for the compound Poisson total cost S.

    ##
    # :singleton-method: pmf_table(lambda_value, cost_table, range)
    # PMF of the total cost over +range+ (an inclusive Range, or an
    # Integer taken as 0..max), with cumulative probabilities and
    # omitted-mass information. Returns a CompoundPoisson::Result.
  end
end
