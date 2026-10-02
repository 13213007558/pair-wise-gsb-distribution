module Distribution
  module CompoundPoisson
    module Ruby_
      # Tolerance used to validate that the single-cost probabilities
      # provided by the caller add up to one.
      PROBABILITY_SUM_TOLERANCE = 1e-9

      class << self
        # P(S = k), where S = X_1 + ... + X_N, N ~ Poisson(lambda_value)
        # and X_i are i.i.d. with the caller-supplied +cost_table+
        # ({cost => probability}).
        def pdf(k, lambda_value, cost_table)
          validate_point(k)
          full_pmf(lambda_value, cost_table, k)[k]
        end

        # P(S <= k).
        def cdf(k, lambda_value, cost_table)
          validate_point(k)
          full_pmf(lambda_value, cost_table, k).sum
        end

        # PMF of the total cost over a given range, together with
        # cumulative probabilities and the probability mass falling
        # outside the range.
        #
        # == Arguments
        #   * +lambda_value+ - mean of the Poisson request count (>= 0)
        #   * +cost_table+   - Hash of {cost => probability} with
        #                      non-negative integer costs, including 0
        #   * +range+        - an inclusive Range of integers, or a single
        #                      Integer taken as 0..max
        #
        # Returns a Distribution::CompoundPoisson::Result.
        def pmf_table(lambda_value, cost_table, range)
          range = 0..range if range.is_a?(Integer)
          validate_range(range)

          probs = full_pmf(lambda_value, cost_table, range.end)

          pmf = {}
          cdf = {}
          cumulative = probs[0...range.begin].sum
          lower_mass = cumulative
          (range.begin..range.end).each do |k|
            cumulative += probs[k]
            pmf[k] = probs[k]
            cdf[k] = cumulative
          end
          upper_mass = [1.0 - cumulative, 0.0].max

          Result.new(range, pmf, cdf, lower_mass, upper_mass)
        end

        private

        # Exact PMF of the compound Poisson sum on 0..max via the
        # Panjer recursion:
        #   g(0) = exp(-lambda * (1 - P(X = 0)))
        #   g(s) = (lambda / s) * sum_j j * P(X = j) * g(s - j)
        # The recursion is exact (up to floating point rounding); the
        # Poisson count is never truncated at an arbitrary bound.
        def full_pmf(lambda_value, cost_table, max)
          validate_lambda(lambda_value)
          costs = validate_cost_table(cost_table)

          probs = Array.new(max + 1, 0.0)
          p_zero = cost_table[0] || 0.0
          probs[0] = Math.exp(-lambda_value * (1.0 - p_zero))
          (1..max).each do |s|
            acc = 0.0
            costs.each do |cost, prob|
              break if cost > s
              acc += cost * prob * probs[s - cost]
            end
            probs[s] = lambda_value * acc / s
          end
          probs
        end

        def validate_lambda(lambda_value)
          unless lambda_value.is_a?(Numeric) && lambda_value >= 0
            raise ArgumentError, 'lambda must be a non-negative number'
          end
        end

        def validate_cost_table(cost_table)
          unless cost_table.is_a?(Hash) && !cost_table.empty?
            raise ArgumentError, 'cost_table must be a non-empty Hash of {cost => probability}'
          end
          sum = 0.0
          cost_table.each do |cost, prob|
            unless cost.is_a?(Integer) && cost >= 0
              raise ArgumentError, 'costs must be non-negative integers'
            end
            unless prob.is_a?(Numeric) && prob >= 0
              raise ArgumentError, 'probabilities must be non-negative numbers'
            end
            sum += prob
          end
          unless (sum - 1.0).abs <= PROBABILITY_SUM_TOLERANCE
            raise ArgumentError, "probabilities must sum to 1 (got #{sum})"
          end
          cost_table.reject { |cost, prob| cost.zero? || prob.zero? }.sort
        end

        def validate_point(k)
          unless k.is_a?(Integer) && k >= 0
            raise ArgumentError, 'k must be a non-negative integer'
          end
        end

        def validate_range(range)
          valid = range.is_a?(Range) && !range.exclude_end? &&
                  range.begin.is_a?(Integer) && range.begin >= 0 &&
                  range.end.is_a?(Integer) && range.begin <= range.end
          raise ArgumentError, 'range must be an inclusive integer range with 0 <= min <= max' unless valid
        end
      end
    end
  end
end
