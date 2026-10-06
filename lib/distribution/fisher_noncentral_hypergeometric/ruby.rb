module Distribution
  module FisherNoncentralHypergeometric
    module Ruby_
      class << self
        def pmf(k, m1, m2, n, odds)
          Calculator.new(m1, m2, n, odds).pmf(k)
        end

        alias_method :pdf, :pmf

        def cdf(k, m1, m2, n, odds)
          Calculator.new(m1, m2, n, odds).cdf(k)
        end

        def sf(k, m1, m2, n, odds)
          Calculator.new(m1, m2, n, odds).sf(k)
        end

        def quantile(pr, m1, m2, n, odds)
          Calculator.new(m1, m2, n, odds).quantile(pr)
        end

        alias_method :p_value, :quantile

        def mean(m1, m2, n, odds)
          Calculator.new(m1, m2, n, odds).mean
        end

        def variance(m1, m2, n, odds)
          Calculator.new(m1, m2, n, odds).variance
        end
      end

      class Calculator
        NEGATIVE_INFINITY = -Float::INFINITY

        attr_reader :lower, :upper, :mode

        def initialize(m1, m2, n, odds)
          validate_parameters(m1, m2, n, odds)

          @m1 = m1
          @m2 = m2
          @n = n
          @lower = [0, n - m2].max
          @upper = [n, m1].min
          @log_odds = log_odds(odds)
          @mode = find_mode
          @log_weights = build_log_weights
          @log_partition = log_sum(@log_weights)
          @probabilities = @log_weights.map { |weight| Math.exp(weight - @log_partition) }
        end

        def pmf(k)
          validate_k(k)
          return 0.0 unless support?(k)

          @probabilities[index(k)]
        end

        def cdf(k)
          validate_point(k)
          point = k.floor
          return 0.0 if point < @lower
          return 1.0 if point >= @upper

          Math.exp(log_sum(@log_weights[0..index(point)]) - @log_partition)
        end

        def sf(k)
          validate_point(k)
          point = k.floor
          return 1.0 if point < @lower
          return 0.0 if point >= @upper

          Math.exp(log_sum(@log_weights[(index(point) + 1)..-1]) - @log_partition)
        end

        def quantile(pr)
          unless pr.is_a?(Numeric) && pr >= 0 && pr <= 1
            fail ArgumentError, 'p must be between 0 and 1'
          end

          return @lower if pr == 0
          return @upper if pr == 1

          target = pr.to_f
          if target <= 0.5
            cumulative = 0.0
            (@lower..@upper).each do |k|
              cumulative += @probabilities[index(k)]
              return k if cumulative >= target
            end
            @upper
          else
            upper_tail = 0.0
            tail_target = 1.0 - target
            @upper.downto(@lower) do |k|
              upper_tail += @probabilities[index(k)]
              return k if upper_tail >= tail_target
            end
            @lower
          end
        end

        def mean
          centered_first_moment
          @mode + @first_moment
        end

        def variance
          centered_first_moment
          [0.0, @second_moment - @first_moment * @first_moment].max
        end

        private

        def validate_parameters(m1, m2, n, odds)
          unless m1.is_a?(Integer) && m1 >= 0
            fail ArgumentError, 'm1 must be a non-negative integer'
          end
          unless m2.is_a?(Integer) && m2 >= 0
            fail ArgumentError, 'm2 must be a non-negative integer'
          end
          unless n.is_a?(Integer) && n >= 0 && n <= m1 + m2
            fail ArgumentError, 'n must be an integer between 0 and m1 + m2'
          end
          unless odds.is_a?(Numeric) && odds >= 0
            fail ArgumentError, 'odds must be non-negative'
          end
        end

        def validate_k(k)
          unless k.is_a?(Integer)
            fail ArgumentError, 'k must be an integer'
          end
        end

        def validate_point(k)
          unless k.is_a?(Numeric)
            fail ArgumentError, 'k must be numeric'
          end
        end

        def support?(k)
          k >= @lower && k <= @upper
        end

        def index(k)
          k - @lower
        end

        def log_odds(odds)
          return -Float::INFINITY if odds == 0
          return Float::INFINITY if odds.respond_to?(:infinite?) && odds.infinite?

          Math.log(odds.to_f)
        end

        # Logarithm of w(k + 1) / w(k), where w(k) is the unnormalized
        # Fisher weight C(m1,k) * C(m2,n-k) * odds**k.
        def log_weight_ratio(k)
          return @log_odds unless @log_odds.finite?

          @log_odds +
            Math.log(@m1 - k) +
            Math.log(@n - k) -
            Math.log(k + 1) -
            Math.log(@m2 - @n + k + 1)
        end

        def find_mode
          low = @lower
          high = @upper

          while low < high
            candidate = (low + high + 1) / 2
            if log_weight_ratio(candidate - 1) >= 0
              low = candidate
            else
              high = candidate - 1
            end
          end

          low
        end

        def build_log_weights
          weights = Array.new(@upper - @lower + 1, NEGATIVE_INFINITY)
          weights[index(@mode)] = 0.0

          (@mode - 1).downto(@lower) do |k|
            weights[index(k)] = weights[index(k + 1)] - log_weight_ratio(k)
          end

          (@mode + 1).upto(@upper) do |k|
            weights[index(k)] = weights[index(k - 1)] + log_weight_ratio(k - 1)
          end

          weights
        end

        def log_sum(values)
          max = values.max
          return NEGATIVE_INFINITY if max == NEGATIVE_INFINITY

          sum = values.reduce(0.0) do |total, value|
            value == NEGATIVE_INFINITY ? total : total + Math.exp(value - max)
          end

          max + Math.log(sum)
        end

        def centered_first_moment
          return if defined?(@first_moment)

          @first_moment = 0.0
          @second_moment = 0.0

          @probabilities.each_with_index do |probability, offset|
            centered = @lower + offset - @mode
            @first_moment += centered * probability
            @second_moment += centered * centered * probability
          end
        end
      end
    end
  end
end
