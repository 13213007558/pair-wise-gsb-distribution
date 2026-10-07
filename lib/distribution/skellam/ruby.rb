module Distribution
  module Skellam
    module Ruby_
      class << self
        def pmf(k, mu1, mu2)
          Math.exp(log_pmf(k, mu1, mu2))
        end
        alias_method :pdf, :pmf

        def log_pmf(k, mu1, mu2)
          mu1 = validate_mean(mu1, :mu1)
          mu2 = validate_mean(mu2, :mu2)
          return -Float::INFINITY unless k.is_a?(Integer)
          return log_pmf(-k, mu2, mu1) if mu1 < mu2
          return log_poisson_pmf(k, mu1) if mu2.zero?

          log_convolution_pmf(k, mu1, mu2)
        end
        alias_method :log_pdf, :log_pmf

        def cdf(k, mu1, mu2)
          mu1 = validate_mean(mu1, :mu1)
          mu2 = validate_mean(mu2, :mu2)
          k = floor_integer(k)

          if mu1 >= mu2
            oriented_cdf(k, mu1, mu2)
          else
            oriented_survival(-k - 1, mu2, mu1)
          end
        end

        def survival(k, mu1, mu2)
          mu1 = validate_mean(mu1, :mu1)
          mu2 = validate_mean(mu2, :mu2)
          k = floor_integer(k)

          if mu1 >= mu2
            oriented_survival(k, mu1, mu2)
          else
            oriented_cdf(-k - 1, mu2, mu1)
          end
        end
        alias_method :sf, :survival

        def quantile(probability, mu1, mu2)
          mu1 = validate_mean(mu1, :mu1)
          mu2 = validate_mean(mu2, :mu2)
          unless probability.is_a?(Numeric) && probability.to_f.nan? == false
            raise ArgumentError, 'probability must be a number'
          end

          probability = probability.to_f
          if probability < 0.0 || probability > 1.0
            raise ArgumentError, 'probability must be between 0 and 1'
          end
          return mu2.zero? ? 0 : -Float::INFINITY if probability.zero?
          return mu1.zero? ? 0 : Float::INFINITY if probability == 1.0

          center = mu2.zero? ? 0 : Integer((mu1 - mu2).floor)

          candidate = center
          lower = nil
          upper = nil
          if cdf(candidate, mu1, mu2) >= probability
            distance = 1
            loop do
              next_candidate = candidate - distance
              next_candidate = -1 if mu2.zero? && next_candidate < -1
              if cdf(next_candidate, mu1, mu2) < probability
                lower = next_candidate
                upper = candidate
                break
              end
              candidate = next_candidate
              distance *= 2
            end
          else
            distance = 1
            loop do
              next_candidate = candidate + distance
              if cdf(next_candidate, mu1, mu2) >= probability
                lower = candidate
                upper = next_candidate
                break
              end
              candidate = next_candidate
              distance *= 2
            end
          end

          while lower + 1 < upper
            middle = (lower + upper) / 2
            if cdf(middle, mu1, mu2) >= probability
              upper = middle
            else
              lower = middle
            end
          end

          upper
        end
        alias_method :p_value, :quantile

        private

        def validate_mean(value, name)
          unless value.is_a?(Numeric) && value >= 0
            raise ArgumentError, "#{name} must be a non-negative number"
          end
          value.to_f
        end

        def floor_integer(value)
          return value if value.is_a?(Integer)
          float_value = value.to_f
          unless float_value.finite?
            raise ArgumentError, 'k must be a finite number'
          end
          Integer(float_value.floor)
        end

        def log_poisson_pmf(k, mean)
          return -Float::INFINITY if k < 0
          return k.zero? ? 0.0 : -Float::INFINITY if mean.zero?
          k * Math.log(mean) - mean - Math.lgamma(k + 1).first
        end

        def poisson_pmf(k, mean)
          Math.exp(log_poisson_pmf(k, mean))
        end

        def poisson_cdf(k, mean)
          return 0.0 if k < 0
          return 1.0 if mean.zero?
          direct_poisson_cdf(k, mean)
        end

        def poisson_survival(k, mean)
          return 1.0 if k < 0
          return 0.0 if mean.zero?
          if k >= mean
            clamp_probability(direct_poisson_tail(k + 1, mean))
          else
            clamp_probability(1.0 - direct_poisson_cdf(k, mean))
          end
        end

        def clamp_probability(value)
          return 0.0 if value < 0.0
          return 1.0 if value > 1.0
          value
        end

        def direct_poisson_cdf(k, mean)
          if k >= mean
            clamp_probability(1.0 - direct_poisson_tail(k + 1, mean))
          else
            sum = 0.0
            probability = poisson_pmf(0, mean)
            0.upto(k) do |index|
              probability = probability * mean / index unless index.zero?
              sum += probability
            end
            clamp_probability(sum)
          end
        end

        def direct_poisson_tail(first, mean)
          sum = probability = poisson_pmf(first, mean)
          index = first + 1
          loop do
            probability *= mean / index
            sum += probability
            break if probability <= 1e-17 * sum
            index += 1
          end
          sum
        end

        def oriented_cdf(k, high_mean, low_mean)
          sum = 0.0
          each_significant_poisson(low_mean) do |n, probability|
            sum += probability * poisson_cdf(n + k, high_mean)
          end
          clamp_probability(sum)
        end

        def oriented_survival(k, high_mean, low_mean)
          sum = 0.0
          each_significant_poisson(low_mean) do |n, probability|
            sum += probability * poisson_survival(k + n, high_mean)
          end
          clamp_probability(sum)
        end

        def each_significant_poisson(mean)
          if mean.zero?
            yield 0, 1.0
            return
          end

          upper = lower = Integer(mean.floor)
          probability = poisson_pmf(upper, mean)
          lower_probability = probability
          yield upper, probability

          loop do
            upper_bound = poisson_upper_tail_bound(upper, probability, mean)
            lower_bound = poisson_lower_tail_bound(lower, lower_probability, mean)
            break if upper_bound + lower_bound <= 1e-15

            if upper_bound >= lower_bound
              upper += 1
              probability = poisson_pmf(upper, mean)
              yield upper, probability
            else
              lower -= 1
              lower_probability = poisson_pmf(lower, mean)
              yield lower, lower_probability
            end
          end
        end

        def poisson_upper_tail_bound(n, probability, mean)
          first = probability * mean / (n + 1)
          ratio = mean / (n + 2)
          ratio >= 1.0 ? first : first / (1.0 - ratio)
        end

        def poisson_lower_tail_bound(n, probability, mean)
          return 0.0 if n.zero?
          first = probability * n / mean
          ratio = n > 1 ? (n - 1) / mean : 0.0
          ratio >= 1.0 ? first : first / (1.0 - ratio)
        end

        def log_convolution_pmf(k, high_mean, low_mean)
          lower_bound = [0, -k].max
          discriminant = k * k + 4.0 * low_mean * high_mean
          mode = Integer(((Math.sqrt(discriminant) - k) / 2.0).floor)
          upper = lower = [mode, lower_bound].max

          log_sum = convolution_log_term(upper, k, low_mean, high_mean)

          loop do
            upper_bound = convolution_upper_log_bound(upper, k, low_mean, high_mean)
            lower_bound = convolution_lower_log_bound(lower, k, low_mean, high_mean)
            bound = log_max(upper_bound, lower_bound)
            break if bound == -Float::INFINITY || bound <= log_sum + Math.log(1e-15)

            if upper_bound == bound
              upper += 1
              log_sum = log_add(log_sum, convolution_log_term(upper, k, low_mean, high_mean))
            else
              lower -= 1
              log_sum = log_add(log_sum, convolution_log_term(lower, k, low_mean, high_mean))
            end
          end

          log_sum
        end

        def convolution_log_term(n, k, low_mean, high_mean)
          return -Float::INFINITY if n < 0 || n + k < 0
          log_poisson_pmf(n, low_mean) + log_poisson_pmf(n + k, high_mean)
        end

        def convolution_upper_log_bound(n, k, low_mean, high_mean)
          next_n = n + 1
          log_term = convolution_log_term(next_n, k, low_mean, high_mean)
          ratio = low_mean / next_n * high_mean / (next_n + k)
          return log_term if ratio >= 1.0
          log_term - Math.log(1.0 - ratio)
        end

        def convolution_lower_log_bound(n, k, low_mean, high_mean)
          return -Float::INFINITY if n.zero? || n + k <= 0
          previous_n = n - 1
          log_term = convolution_log_term(previous_n, k, low_mean, high_mean)
          ratio = previous_n / low_mean * (previous_n + k) / high_mean
          return log_term if ratio >= 1.0
          log_term - Math.log(1.0 - ratio)
        end

        def log_add(first, second)
          return second if first == -Float::INFINITY
          return first if second == -Float::INFINITY
          if first > second
            first + Math.log(1.0 + Math.exp(second - first))
          else
            second + Math.log(1.0 + Math.exp(first - second))
          end
        end

        def log_max(first, second)
          [first, second].max
        end
      end
    end
  end
end
