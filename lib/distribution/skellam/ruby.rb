module Distribution
  module Skellam
    module Ruby_
      class << self
        def log_pmf(k, mean1, mean2)
          mean1, mean2 = validate_means(mean1, mean2)
          return -Float::INFINITY unless integer_value?(k)

          k = k.to_i

          return log_poisson_pmf(k, mean1) if mean2.zero?
          return log_poisson_pmf(-k, mean2) if mean1.zero?

          order = k.abs
          z = 2.0 * Math.sqrt(mean1 * mean2)
          -mean1 - mean2 +
            (k * 0.5 * (Math.log(mean1) - Math.log(mean2))) +
            log_modified_bessel(order, z)
        end

        def pmf(k, mean1, mean2)
          Math.exp(log_pmf(k, mean1, mean2))
        end

        alias_method :pdf, :pmf
        alias_method :log_pdf, :log_pmf

        def cdf(k, mean1, mean2)
          mean1, mean2 = validate_means(mean1, mean2)
          k = k.floor

          return k < 0 ? 0.0 : 1.0 if mean1.zero? && mean2.zero?
          return poisson_cdf(k, mean1) if mean2.zero?
          if mean1.zero?
            return 1.0 if k >= 0
            Math.gammq(-k, mean2)
          else
            weighted_poisson_sum(mean1, mean2, k, :cdf)
          end
        end

        def survival(k, mean1, mean2)
          mean1, mean2 = validate_means(mean1, mean2)
          k = k.floor

          return k < 0 ? 1.0 : 0.0 if mean1.zero? && mean2.zero?
          return poisson_survival(k, mean1) if mean2.zero?
          if mean1.zero?
            return 0.0 if k >= 0
            Math.gammp(-k, mean2)
          else
            weighted_poisson_sum(mean1, mean2, k, :survival)
          end
        end

        alias_method :sf, :survival

        def quantile(probability, mean1, mean2)
          mean1, mean2 = validate_means(mean1, mean2)
          fail ArgumentError, 'probability must be between 0 and 1' if probability < 0 || probability > 1
          return -Float::INFINITY if probability.zero?
          return Float::INFINITY if probability == 1
          return 0 if mean1.zero? && mean2.zero?

          center = mean1 - mean2
          spread = Math.sqrt(mean1 + mean2)
          z = Distribution::Normal.p_value(probability)
          point = (center + z * spread).floor
          low = point - 4
          high = point + 4

          low -= 8 while cdf(low, mean1, mean2) >= probability
          high += 8 while cdf(high, mean1, mean2) < probability

          while low + 1 < high
            middle = (low + high).div(2)
            if cdf(middle, mean1, mean2) >= probability
              high = middle
            else
              low = middle
            end
          end

          high
        end

        alias_method :p_value, :quantile

        private

        def validate_means(mean1, mean2)
          [mean1, mean2].each do |mean|
            unless mean.is_a?(Numeric) && mean.finite? && mean >= 0
              fail ArgumentError, 'Poisson means must be finite non-negative numbers'
            end
          end
          [mean1.to_f, mean2.to_f]
        end

        def integer_value?(value)
          value.is_a?(Numeric) && !value.is_a?(Complex) && value == value.to_i
        end

        def log_poisson_pmf(k, mean)
          return -Float::INFINITY if k < 0
          return k.zero? ? 0.0 : -Float::INFINITY if mean.zero?
          -mean + k * Math.log(mean) - Math.lgamma(k + 1).first
        end

        def poisson_cdf(k, mean)
          return 0.0 if k < 0
          return 1.0 if k >= summation_limit(mean)
          return 0.0 if k > mean && log_poisson_pmf(k, mean) < -80.0
          discrete_poisson_cdf(k, mean)
        end

        def poisson_survival(k, mean)
          lower = k + 1
          poisson_probability_at_least(lower, mean)
        end

        def poisson_probability_at_least(lower, mean)
          return 1.0 if lower <= 0
          return 0.0 if lower >= summation_limit(mean)
          return 0.0 if lower > mean && log_poisson_pmf(lower, mean) < -80.0
          discrete_poisson_survival(lower, mean)
        end

        def weighted_poisson_sum(mean1, mean2, k, kind)
          if mean2 <= mean1
            summation_mean = mean2
            lower = 0
            summation_factor = if kind == :cdf
                                 proc { |j| poisson_cdf(k + j, mean1) }
                               else
                                 proc { |j| poisson_survival(k + j, mean1) }
                                end
          else
            summation_mean = mean1
            lower = kind == :cdf ? 0 : [0, k + 1].max
            summation_factor = if kind == :cdf
                                 proc { |i| poisson_probability_at_least(i - k, mean2) }
                               else
                                 proc { |i| poisson_cdf(i - k - 1, mean2) }
                               end
          end

          upper = [summation_limit(summation_mean), lower].max
          if kind == :cdf && mean2 > mean1
            required_upper = [0, -k].max
            upper = [upper, required_upper].max if log_poisson_pmf(required_upper, summation_mean) > -700.0
          end
          maximum_log_term = -Float::INFINITY
          log_terms = []

          lower.upto(upper) do |count|
            factor = summation_factor.call(count)
            next if factor <= 0
            log_weight = log_poisson_pmf(count, summation_mean)
            next if log_weight == -Float::INFINITY
            log_term = log_weight + Math.log(factor)
            log_terms << log_term
            maximum_log_term = log_term if log_term > maximum_log_term
          end

          return 0.0 if maximum_log_term == -Float::INFINITY

          sum = log_terms.inject(0.0) do |total, log_term|
            total + Math.exp(log_term - maximum_log_term)
          end
          Math.exp(maximum_log_term + Math.log(sum))
        end

        def summation_limit(mean)
          return 0 if mean.zero?
          (mean + 40.0 * Math.sqrt(mean) + 50.0).ceil
        end

        def discrete_poisson_cdf(k, mean)
          return 0.0 if k > mean && log_poisson_pmf(k, mean) < -740.0
          log_terms = 0.upto(k).map { |count| log_poisson_pmf(count, mean) }
          Math.exp(log_sum(log_terms))
        end

        def discrete_poisson_survival(lower, mean)
          first_log_term = log_poisson_pmf(lower, mean)
          return 0.0 if first_log_term < -80.0
          log_terms = [first_log_term]
          log_term = first_log_term
          count = lower

          maximum_log_term = first_log_term
          while count < summation_limit(mean) &&
                log_term > [maximum_log_term - 80.0, -700.0].max
            count += 1
            log_term += Math.log(mean) - Math.log(count)
            log_terms << log_term
            maximum_log_term = log_term if log_term > maximum_log_term
          end

          Math.exp(log_sum(log_terms))
        end

        def log_modified_bessel(order, z)
          return 0.0 if z.zero? && order.zero?
          return -Float::INFINITY if z.zero?

          log_zero = if z <= 100.0
                       log_modified_bessel_zero_series(z)
                     else
                       log_modified_bessel_zero_asymptotic(z)
                     end
          log_zero + log_modified_bessel_ratios(order, z)
        end

        def log_modified_bessel_zero_series(z)
          x = z * z / 4.0
          terms = [0.0]
          1.upto(summation_limit(z)) do |j|
            terms << j * Math.log(x) - 2.0 * Math.lgamma(j + 1).first
          end
          log_sum(terms)
        end

        def log_modified_bessel_zero_asymptotic(z)
          factor = 1.0 + 1.0 / (8.0 * z) +
            9.0 / (128.0 * z * z) +
            225.0 / (3072.0 * z * z * z) +
            11025.0 / (98304.0 * z**4) +
            893025.0 / (2_359_296.0 * z**5)
          z - Math.log(2.0 * Math::PI * z) / 2.0 + Math.log(factor)
        end

        def log_modified_bessel_ratios(order, z)
          start = [order + 50, (z + 20.0 * Math.sqrt(z) + 50.0).ceil].max
          ratio = 0.0
          total = 0.0

          start.downto(1) do |j|
            denominator = (2.0 * j) / z + ratio
            ratio = denominator.zero? ? 0.0 : 1.0 / denominator
            total += Math.log(ratio) if j <= order && ratio > 0
          end

          total
        end

        def log_sum(log_terms)
          maximum = log_terms.max
          sum = log_terms.inject(0.0) do |total, log_term|
            total + Math.exp(log_term - maximum)
          end
          maximum + Math.log(sum)
        end
      end
    end
  end
end
