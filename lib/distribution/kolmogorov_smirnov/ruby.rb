require 'bigdecimal'

module Distribution
  module KolmogorovSmirnov
    module Ruby_
      class << self
        def cdf(d, n)
          exact_cdf(d, n)
        end

        def sf(d, n)
          exact_sf(d, n)
        end

        def quantile(p, n)
          n = validate_n(n)
          validate_probability(p)
          return lower_support(n) if p.zero?
          return 1.0 if p == 1.0

          candidate = if p <= 0.5
            inverse_with_cdf(p, n)
          else
            inverse_with_sf(1.0 - p, n)
          end
          [[lower_support(n), candidate.next_float].max, 1.0].min
        end
        alias p_value quantile

        private

        def exact_cdf(d, n)
          n = validate_n(n)
          d = d.to_f
          return Float::NAN if d.nan?
          return 0.0 if d <= lower_support(n)
          return 1.0 if d >= 1.0

          t = n * d
          return lower_boundary_cdf(t, n) if t <= 1.0

          if t >= n - 1
            survival = 2.0 * Math.exp(n * Math.log1p(-d))
            return 1.0 - survival
          end

          dmtw_cdf(n, d)
        end

        def exact_sf(d, n)
          n = validate_n(n)
          d = d.to_f
          return Float::NAN if d.nan?
          return 1.0 if d <= lower_support(n)
          return 0.0 if d >= 1.0

          t = n * d
          if t <= 1.0
            return 1.0 - lower_boundary_cdf(t, n)
          end

          if t >= n - 1
            return 2.0 * Math.exp(n * Math.log1p(-d))
          end

          if n <= 100
            high_precision_sf(d, n)
          else
            1.0 - dmtw_cdf(n, d)
          end
        end

        def lower_support(n)
          0.5 / n
        end

        def validate_n(n)
          value = Integer(n)
          if value <= 0 || value.to_f != n.to_f
            fail ArgumentError, 'n must be a positive integer'
          end
          value
        rescue ArgumentError, TypeError
          raise ArgumentError, 'n must be a positive integer'
        end

        def validate_probability(p)
          unless p.is_a?(Numeric) && p >= 0.0 && p <= 1.0
            fail RangeError, 'probability must be in [0, 1]'
          end
        end

        def lower_boundary_cdf(t, n)
          factor = n <= 140 ? factorial_ratio(n) : Math.exp(log_factorial_ratio(n))
          factor * (2.0 * t - 1.0)**n
        end

        def factorial_ratio(n)
          (1..n).reduce(1.0) { |product, i| product * i / n }
        end

        def log_factorial_ratio(n)
          Math::Log.lgamma(n + 1).first - n * Math.log(n)
        end

        def dmtw_cdf(n, d)
          probability = dmtw_probability(n, d).to_f
          [[probability, 0.0].max, 1.0].min
        end

        def high_precision_sf(d, n)
          BigDecimal.limit(50)
          key = [n, d]
          @high_precision_sf_cache ||= {}
          return @high_precision_sf_cache[key] if @high_precision_sf_cache.key?(key)

          probability = 1 - dmtw_probability(n, BigDecimal(d.to_s))
          probability = 0 if probability < 0
          @high_precision_sf_cache[key] = probability.to_f
        end

        def dmtw_probability(n, d)
          @probability_cache ||= {}
          key = [n, d]
          return @probability_cache[key] if @probability_cache.key?(key)

          decimal = d.is_a?(BigDecimal)
          t = n * d
          k = t.ceil
          h = k - t
          m = 2 * k - 1
          matrix = build_dmtw_matrix(m, h, decimal)

          power = identity_matrix(m, decimal)
          remaining = n
          exponent = 0
          matrix_exponent = 0

          while remaining > 0
            if remaining.odd?
              power = multiply_matrices(power, matrix)
              exponent += matrix_exponent
            end
            matrix = multiply_matrices(matrix, matrix)
            unless decimal
              if matrix[k - 1][k - 1].abs > 2.0**128
                scale_matrix(matrix, 2.0**-128)
                matrix_exponent = matrix_exponent * 2 + 128
              else
                matrix_exponent *= 2
              end
            end
            remaining >>= 1
          end

          p = power[k - 1][k - 1]
          (1..n).each { |i| p = i * p / n }
          unless decimal
            p = Math.ldexp(p, exponent)
          end
          @probability_cache[key] = p
        end

        def build_dmtw_matrix(m, h, decimal = false)
          zero = decimal ? BigDecimal('0') : 0.0
          one = decimal ? BigDecimal('1') : 1.0
          matrix = Array.new(m) { Array.new(m, zero.dup) }
          v = Array.new(m)
          w = Array.new(m)
          factorial = one.dup

          (1..m).each do |j|
            w[j - 1] = factorial
            factorial /= j
            v[j - 1] = (one - h**j) * factorial
          end

          tail = [2 * h - one, zero].max**m - 2 * h**m
          v[-1] = (one + tail) * factorial

          1.upto(m - 1) do |column|
            column.upto(m) do |row|
              matrix[row - 1][column] = w[row - column]
            end
          end

          m.times do |i|
            matrix[i][0] = v[i]
            matrix[-1][i] = v[m - 1 - i]
          end

          matrix
        end

        def inverse_with_cdf(p, n)
          low = lower_support(n)
          high = 1.0
          64.times do
            middle = (low + high) / 2.0
            if cdf(middle, n) <= p
              low = middle
            else
              high = middle
            end
          end
          high
        end

        def inverse_with_sf(q, n)
          low = lower_support(n)
          high = 1.0
          64.times do
            middle = (low + high) / 2.0
            if sf(middle, n) >= q
              low = middle
            else
              high = middle
            end
          end
          low
        end

        def identity_matrix(size, decimal = false)
          one = decimal ? BigDecimal('1') : 1.0
          zero = decimal ? BigDecimal('0') : 0.0
          Array.new(size) do |row|
            Array.new(size) { |column| row == column ? one.dup : zero.dup }
          end
        end

        def scale_matrix(matrix, scale)
          matrix.each { |row| row.map! { |value| value * scale } }
        end

        def multiply_matrices(left, right)
          size = left.size
          result = Array.new(size) { Array.new(size, 0.0) }
          size.times do |row|
            size.times do |k|
              next if left[row][k].zero?

              value = left[row][k]
              size.times do |column|
                result[row][column] += value * right[k][column]
              end
            end
          end
          result
        end

      end
    end
  end
end
