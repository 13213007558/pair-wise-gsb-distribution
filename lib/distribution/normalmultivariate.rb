module Distribution
  # Monte Carlo rectangle probabilities for finite-dimensional normal distributions.
  module NormalMultivariate
    INVERSE_NORMAL_A = [-3.969683028665376e+01, 2.209460984245205e+02,
                        -2.759285104469687e+02, 1.383577518672690e+02,
                        -3.066479806614716e+01, 2.506628277459239e+00].freeze
    INVERSE_NORMAL_B = [-5.447609879822406e+01, 1.615858368580409e+02,
                        -1.556989798598866e+02, 6.680131188771972e+01,
                        -1.328068155288572e+01].freeze
    INVERSE_NORMAL_C = [-7.784894002430293e-03, -3.223964580411365e-01,
                        -2.400758277161838e+00, -2.549732539343734e+00,
                        4.374664141464968e+00, 2.938163982698783e+00].freeze
    INVERSE_NORMAL_D = [7.784695709041462e-03, 3.224671290700398e-01,
                        2.445134137142996e+00, 3.754408661907416e+00].freeze

    class << self
      # Estimates P(lower_i <= X_i <= upper_i for every i), where X is normal
      # with the supplied mean and symmetric positive-definite covariance.
      #
      # +random+ must be an explicit random source responding to #rand, such
      # as Random.new(seed). Bounds may be -Float::INFINITY or
      # Float::INFINITY. The returned error is a statistical Monte Carlo error
      # estimate (an approximate confidence-interval half-width), not a
      # deterministic integration error or an exact bound.
      #
      # Options:
      # * +epsilon+ - requested error estimate
      # * +alpha+ - multiplier applied to the estimated standard error
      # * +max_iterations+ - maximum number of Monte Carlo samples
      def cdf(lower, upper, mean:, covariance:, random:, epsilon: 0.0001,
              alpha: 2.5, max_iterations: 10_000)
        epsilon = Float(epsilon)
        alpha = Float(alpha)
        iterations = Integer(max_iterations)
        raise ArgumentError, 'epsilon must be positive and finite' unless epsilon.positive? && epsilon.finite?
        raise ArgumentError, 'alpha must be positive and finite' unless alpha.positive? && alpha.finite?
        raise ArgumentError, 'max_iterations must be positive' unless iterations.positive?
        raise ArgumentError, 'random must respond to #rand' unless random.respond_to?(:rand)

        lower = to_vector(lower, :lower)
        upper = to_vector(upper, :upper)
        mean = to_vector(mean, :mean)
        matrix = to_matrix(covariance, mean.size)
        bounds = validate_bounds!(lower, upper, mean)
        lower = bounds.fetch(:lower)
        upper = bounds.fetch(:upper)
        cholesky = positive_definite_cholesky(matrix)

        estimate_rectangle(lower, upper, cholesky, random, epsilon, alpha, iterations)
      end

      private

      def estimate_rectangle(lower, upper, cholesky, random, epsilon, alpha, iterations)
        dimension = mean_size(cholesky)
        return exact_rectangle(0.0) if lower.each_with_index.any? { |bound, i| bound == upper[i] }
        return exact_rectangle(1.0) if full_space?(lower, upper)

        mean_estimate = 0.0
        mean_square_delta = 0.0
        standard_error = 0.0
        used_samples = 0
        minimum_samples = [30, iterations].min

        iterations.times do |index|
          weight = rectangle_sample_weight(lower, upper, cholesky, random, dimension)
          sample_count = index + 1
          used_samples = sample_count
          delta = weight - mean_estimate
          mean_estimate += delta / sample_count
          mean_square_delta += delta * (weight - mean_estimate)

          if sample_count >= minimum_samples
            variance = mean_square_delta / (sample_count - 1)
            standard_error = Math.sqrt([0.0, variance].max / sample_count)
            break if alpha * standard_error <= epsilon
          end
        end

        error_estimate = alpha * standard_error
        variance_estimate = used_samples > 1 ? mean_square_delta / (used_samples - 1) : 0.0
        deterministic = diagonal?(cholesky)
        error_kind = :monte_carlo_confidence_interval_half_width
        if !deterministic && variance_estimate.zero?
          standard_error = [standard_error, 0.5 / Math.sqrt(used_samples)].max
          error_estimate = alpha * standard_error
          error_kind = :monte_carlo_conservative_standard_error_bound
        end
        target_reached = error_estimate <= epsilon && (deterministic || variance_estimate.positive?)
        {
          estimate: mean_estimate,
          error_estimate: error_estimate,
          standard_error: standard_error,
          target_reached: target_reached,
          samples: used_samples,
          error_kind: error_kind,
          error_description: 'Statistical Monte Carlo error estimate; not an exact probability or deterministic integration bound.'
        }
      end

      def exact_rectangle(estimate)
        {
          estimate: estimate,
          error_estimate: 0.0,
          standard_error: 0.0,
          target_reached: true,
          samples: 0,
          error_kind: :deterministic_degenerate_rectangle,
          error_description: 'Exact result for a degenerate empty or complete rectangular region.'
        }
      end

      def full_space?(lower, upper)
        lower.each_with_index.all? { |bound, i| bound == -Float::INFINITY && upper[i] == Float::INFINITY }
      end

      def diagonal?(matrix)
        matrix.each_with_index.all? do |row, i|
          row.each_with_index.all? { |entry, j| i == j || entry == 0.0 }
        end
      end

      def rectangle_sample_weight(lower, upper, cholesky, random, dimension)
        normal_values = Array.new(dimension, 0.0)
        weight = 1.0

        dimension.times do |i|
          location = 0.0
          i.times { |j| location += cholesky[i][j] * normal_values[j] }
          scale = cholesky[i][i]
          lower_probability = bound_cdf(lower[i], location, scale)
          upper_probability = bound_cdf(upper[i], location, scale)
          conditional_probability = [0.0, upper_probability - lower_probability].max
          weight *= conditional_probability

          unless i == dimension - 1
            if weight.zero?
              normal_values[i] = 0.0
            else
              uniform = lower_probability + random.rand * (upper_probability - lower_probability)
              uniform = [[uniform, 1e-300].max, 1.0 - 1e-300].min
              normal_values[i] = normal_quantile(uniform)
            end
          end
        end

        weight
      end

      def bound_cdf(bound, location, scale)
        return 0.0 if bound == -Float::INFINITY
        return 1.0 if bound == Float::INFINITY

        normal_cdf((bound - location) / scale)
      end

      def normal_cdf(value)
        0.5 * Math.erfc(-value / Math.sqrt(2.0))
      end

      def normal_quantile(probability)
        return -Float::INFINITY if probability.zero?
        return Float::INFINITY if probability == 1.0

        lower_breakpoint = 0.02425
        upper_breakpoint = 1.0 - lower_breakpoint

        if probability < lower_breakpoint
          q = Math.sqrt(-2.0 * Math.log(probability))
          (((((INVERSE_NORMAL_C[0] * q + INVERSE_NORMAL_C[1]) * q + INVERSE_NORMAL_C[2]) * q + INVERSE_NORMAL_C[3]) * q + INVERSE_NORMAL_C[4]) * q + INVERSE_NORMAL_C[5]) /
            ((((INVERSE_NORMAL_D[0] * q + INVERSE_NORMAL_D[1]) * q + INVERSE_NORMAL_D[2]) * q + INVERSE_NORMAL_D[3]) * q + 1.0)
        elsif probability <= upper_breakpoint
          q = probability - 0.5
          r = q * q
          (((((INVERSE_NORMAL_A[0] * r + INVERSE_NORMAL_A[1]) * r + INVERSE_NORMAL_A[2]) * r + INVERSE_NORMAL_A[3]) * r + INVERSE_NORMAL_A[4]) * r + INVERSE_NORMAL_A[5]) * q /
            (((((INVERSE_NORMAL_B[0] * r + INVERSE_NORMAL_B[1]) * r + INVERSE_NORMAL_B[2]) * r + INVERSE_NORMAL_B[3]) * r + INVERSE_NORMAL_B[4]) * r + 1.0)
        else
          q = Math.sqrt(-2.0 * Math.log(1.0 - probability))
          -(((((INVERSE_NORMAL_C[0] * q + INVERSE_NORMAL_C[1]) * q + INVERSE_NORMAL_C[2]) * q + INVERSE_NORMAL_C[3]) * q + INVERSE_NORMAL_C[4]) * q + INVERSE_NORMAL_C[5]) /
            ((((INVERSE_NORMAL_D[0] * q + INVERSE_NORMAL_D[1]) * q + INVERSE_NORMAL_D[2]) * q + INVERSE_NORMAL_D[3]) * q + 1.0)
        end
      end

      def mean_size(matrix)
        matrix.size
      end

      def to_vector(value, name)
        vector = if value.respond_to?(:to_ary)
                   value.to_ary
                 elsif value.respond_to?(:to_a)
                   value.to_a
                 end
        raise ArgumentError, "#{name} must be an array" unless vector.is_a?(Array)

        vector
      end

      def to_matrix(value, dimension)
        rows = value.respond_to?(:to_a) ? value.to_a : nil
        raise ArgumentError, 'covariance must be a matrix' unless rows.is_a?(Array) && rows.size == dimension
        raise ArgumentError, 'covariance must not be empty' if dimension.zero?

        converted_rows = rows.map do |row|
          converted_row = if row.respond_to?(:to_ary)
                            row.to_ary
                          elsif row.respond_to?(:to_a)
                            row.to_a
                          end
          unless converted_row.is_a?(Array) && converted_row.size == dimension
            raise ArgumentError, 'covariance must be a square matrix'
          end

          converted_row.map do |entry|
            converted_entry = Float(entry)
            unless converted_entry.finite?
              raise ArgumentError, 'covariance entries must be finite numeric values'
            end

            converted_entry
          end
        end

        converted_rows
      rescue TypeError
        raise ArgumentError, 'covariance entries must be finite numeric values'
      end

      def validate_bounds!(lower, upper, mean)
        dimension = mean.size
        raise ArgumentError, 'mean must not be empty' if dimension.zero?
        raise ArgumentError, 'lower and upper must match mean dimension' unless lower.size == dimension && upper.size == dimension

        mean = mean.map do |value|
          converted = Float(value)
          raise ArgumentError, 'mean entries must be finite numeric values' unless converted.finite?

          converted
        end

        normalized_lower = lower.each_with_index.map { |value, index| normalize_bound(value, :lower, index, mean[index]) }
        normalized_upper = upper.each_with_index.map { |value, index| normalize_bound(value, :upper, index, mean[index]) }

        dimension.times do |i|
          next if normalized_lower[i] <= normalized_upper[i]

          raise ArgumentError, "lower bound at dimension #{i} must be <= upper bound"
        end

        { lower: normalized_lower, upper: normalized_upper }
      end

      def normalize_bound(value, name, index, mean)
        value = -Float::INFINITY if value.nil? && name == :lower
        value = Float::INFINITY if value.nil? && name == :upper
        converted = Float(value)
        converted - mean
      rescue TypeError, ArgumentError
        raise ArgumentError, "#{name} bound at dimension #{index} must be numeric or infinite"
      end

      def positive_definite_cholesky(matrix)
        dimension = matrix.size
        scale = dimension * matrix.map { |row| row.map { |entry| entry.abs }.max }.max
        tolerance = 1e-14 * [1.0, scale].max

        dimension.times do |i|
          dimension.times do |j|
            next if (matrix[i][j] - matrix[j][i]).abs <= tolerance

            raise ArgumentError, 'covariance must be symmetric'
          end
        end

        cholesky = Array.new(dimension) { Array.new(dimension, 0.0) }

        dimension.times do |i|
          diagonal = matrix[i][i]
          i.times { |j| diagonal -= cholesky[i][j]**2 }
          raise ArgumentError, 'covariance must be positive definite' unless diagonal.positive? && diagonal.finite?

          cholesky[i][i] = Math.sqrt(diagonal)

          (i + 1...dimension).each do |row|
            dot_product = 0.0
            i.times { |j| dot_product += cholesky[row][j] * cholesky[i][j] }
            cholesky[row][i] = (matrix[row][i] - dot_product) / cholesky[i][i]
          end
        end

        cholesky
      end
    end
  end
end
