module Distribution
  module Stable
    module Ruby_
      extend self

      TWO_OVER_PI = 2.0 / Math::PI
      DEFAULT_TOLERANCE = 1e-10
      MAX_DEPTH = 24
      KRONROD_NODES = [-0.9914553711208126, -0.9491079123427585,
                       -0.8648644233597691, -0.7415311855993944,
                       -0.5860872354676911, -0.4058451513773972,
                       -0.2077849550078985, 0.0, 0.2077849550078985,
                       0.4058451513773972, 0.5860872354676911,
                       0.7415311855993944, 0.8648644233597691,
                       0.9491079123427585, 0.9914553711208126].freeze
      KRONROD_WEIGHTS = [0.022935322010529225, 0.06309209262997855,
                         0.10479001032225018, 0.14065325971552592,
                         0.1690047266392679, 0.1903505780647854,
                         0.2044329400752989, 0.20948214108472783,
                         0.2044329400752989, 0.1903505780647854,
                         0.1690047266392679, 0.14065325971552592,
                         0.10479001032225018, 0.06309209262997855,
                         0.022935322010529225].freeze
      GAUSS_WEIGHTS = [0.0, 0.1294849661688697, 0.0, 0.27970539148927664,
                       0.0, 0.38183005050511896, 0.0, 0.4179591836734694,
                       0.0, 0.38183005050511896, 0.0, 0.27970539148927664,
                       0.0, 0.1294849661688697, 0.0].freeze

      def parameters_from_s1(alpha:, beta:, delta:, gamma:)
        validate_parameters(alpha, beta, 1.0, gamma)
        location =
          if alpha == 1.0
            delta + TWO_OVER_PI * beta * gamma * Math.log(gamma)
          else
            delta + beta * gamma * stable_tan(alpha)
          end
        { alpha: alpha.to_f, beta: beta.to_f, delta: location.to_f,
          gamma: gamma.to_f }
      end

      def parameters_to_s1(alpha:, beta:, delta:, gamma:)
        validate_parameters(alpha, beta, 1.0, gamma)
        location =
          if alpha == 1.0
            delta - TWO_OVER_PI * beta * gamma * Math.log(gamma)
          else
            delta - beta * gamma * stable_tan(alpha)
          end
        { alpha: alpha.to_f, beta: beta.to_f, delta: location.to_f,
          gamma: gamma.to_f }
      end

      def pdf(x, alpha:, beta: 0, delta: 0, gamma: 1,
              diagnostics: false, tolerance: DEFAULT_TOLERANCE)
        result = evaluate(:pdf, x, alpha, beta, delta, gamma, tolerance)
        diagnostics ? result : fetch_value(result)
      end

      def cdf(x, alpha:, beta: 0, delta: 0, gamma: 1,
              diagnostics: false, tolerance: DEFAULT_TOLERANCE)
        result = evaluate(:cdf, x, alpha, beta, delta, gamma, tolerance)
        diagnostics ? result : fetch_value(result)
      end

      def rng(alpha:, beta: 0, delta: 0, gamma: 1, source: nil)
        validate_parameters(alpha, beta, 0.0, gamma)
        uniform = uniform_source(source)
        spare = []

        lambda do
          value =
            if alpha == 2.0
              normal_sample(uniform, spare) * Math.sqrt(2.0) * gamma
            elsif alpha == 1.0 && beta == 0.0
              uniform_value = 0.0
              loop do
                uniform_value = uniform.call
                break if uniform_value > 0.0 && uniform_value < 1.0
              end
              Math.tan(Math::PI * (uniform_value - 0.5)) * gamma
            else
              standard_s1_sample(alpha, beta, uniform) * gamma
            end
          value + delta
        end
      end

      private

      def validate_parameters(alpha, beta, delta, gamma)
        unless alpha.is_a?(Numeric) && alpha > 0.0 && alpha <= 2.0
          fail ArgumentError, 'alpha must be in (0, 2]'
        end
        unless beta.is_a?(Numeric) && beta >= -1.0 && beta <= 1.0
          fail ArgumentError, 'beta must be in [-1, 1]'
        end
        unless gamma.is_a?(Numeric) && gamma > 0.0
          fail ArgumentError, 'gamma must be positive'
        end
        unless delta.is_a?(Numeric) && delta.finite?
          fail ArgumentError, 'delta must be finite'
        end
      end

      def stable_tan(alpha)
        return 0.0 if alpha == 2.0
        Math.tan(Math::PI * alpha / 2.0)
      end

      def fetch_value(result)
        fail RuntimeError, result[:message] unless result[:converged]
        result[:value]
      end

      def evaluate(kind, x, alpha, beta, delta, gamma, tolerance)
        validate_parameters(alpha, beta, delta, gamma)
        unless tolerance.is_a?(Numeric) && tolerance > 0.0 && tolerance.finite?
          fail ArgumentError, 'tolerance must be positive and finite'
        end

        z = (x.to_f - delta.to_f) / gamma.to_f
        value, error, method =
          if alpha == 2.0
            normal_result(kind, z)
          elsif alpha == 1.0 && beta == 0.0
            cauchy_result(kind, z)
          elsif alpha < 1.0 && beta == 1.0 && z <= -stable_tan(alpha)
            boundary_result(kind, gamma)
          elsif alpha < 1.0 && beta == -1.0 && z >= stable_tan(alpha)
            boundary_result(kind, gamma, true)
          else
            numerical_result(kind, z, alpha, beta, gamma, tolerance)
          end

        { value: value, error: error, requested_tolerance: tolerance,
          converged: error <= tolerance, method: method, x: x.to_f,
          alpha: alpha.to_f, beta: beta.to_f, delta: delta.to_f,
          gamma: gamma.to_f,
          message: "stable #{kind} integration did not reach requested "\
                   "tolerance #{tolerance}: estimated error #{error}" }
      end

      def normal_result(kind, z)
        if kind == :pdf
          [normal_pdf(z / Math.sqrt(2.0)) / Math.sqrt(2.0), 0.0, :normal]
        else
          [normal_cdf(z / Math.sqrt(2.0)), 0.0, :normal]
        end
      end

      def cauchy_result(kind, z)
        if kind == :pdf
          [1.0 / (Math::PI * (1.0 + z * z)), 0.0, :cauchy]
        else
          [0.5 + Math.atan(z) / Math::PI, 0.0, :cauchy]
        end
      end

      def boundary_result(kind, _gamma, upper = false)
        if kind == :pdf
          [0.0, 0.0, :support_boundary]
        else
          [(upper ? 1.0 : 0.0), 0.0, :support_boundary]
        end
      end

      def normal_pdf(z)
        Math.exp(-0.5 * z * z) / Distribution::SQ2PI
      end

      def normal_cdf(z)
        return 0.0 if z < -12.0
        return 1.0 if z > 12.0
        0.5 * (1.0 + Math.erf(z / Math.sqrt(2.0)))
      end

      def numerical_result(kind, z, alpha, beta, gamma, tolerance)
        integrand = integral_function(kind, z, alpha.to_f, beta.to_f)
        head_power = alpha < 1.0 ? 2.0 / alpha : 2.0
        pieces = []
        if alpha > 1.0
          head_segments = [64, (z.abs * 8.0).to_i + 1].max
          head_segments.times do |index|
            power = 2.0 / alpha
            left = index / head_segments.to_f
            right = (index + 1) / head_segments.to_f
            pieces << [left, right, power, :head, 1.0]
          end
          tail_points = [0.0] + 8.downto(1).map { |index| 1.0 / (2**index) } + [1.0]
          tail_points.each_cons(2) do |left, right|
            pieces << [left, right, 1.0, :tail]
          end
        else
          pieces << [0.0, 1.0, head_power, :head]
          tail_points = [0.0] + 8.downto(1).map { |index| 1.0 / (2**index) } + [1.0]
          tail_points.each_cons(2) do |left, right|
            pieces << [left, right, 1.0, :tail]
          end
        end
        scale = kind == :pdf ? 1.0 / (Math::PI * gamma) : 1.0 / Math::PI
        total = 0.0
        error = 0.0
        converged = true
        evaluations = 0
        method = kind == :pdf ? :characteristic_function_pdf\
                              : :characteristic_function_cdf

        pieces.each do |left, right, power, segment, segment_scale|
          part = adaptive_integral(integrand, left, right, power, segment,
                                   tolerance / pieces.count.to_f, evaluations,
                                   segment_scale || 1.0)
          evaluations = part[:evaluations]
          total += part[:value]
          error += part[:error]
          converged &&= part[:converged]
        end

        value = total.abs * scale
        if kind == :cdf
          value = 0.5 + (total >= 0.0 ? value : -value)
          if (value < 0.0 && value > -tolerance) ||
             (value > 1.0 && value < 1.0 + tolerance)
            value = [[value, 0.0].max, 1.0].min
          end
        end
        [value, converged ? error.abs * scale : Float::INFINITY, method]
      end

      def integral_function(kind, z, alpha, beta)
        lambda do |t|
          return 0.0 if t == 0.0
          phase = stable_phase(alpha, beta, t)
          amplitude = Math.exp(-t**alpha)
          return 0.0 if amplitude == 0.0
          angle = z * t - phase
          if kind == :pdf
            amplitude * Math.cos(angle)
          else
            amplitude * Math.sin(angle) / t
          end
        end
      end

      def stable_phase(alpha, beta, t)
        return -TWO_OVER_PI * beta * t * Math.log(t) if alpha == 1.0

        epsilon = alpha - 1.0
        y = epsilon * Math.log(t)
        difference =
          if y.abs < 0.05
            t * expm1(y)
          else
            t * (t**epsilon - 1.0)
          end
        beta * stable_tan(alpha) * difference
      end

      def expm1(value)
        if value.abs < 1e-4
          value * (1.0 + value * (0.5 + value / 6.0))
        else
          Math.exp(value) - 1.0
        end
      end

      def adaptive_integral(function, left, right, power, segment,
                            tolerance, evaluations, scale = 1.0)
        state = { evaluations: evaluations }
        result = adaptive_kronrod(function, left, right, power, segment,
                                  scale, tolerance, 0, state)
        result[:evaluations] = state[:evaluations]
        result
      end

      def adaptive_kronrod(function, left, right, power, segment, scale,
                           tolerance, depth, state)
        kronrod = 0.0
        gauss = 0.0
        width = right - left

        KRONROD_NODES.each_with_index do |node, index|
          point = (left + right + node * width) / 2.0
          value = transformed_value(function, point, power, segment, scale)
          state[:evaluations] += 1
          kronrod += KRONROD_WEIGHTS[index] * value * width / 2.0
          gauss += GAUSS_WEIGHTS[index] * value * width / 2.0
        end

        error = (kronrod - gauss).abs
        return { value: kronrod, error: error, converged: true } if error <= tolerance
        if depth >= MAX_DEPTH || state[:evaluations] > 200_000
          return { value: kronrod, error: Float::INFINITY, converged: false }
        end

        middle = (left + right) / 2.0
        left_result = adaptive_kronrod(function, left, middle, power, segment,
                                       scale, tolerance, depth + 1, state)
        right_result = adaptive_kronrod(function, middle, right, power, segment,
                                        scale, tolerance, depth + 1, state)
        { value: left_result[:value] + right_result[:value],
          error: left_result[:error] + right_result[:error],
          converged: left_result[:converged] && right_result[:converged] }
      end

      def finite_number?(value)
        value.is_a?(Numeric) && value.finite?
      end

      def transformed_value(function, point, power, segment, scale)
        if segment == :tail
          t = point**(-power)
          jacobian = power * point**(-power - 1.0)
          return 0.0 if t.finite? && Math.exp(-t**power) == 0.0
          return Float::NAN unless t.finite? && jacobian.finite?
        else
          return 0.0 if point.zero? && power < 1.0
          t = scale * point**power
          jacobian = scale * power * point**(power - 1.0)
        end
        value = function.call(t) * jacobian
        finite_number?(value) ? value : 0.0
      end

      def uniform_source(source)
        if source.nil?
          random = Random.new
          -> { random.rand }
        elsif source.respond_to?(:rand)
          -> { source.rand }
        elsif source.respond_to?(:call)
          source
        elsif source.respond_to?(:next)
          -> { source.next }
        else
          fail ArgumentError, 'source must respond to rand or call'
        end
      end

      def normal_sample(uniform, spare)
        return spare.pop unless spare.empty?
        radius_squared = 0.0
        x = y = 0.0
        loop do
          x = 2.0 * uniform.call - 1.0
          y = 2.0 * uniform.call - 1.0
          radius_squared = x * x + y * y
          break if radius_squared > 0.0 && radius_squared < 1.0
        end
        multiplier = Math.sqrt(-2.0 * Math.log(radius_squared) / radius_squared)
        spare.push(y * multiplier)
        x * multiplier
      end

      def standard_s1_sample(alpha, beta, uniform)
        theta = 0.0
        w = 0.0
        loop do
          theta = Math::PI * (uniform.call - 0.5)
          w = -Math.log(uniform.call)
          break if theta.abs < Math::PI / 2.0 && w > 0.0 && w.finite?
        end

        if alpha == 1.0
          TWO_OVER_PI * ((Math::PI / 2.0 + beta * theta) * Math.tan(theta) -
            beta * Math.log((Math::PI / 2.0 * w * Math.cos(theta)) /
                           (Math::PI / 2.0 + beta * theta)))
        else
          tangent = stable_tan(alpha)
          zeta = beta * tangent
          theta_zero = Math.atan(zeta) / alpha
          adjusted_theta = alpha * theta
          denominator = Math.cos(theta) / Math.tan(alpha * (theta_zero + theta)) +
                        Math.sin(theta)
          numerator = Math.cos(adjusted_theta) +
                      Math.sin(adjusted_theta) * Math.tan(theta) -
                      zeta * (Math.sin(adjusted_theta) -
                              Math.cos(adjusted_theta) * Math.tan(theta))
          base = numerator / w
          base = Float::EPSILON if base <= 0.0 && base > -1e-12
          (w / denominator) * base**(1.0 / alpha) - zeta
        end
      end
    end
  end
end
