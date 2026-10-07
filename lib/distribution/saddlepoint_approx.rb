module Distribution
  # Saddlepoint tail approximations built only from a cumulant generating
  # function and its derivatives.
  class SaddlepointApproximation
    class DomainError < StandardError; end

    Result = Struct.new(:probability, :approximate, :distribution, :status,
                        :saddlepoint, :requested_x, :lattice_point, :diagnostics) do
      alias_method :value, :probability
      alias_method :approximate?, :approximate

      def to_h
        members.each_with_object({}) { |member, hash| hash[member] = self[member] }
      end
    end

    TWO_PI = 2.0 * Math::PI
    SQRT_TWO_PI = Math.sqrt(TWO_PI)
    SERIES_LIMIT = 1.0e-3
    MAX_BRACKET_STEPS = 80
    MAX_BISECTION_STEPS = 80

    def initialize(cgf:, cgf_derivative:, second_derivative:,
                   third_derivative: nil, kind: :continuous,
                   lattice_spacing: 1.0, lattice_origin: 0.0)
      unless cgf.respond_to?(:call) && cgf_derivative.respond_to?(:call) &&
             second_derivative.respond_to?(:call)
        raise ArgumentError, 'cgf, cgf_derivative and second_derivative must respond to #call'
      end

      if kind == :continuous
        @continuous = true
        @lattice_spacing = nil
      else
        unless kind == :lattice
          raise ArgumentError, "kind must be :continuous or :lattice"
        end
        @continuous = false
        @lattice_spacing = Float(lattice_spacing)
        unless @lattice_spacing.finite? && @lattice_spacing > 0
          raise ArgumentError, 'lattice_spacing must be a finite positive number'
        end
      end

      @cgf = cgf
      @cgf_derivative = cgf_derivative
      @second_derivative = second_derivative
      @third_derivative = third_derivative
      @lattice_origin = Float(lattice_origin)
      unless @lattice_origin.finite?
        raise ArgumentError, 'lattice_origin must be finite'
      end

      @cgf_at_zero = call_callback(:cgf, 0.0)
      @mean = call_callback(:cgf_derivative, 0.0)
      @variance = call_callback(:second_derivative, 0.0)

      unless finite?(@cgf_at_zero) && finite?(@mean) && finite?(@variance) &&
             @variance > 0
        raise ArgumentError, 'K(0), K\'(0) and a positive K\'\'(0) are required'
      end
    end

    def continuous?
      @continuous
    end

    attr_reader :mean, :variance

    def cdf(x)
      x = Float(x)
      if continuous?
        if x >= @mean
          with_request(complement(direct_upper_tail(x)), x)
        else
          with_request(reflected_upper_tail(-x), x)
        end
      else
        point = lattice_floor(x)
        if point <= @mean
          with_request(reflected_upper_tail(-point, point == @mean, point), x, point)
        else
          upper = direct_upper_tail(point + @lattice_spacing, normalized_cgf,
                                    @cgf_derivative, @second_derivative,
                                    @third_derivative, false, false, point)
          with_request(complement(upper), x, point)
        end
      end
    end

    def upper_tail(x)
      x = Float(x)
      if continuous?
        if x >= @mean
          with_request(direct_upper_tail(x), x)
        else
          lower = reflected_upper_tail(-x)
          with_request(complement(lower), x)
        end
      else
        point = lattice_ceil(x)
        if point >= @mean
          upper = direct_upper_tail(point, normalized_cgf, @cgf_derivative,
                                    @second_derivative, @third_derivative, false,
                                    point == @mean, point)
          with_request(upper, x, point)
        else
          threshold = point - @lattice_spacing
          lower = reflected_upper_tail(-threshold, threshold == @mean, threshold)
          with_request(complement(lower), x, point)
        end
      end
    end

    alias_method :survival, :upper_tail

    private

    def with_request(result, x, lattice_point = nil)
      annotated = result.dup
      annotated.requested_x = x
      annotated.lattice_point = lattice_point unless continuous?
      annotated
    end

    def reflected_upper_tail(target, mean_limit = false, reported_lattice_point = target)
      reflection = lambda { |t| @cgf.call(-t) }
      derivative = lambda { |t| -@cgf_derivative.call(-t) }
      second = lambda { |t| @second_derivative.call(-t) }
      third = @third_derivative.nil? ? nil : lambda { |t| -@third_derivative.call(-t) }

      direct_upper_tail(target, reflection, derivative, second, third, true,
                        mean_limit, reported_lattice_point)
    end

    def direct_upper_tail(x, cgf = normalized_cgf, derivative = @cgf_derivative,
                          second = @second_derivative, third = @third_derivative,
                          reflected = false, mean_limit = x == @mean,
                          reported_lattice_point = x)
      point = continuous? ? nil : reported_lattice_point
      if mean_limit
        t = 0.0
        variance = @variance
        exponent = 0.0
        root = {
          iterations: 0,
          bracket_iterations: 0,
          root_iterations: 0,
          residual: 0.0,
          variance: variance,
          bracket: [0.0, 0.0],
          reflected: reflected
        }
      else
        root = solve_saddlepoint(x, derivative, second, reflected)
        return failure_result(x, point, root) unless root[:status] == :ok

        t = root[:saddlepoint]
        variance = root[:variance]
        k_value = evaluate(cgf, t)
        unless k_value[:status] == :ok && finite?(variance) && variance > 0
          return failure_result(x, point, root.merge(status: :domain_error,
                                                    message: 'The saddlepoint is outside the CGF domain'))
        end

        exponent = t * x - k_value[:value]
        if exponent < -Float::EPSILON
          return failure_result(x, point, root.merge(status: :numerical_error,
                                                    message: 'Negative saddlepoint exponent encountered'))
        end
      end

      probability, method_name =
        if continuous?
          continuous_probability(t, variance, exponent, third, second)
        else
          lattice_probability(t, variance, exponent, third, second)
        end

      raw_probability = probability
      clipped = probability < 0 || probability > 1
      probability = [[probability, 0.0].max, 1.0].min
      diagnostics = root.reject { |key| [:status, :message].include?(key) }.merge(
        approximation: method_name,
        exponent: exponent,
        raw_probability: raw_probability,
        clipped: clipped
      )

      Result.new(probability, true, distribution_description, :ok, t, x, point,
                 diagnostics)
    end

    def continuous_probability(t, variance, exponent, third, second)
      if t.abs <= SERIES_LIMIT
        k3 = third_cumulant_at_zero(third, second)
        correction = Math.sqrt(variance) * t.abs + k3 / (6.0 * variance**1.5)
        probability = 0.5 - normal_pdf(0.0) * correction
        return probability, :lugannani_rice_mean_expansion if t != 0
        return probability, :lugannani_rice_mean_limit
      end

      w = Math.sqrt([2.0 * exponent, 0.0].max)
      u = t * Math.sqrt(variance)
      probability = normal_survival(w) + normal_pdf(w) * (1.0 / u - 1.0 / w)
      [probability, :lugannani_rice]
    end

    def lattice_probability(t, variance, exponent, third, second)
      if t.abs <= SERIES_LIMIT
        k3 = third_cumulant_at_zero(third, second)
        local_mass = @lattice_spacing / (SQRT_TWO_PI * Math.sqrt(variance))
        correction = Math.sqrt(variance) * t.abs + k3 / (6.0 * variance**1.5)
        probability = 0.5 + local_mass / 2.0 - normal_pdf(0.0) * correction
        return probability, :daniels_mean_expansion if t != 0
        return probability, :daniels_mean_limit
      end

      w = Math.sqrt([2.0 * exponent, 0.0].max)
      lattice_factor = (1.0 - Math.exp(-@lattice_spacing * t)) / @lattice_spacing
      u = lattice_factor * Math.sqrt(variance)
      probability = normal_survival(w) + normal_pdf(w) * (1.0 / u - 1.0 / w)
      [probability, :daniels_lattice]
    end

    def third_cumulant_at_zero(third, second)
      value = if third.respond_to?(:call)
                third.call(0.0)
              else
                epsilon = 1.0e-5
                (second.call(epsilon) - second.call(-epsilon)) / (2.0 * epsilon)
              end
      raise DomainError, 'Third cumulant is not finite' unless finite?(value)
      value
    rescue StandardError
      0.0
    end

    def solve_saddlepoint(target, derivative, second, reflected)
      bracket = build_bracket(target, derivative)
      return bracket unless bracket[:status] == :ok

      low = bracket[:low]
      high = bracket[:high]
      iterations = 0
      residual = Float::INFINITY

      MAX_BISECTION_STEPS.times do
        iterations += 1
        middle = low + (high - low) / 2.0
        middle_value = evaluate(derivative, middle)
        if middle_value[:status] != :ok
          high = middle
          next
        end

        residual = middle_value[:value] - target
        if residual < 0
          low = middle
        else
          high = middle
        end

        break if residual.abs <= [1.0, target.abs].max * 1.0e-11 &&
                 high - low <= [1.0, target.abs].max * 1.0e-10
      end

      bracket_width = high - low
      if residual.abs > [1.0, target.abs].max * 1.0e-8 ||
         bracket_width > [1.0, target.abs].max * 1.0e-8
        return {
          status: :root_did_not_converge,
          requested_target: target,
          bracket: { low: low, high: high },
          iterations: bracket[:bracket_iterations] + iterations,
          message: 'Saddlepoint root search did not converge to requested tolerance'
        }
      end

      saddlepoint = low + (high - low) / 2.0
      first = evaluate(derivative, saddlepoint)
      variance_value = evaluate(second, saddlepoint)
      if first[:status] != :ok || variance_value[:status] != :ok
        return not_found(target, bracket, iterations, 'The saddlepoint is outside the CGF domain')
      end

      variance = variance_value[:value]
      unless finite?(variance) && variance > 0
        return not_found(target, bracket, iterations, 'K\'\'(t) must be positive at the saddlepoint')
      end

      residual = first[:value] - target
      if variance < 1.0e-10 * [1.0, target.abs].max**2
        return not_found(target, bracket, iterations, 'No finite interior saddlepoint exists')
      end

      {
        status: :ok,
        saddlepoint: saddlepoint,
        iterations: bracket[:bracket_iterations] + iterations,
        bracket_iterations: bracket[:bracket_iterations],
        root_iterations: iterations,
        residual: residual,
        variance: variance,
        bracket: [low, high],
        reflected: reflected
      }
    end

    def build_bracket(target, derivative)
      low = 0.0
      high = 1.0 / 8.0
      low_value = evaluate(derivative, low)
      return not_found(target, {}, 0, 'K\'(0) could not be evaluated') if low_value[:status] != :ok

      iterations = 0
      MAX_BRACKET_STEPS.times do
        iterations += 1
        high_value = evaluate(derivative, high)

        if high_value[:status] == :ok
          if high_value[:value] >= target
            return { status: :ok, low: low, high: high, bracket_iterations: iterations }
          end

          low = high
          high *= 2.0
        else
          barrier = high
          MAX_BISECTION_STEPS.times do
            iterations += 1
            middle = low + (barrier - low) / 2.0
            middle_value = evaluate(derivative, middle)

            if middle_value[:status] == :ok && middle_value[:value] >= target
              return { status: :ok, low: low, high: middle, bracket_iterations: iterations }
            end

            if middle_value[:status] == :ok
              low = middle
            else
              barrier = middle
            end
          end

          return not_found(target, { low: low, high: barrier }, iterations,
                           'No finite interior saddlepoint exists before the CGF boundary')
        end
      end

      not_found(target, { low: low, high: high }, iterations,
               'Unable to bracket a finite saddlepoint')
    end

    def evaluate(callback, point)
      value = callback.call(point)
      value = Float(value)
      finite?(value) ? { status: :ok, value: value } : { status: :domain_error }
    rescue StandardError
      { status: :domain_error }
    end

    def call_callback(name, point)
      value = instance_variable_get("@#{name}").call(point)
      Float(value)
    rescue StandardError => error
      raise DomainError, "#{name} failed at t=#{point}: #{error.message}"
    end

    def normalized_cgf
      @normalized_cgf ||= lambda { |t| @cgf.call(t) - @cgf_at_zero }
    end

    def lattice_floor(x)
      index = ((x - @lattice_origin) / @lattice_spacing).floor
      @lattice_origin + index * @lattice_spacing
    end

    def lattice_ceil(x)
      distance = x - @lattice_origin
      index = (distance / @lattice_spacing).ceil
      @lattice_origin + index * @lattice_spacing
    end

    def complement(result)
      return result unless result.status == :ok

      complemented = result.dup
      complemented.probability = 1.0 - result.probability
      complemented.diagnostics = result.diagnostics.merge(complemented: true)
      complemented
    end

    def failure_result(x, point, diagnostics)
      status = diagnostics[:status] || :saddlepoint_not_found
      Result.new(nil, false, distribution_description, status, nil, x, point,
                 diagnostics)
    end

    def not_found(target, bracket, iterations, message)
      {
        status: :saddlepoint_not_found,
        requested_target: target,
        bracket: bracket,
        iterations: iterations,
        message: message
      }
    end

    def normal_pdf(x)
      Math.exp(-x * x / 2.0) / SQRT_TWO_PI
    end

    def normal_survival(x)
      0.5 * Math.erfc(x / Math.sqrt(2.0))
    end

    def finite?(value)
      value.is_a?(Numeric) && value.finite?
    end

    def distribution_description
      continuous? ? :continuous : { lattice_spacing: @lattice_spacing,
                                    lattice_origin: @lattice_origin }
    end
  end
end
