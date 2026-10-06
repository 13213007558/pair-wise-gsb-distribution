module Distribution
  module Stable
    module Ruby_
      GK15_X = [0.9914553711208126392069,
                0.9491079123427585245262,
                0.8648644233597690727897,
                0.7415311855993944398639,
                0.5860872354676911302941,
                0.4058451513773971669066,
                0.2077849550078984676007].freeze
      GK15_W = [0.0229353220105292249637,
                0.0630920926299785532907,
                0.1047900103222501838399,
                0.1406532597155259187452,
                0.1690047266392679028266,
                0.1903505780647854099133,
                0.2044329400752988924142,
                0.2094821410847278280130].freeze
      GK7_W = [0.0,
               0.1294849661688696932706,
               0.0,
               0.2797053914892766679015,
               0.0,
               0.3818300505051189449504,
               0.0,
               0.4179591836734693877551].freeze

      class << self
        def pdf(x, alpha, beta = 0.0, scale = 1.0, location = 0.0, options = {})
          options = {} unless options.is_a?(Hash)
          alpha, beta, scale, location = validate_parameters(alpha, beta, scale, location)
          z = (x.to_f - location) / scale
          result = standard_pdf(z, alpha, beta, options)
          options[:diagnostics] ? add_scaled_diagnostics(result, scale) : result[:value] / scale
        end

        def cdf(x, alpha, beta = 0.0, scale = 1.0, location = 0.0, options = {})
          options = {} unless options.is_a?(Hash)
          alpha, beta, scale, location = validate_parameters(alpha, beta, scale, location)
          z = (x.to_f - location) / scale
          result = standard_cdf(z, alpha, beta, options)
          options[:diagnostics] ? result : result[:value]
        end

        def rng(alpha, beta = 0.0, scale = 1.0, location = 0.0, random_source = nil)
          alpha, beta, scale, location = validate_parameters(alpha, beta, scale, location)
          source = random_source || Random.new(Random.new_seed)
          if alpha == 2.0
            normal_state = { :available => false, :value => 0.0 }
            return lambda do
              if normal_state[:available]
                normal_state[:available] = false
                location + scale * normal_state[:value]
              else
                radius = Math.sqrt(-2.0 * Math.log(draw_uniform(source)))
                angle = 2.0 * Math::PI * draw_uniform(source)
                normal_state[:available] = true
                normal_state[:value] = radius * Math.sin(angle)
                location + scale * radius * Math.cos(angle)
              end
            end
          end
          lambda do
            uniform = draw_uniform(source)
            exponential = -Math.log(draw_uniform(source))
            location + scale * standard_sample(alpha, beta, uniform, exponential)
          end
        end

        def s1_parameters(alpha, beta, scale, location)
          validate_parameters(alpha, beta, scale, location)
          [alpha.to_f, beta.to_f, scale.to_f, location.to_f - location_shift(alpha, beta, scale)]
        end

        def from_s1_parameters(alpha, beta, scale, location_s1)
          validate_parameters(alpha, beta, scale, location_s1)
          [alpha.to_f, beta.to_f, scale.to_f, location_s1.to_f + location_shift(alpha, beta, scale)]
        end

        private

        def validate_parameters(alpha, beta, scale, location)
          alpha = alpha.to_f
          beta = beta.to_f
          scale = scale.to_f
          location = location.to_f
          unless alpha > 0.0 && alpha <= 2.0
            raise ArgumentError, 'alpha must satisfy 0 < alpha <= 2'
          end
          raise ArgumentError, 'beta must satisfy -1 <= beta <= 1' unless beta.between?(-1.0, 1.0)
          raise ArgumentError, 'scale must be positive' unless scale > 0.0
          [alpha, beta, scale, location]
        end

        def draw_uniform(source)
          if source.respond_to?(:rand)
            source.rand
          elsif source.respond_to?(:call)
            source.call
          else
            source.next
          end
        end

        def location_shift(alpha, beta, scale)
          if alpha == 1.0
            2.0 * beta * scale * Math.log(scale) / Math::PI
          else
            beta * scale * tan_pi_alpha_half(alpha)
          end
        end

        def tan_pi_alpha_half(alpha)
          epsilon = alpha - 1.0
          if epsilon.abs < 1.0e-3
            angle = Math::PI * epsilon / 2.0
            -1.0 / Math.tan(angle)
          else
            Math.tan(Math::PI * alpha / 2.0)
          end
        end

        def standard_sample(alpha, beta, theta, exponential)
          theta = theta * Math::PI - Math::PI / 2.0
          if alpha == 1.0
            beta_theta = beta * theta
            (2.0 / Math::PI) * (((Math::PI / 2.0) + beta_theta) * Math.tan(theta) -
              beta * Math.log((Math::PI * exponential * Math.cos(theta) / 2.0) /
              ((Math::PI / 2.0) + beta_theta)))
          elsif beta.zero? || alpha == 2.0
            angle = alpha * theta
            (Math.sin(angle) / Math.cos(theta)**(1.0 / alpha)) *
              (Math.cos((1.0 - alpha) * theta) / exponential)**((1.0 - alpha) / alpha)
          else
            q = beta * tan_pi_alpha_half(alpha)
            theta_zero = Math.atan(q) / alpha
            alpha_theta = alpha * theta
            tan_theta = Math.tan(theta)
            val = exponential / (Math.cos(theta) / Math.tan(alpha * (theta_zero + theta)) + Math.sin(theta))
            inner = (Math.cos(alpha_theta) + Math.sin(alpha_theta) * tan_theta -
              q * (Math.sin(alpha_theta) - Math.cos(alpha_theta) * tan_theta)) / exponential
            val * inner**(1.0 / alpha) - q
          end
        end

        def standard_pdf(x, alpha, beta, options)
          return success(exact_pdf(x, alpha, beta)) if exact_case?(alpha, beta)
          fourier_pdf(x, alpha, beta, options)
        end

        def standard_cdf(x, alpha, beta, options)
          return success(exact_cdf(x, alpha, beta)) if exact_case?(alpha, beta)
          fourier_cdf(x, alpha, beta, options)
        end

        def exact_case?(alpha, beta)
          alpha == 2.0 || (alpha == 1.0 && beta.zero?) || (alpha == 0.5 && beta == 1.0)
        end

        def exact_pdf(x, alpha, beta)
          if alpha == 2.0
            Normal::Ruby_.pdf(x / Math.sqrt(2.0)) / Math.sqrt(2.0)
          elsif alpha == 1.0
            1.0 / (Math::PI * (1.0 + x * x))
          else
            shifted = x + 1.0
            return 0.0 if shifted <= 0.0
            Math.exp(-1.0 / (2.0 * shifted)) /
              (Math.sqrt(2.0 * Math::PI) * shifted**1.5)
          end
        end

        def exact_cdf(x, alpha, beta)
          if alpha == 2.0
            Normal::Ruby_.cdf(x / Math.sqrt(2.0))
          elsif alpha == 1.0
            0.5 + Math.atan(x) / Math::PI
          else
            return 0.0 if x <= -1.0
            2.0 * (1.0 - Normal::Ruby_.cdf(1.0 / Math.sqrt(x + 1.0)))
          end
        end

        def fourier_pdf(x, alpha, beta, options)
          cutoff = options.fetch(:cutoff, fourier_cutoff(alpha, options))
          local_options = integration_options(options, alpha)
          breakpoints = fourier_breakpoints(cutoff)
          integral = integrate_adaptive(breakpoints.each_cons(2).to_a, local_options) do |frequency|
            phase = fourier_phase(frequency, alpha, beta)
            Math.exp(-frequency**alpha) * Math.cos(phase + frequency * x)
          end
          value = [0.0, integral[:value] / Math::PI].max
          success(value, integral, integral[:converged])
        end

        def fourier_cdf(x, alpha, beta, options)
          cutoff = options.fetch(:cutoff, fourier_cutoff(alpha, options))
          local_options = integration_options(options, alpha)
          breakpoints = fourier_breakpoints(cutoff)
          integral = integrate_adaptive(breakpoints.each_cons(2).to_a, local_options) do |frequency|
            phase = fourier_phase(frequency, alpha, beta)
            if frequency.zero?
              x + cdf_frequency_slope(alpha, beta)
            else
              Math.exp(-frequency**alpha) * Math.sin(phase + frequency * x) / frequency
            end
          end
          value = 0.5 + integral[:value] / Math::PI
          value = [[value, 0.0].max, 1.0].min
          success(value, integral, integral[:converged])
        end

        def fourier_phase(frequency, alpha, beta)
          if alpha == 1.0
            2.0 * beta * frequency * Math.log(frequency) / Math::PI
          else
            beta * tan_pi_alpha_half(alpha) *
              (frequency - frequency**alpha)
          end
        end

        def cdf_frequency_slope(alpha, beta)
          if alpha == 1.0
            2.0 * beta / Math::PI
          else
            beta * tan_pi_alpha_half(alpha)
          end
        end

        def fourier_cutoff(alpha, options)
          tolerance = options.fetch(:epsilon, 1.0e-10)
          target = -Math.log([tolerance * 0.01, 1.0e-300].max)
          [target**(1.0 / alpha), 10000.0].min
        end

        def fourier_breakpoints(cutoff)
          points = [0.0, 1.0e-9, 1.0e-7, 1.0e-5, 1.0e-3, 1.0e-2, 0.1, 1.0 / Math::E, 1.0]
          boundary = 2.0
          while boundary < cutoff
            points << boundary
            boundary *= 2.0
          end
          points << cutoff
          points
        end

        def integration_options(options, alpha)
          options.merge(
            :max_intervals => options.fetch(:max_intervals, 20000),
            :gk7_unreliable => alpha == 1.0
          )
        end

        def success(value, parent = nil, converged = true)
          result = { :value => value, :converged => converged }
          if parent
            result[:error] = parent[:error]
            result[:evaluations] = parent[:evaluations]
            result[:subdivisions] = parent[:subdivisions]
            result[:segments] = parent[:segments]
          end
          result
        end

        def add_scaled_diagnostics(result, scale)
          scaled = result.merge(:value => result[:value] / scale)
          scaled[:error] = scaled[:error] / scale if scaled[:error]
          scaled
        end

        def integrate_adaptive(segments, options)
          tolerance = options.fetch(:epsilon, 1.0e-10)
          max_intervals = options.fetch(:max_intervals, 1000)
          total = 0.0
          error = 0.0
          evaluations = 0
          subdivisions = 0
          pending = []

          segments.each do |left, right|
            estimate, segment_error, segment_evaluations = gk15_integrate(left, right) { |theta| yield(theta) }
            total += estimate
            error += segment_error.abs unless options[:gk7_unreliable]
            evaluations += segment_evaluations
            pending << [(options[:gk7_unreliable] ? 0.0 : segment_error.abs), left, right, estimate, segment_error]
          end

          target = [tolerance, tolerance * total.abs].max
          while error > target && subdivisions < max_intervals && !pending.empty?
            pending.sort_by! { |segment| -segment.first }
            _, left, right, old_estimate, old_error = pending.shift
            middle = (left + right) / 2.0
            left_estimate, left_error, left_evaluations = gk15_integrate(left, middle) { |theta| yield(theta) }
            right_estimate, right_error, right_evaluations = gk15_integrate(middle, right) { |theta| yield(theta) }
            new_estimate = left_estimate + right_estimate
            new_error = (new_estimate - old_estimate).abs
            total += new_estimate - old_estimate
            error += new_error if options[:gk7_unreliable]
            error += new_error - old_error.abs unless options[:gk7_unreliable]
            evaluations += left_evaluations + right_evaluations
            subdivisions += 1
            target = [tolerance, tolerance * total.abs].max
            pending << [new_error / 2.0, left, middle, left_estimate, new_error / 2.0]
            pending << [new_error / 2.0, middle, right, right_estimate, new_error / 2.0]
          end

          converged = subdivisions < max_intervals &&
            error <= [tolerance, tolerance * total.abs].max
          diagnostics = {
            :value => total,
            :error => error.abs,
            :converged => converged,
            :evaluations => evaluations,
            :subdivisions => subdivisions,
            :segments => segments.length,
            :tolerance => tolerance,
            :max_intervals => max_intervals
          }
          return finish_integration(diagnostics, options) unless converged
          diagnostics
        end

        def finish_integration(diagnostics, options)
          return diagnostics if options[:diagnostics]
          message = 'stable distribution numerical integration did not reach requested precision'
          raise IntegrationError.new(message, diagnostics)
        end

        def gk15_integrate(left, right)
          center = (left + right) / 2.0
          half_width = (right - left) / 2.0
          kronrod = 0.0
          gauss = 0.0
          GK15_X.each_with_index do |node, index|
            value_left = safe_integrand(center - half_width * node) { |theta| yield(theta) }
            value_right = safe_integrand(center + half_width * node) { |theta| yield(theta) }
            sum = value_left + value_right
            kronrod += GK15_W[index] * sum
            gauss += GK7_W[index] * sum
          end
          center_value = safe_integrand(center) { |theta| yield(theta) }
          kronrod = (right - left) * 0.5 * (kronrod + GK15_W.last * center_value)
          gauss = (right - left) * 0.5 * (gauss + GK7_W.last * center_value)
          [kronrod, kronrod - gauss, 15]
        end

        def safe_integrand(theta)
          value = yield(theta)
          value.finite? ? value : 0.0
        rescue NoMethodError
          0.0
        rescue Math::DomainError
          0.0
        end
      end
    end
  end
end
