module Distribution
  class NoncentralFConvergenceError < StandardError
    attr_reader :result, :error, :iterations, :budget

    def initialize(message, result:, error:, iterations:, budget:)
      @result = result
      @error = error
      @iterations = iterations
      @budget = budget
      super(message)
    end
  end

  module NoncentralF
    module Ruby_
      DEFAULT_TOLERANCE = 1e-12
      DEFAULT_ITERATIONS = 50_000

      class << self
        def pdf(x, d1, d2, ncp, tolerance: DEFAULT_TOLERANCE,
                             max_iterations: DEFAULT_ITERATIONS)
          x = Float(x)
          d1, d2, ncp = validate_parameters!(d1, d2, ncp, tolerance, max_iterations)
          return central_f(:pdf, x, d1, d2) if ncp.zero?
          return 0.0 if x.infinite? && x.positive?
          return boundary_pdf(x, d1, d2, ncp) if x <= 0.0

          value, error, iterations = pdf_series(x, d1, d2, ncp, tolerance, max_iterations)
          ensure_convergence!('PDF', value, error, iterations, max_iterations, tolerance)
          value
        end

        def cdf(x, d1, d2, ncp, tolerance: DEFAULT_TOLERANCE,
                             max_iterations: DEFAULT_ITERATIONS)
          x = Float(x)
          d1, d2, ncp = validate_parameters!(d1, d2, ncp, tolerance, max_iterations)
          return 0.0 if x <= 0.0
          return 1.0 if x.infinite? && x.positive?
          return central_f(:cdf, x, d1, d2) if ncp.zero?

          lower, lower_error, lower_iterations = probability_series(:cdf, x, d1, d2, ncp,
                                                                   tolerance, max_iterations)
          if lower > 0.5
            upper, upper_error, upper_iterations = probability_series(:survival, x, d1, d2, ncp,
                                                                      tolerance, max_iterations)
            ensure_convergence!('CDF', 1.0 - upper, upper_error,
                                lower_iterations + upper_iterations, max_iterations, tolerance)
            1.0 - upper
          else
            ensure_convergence!('CDF', lower, lower_error, lower_iterations,
                                max_iterations, tolerance)
            lower
          end
        end

        def survival(x, d1, d2, ncp, tolerance: DEFAULT_TOLERANCE,
                                max_iterations: DEFAULT_ITERATIONS)
          x = Float(x)
          d1, d2, ncp = validate_parameters!(d1, d2, ncp, tolerance, max_iterations)
          return 1.0 if x <= 0.0
          return 0.0 if x.infinite? && x.positive?
          return central_f(:q, x, d1, d2) if ncp.zero?

          upper, upper_error, upper_iterations = probability_series(:survival, x, d1, d2, ncp,
                                                                    tolerance, max_iterations)
          if upper > 0.5
            lower, lower_error, lower_iterations = probability_series(:cdf, x, d1, d2, ncp,
                                                                     tolerance, max_iterations)
            ensure_convergence!('survival', 1.0 - lower, lower_error,
                                lower_iterations + upper_iterations, max_iterations, tolerance)
            1.0 - lower
          else
            ensure_convergence!('survival', upper, upper_error, upper_iterations,
                                max_iterations, tolerance)
            upper
          end
        end

        alias_method :q, :survival

        def quantile(probability, d1, d2, ncp, tolerance: DEFAULT_TOLERANCE,
                                        max_iterations: DEFAULT_ITERATIONS)
          p = Float(probability)
          d1, d2, ncp = validate_parameters!(d1, d2, ncp, tolerance, max_iterations)
          raise ArgumentError, 'probability must be in [0, 1]' unless p.between?(0.0, 1.0)
          return 0.0 if p.zero?
          return Float::INFINITY if p == 1.0
          return central_f(:quantile, p, d1, d2) if ncp.zero?

          inverse_cdf(p, d1, d2, ncp, tolerance, max_iterations)
        end

        alias_method :p_value, :quantile

        private

        def validate_parameters!(d1, d2, ncp, tolerance, max_iterations)
          d1 = Float(d1)
          d2 = Float(d2)
          ncp = Float(ncp)
          unless d1.finite? && d1 > 0.0 && d2.finite? && d2 > 0.0
            raise ArgumentError, 'degrees of freedom must be positive and finite'
          end
          unless ncp.finite? && ncp >= 0.0
            raise ArgumentError, 'noncentrality parameter must be finite and non-negative'
          end
          unless tolerance.is_a?(Numeric) && tolerance > 0 && Float(tolerance).finite?
            raise ArgumentError, 'tolerance must be a positive finite number'
          end
          unless max_iterations.is_a?(Integer) && max_iterations >= 1
            raise ArgumentError, 'max_iterations must be a positive Integer'
          end
          [d1, d2, ncp]
        end

        def central_f(method, *arguments)
          Distribution::F::Ruby_.public_send(method, *arguments)
        end

        def boundary_pdf(x, d1, d2, ncp)
          return 0.0 if x < 0.0
          return Float::INFINITY if d1 < 2.0
          return 0.0 if d1 > 2.0
          Math.exp(-ncp / 2.0) * d1 / d2 / Math.beta(d1 / 2.0, d2 / 2.0)
        end

        def beta_cdf(y, a, b)
          return 0.0 if y <= 0.0
          return 1.0 if y >= 1.0
          [0.0, [1.0, Math.regularized_beta(y, a, b)].min].max
        end

        def beta_sf(y, a, b)
          return 1.0 if y <= 0.0
          return 0.0 if y >= 1.0
          return 1.0 - beta_cdf(y, a, b) if y < 0.999999
          z = 1.0 - y
          log_prefactor = -Math.logbeta(a, b) + b * Math.log(z) + a * Math.log(y)
          prefactor = y * Math.exp(log_prefactor)
          return 0.0 if prefactor.zero?
          value = prefactor * Math::IncompleteBeta.continued_fraction(b, a, z) / b
          [0.0, [1.0, value].min].max
        end

        def beta_value(kind, y, a, b)
          kind == :cdf ? beta_cdf(y, a, b) : beta_sf(y, a, b)
        end

        def poisson_mode(mu)
          mu < 1.0 ? 0 : mu.floor
        end

        def log_poisson_weight(j, mu)
          -mu + j * Math.log(mu) - Math.lgamma(j + 1).first
        end

        def central_log_pdf(x, d1, d2, j)
          y = d1 * x / (d1 * x + d2)
          a = d1 / 2.0 + j
          b = d2 / 2.0
          -Math.logbeta(a, b) + (a - 1.0) * Math.log(y) + (b - 1.0) * Math.log(1.0 - y) +
            Math.log(d1) + Math.log(d2) - 2.0 * Math.log(d1 * x + d2)
        end

        def pdf_up_ratio(j, mu, x, d1, d2, a0)
          a = a0 + j
          denominator = d1 * x + d2
          y = d1 * x / denominator
          mu / (j + 1) * y * (a + d2 / 2.0) / a
        end

        def pdf_down_ratio(j, mu, x, d1, d2, a0)
          1.0 / pdf_up_ratio(j - 1, mu, x, d1, d2, a0)
        end

        def pdf_mode(mode, mu, x, d1, d2, a0, budget)
          j = mode
          budget.times do
            up = pdf_up_ratio(j, mu, x, d1, d2, a0)
            down = j.zero? ? 0.0 : pdf_down_ratio(j, mu, x, d1, d2, a0)
            if up >= 1.0
              j += 1
            elsif down > 1.0
              j -= 1
            else
              return j
            end
          end
          raise NoncentralFConvergenceError.new(
            'PDF mode search did not converge',
            result: nil, error: Float::INFINITY, iterations: budget, budget: budget
          )
        end

        def pdf_series(x, d1, d2, ncp, tolerance, budget)
          mu = Float(ncp) / 2.0
          a0 = d1 / 2.0
          mode = pdf_mode(poisson_mode(mu), mu, x, d1, d2, a0, budget)
          term = Math.exp(log_poisson_weight(mode, mu) + central_log_pdf(x, d1, d2, mode))
          value = term
          iterations = 1
          errors = [0.0, 0.0]

          [[:up, mode], [:down, mode]].each_with_index do |(direction, start), side|
            if direction == :down
              break if start.zero?
              j = start - 1
              current_term = term * pdf_down_ratio(start, mu, x, d1, d2, a0)
              value += current_term
              iterations += 1 unless iterations >= budget
              future_ratio = pdf_future_ratio(:down, j, mu, x, d1, d2, a0)
              errors[side] = geometric_bound(current_term, future_ratio)
              break if errors[side] <= tolerance / 2.0 && errors[side] <= value.abs * tolerance
            else
              j = start
              current_term = term
            end
            loop do
              break if iterations >= budget
              next_j = direction == :up ? j + 1 : j - 1
              break if next_j < 0
              ratio = if direction == :up
                        pdf_up_ratio(j, mu, x, d1, d2, a0)
                      else
                        pdf_down_ratio(j, mu, x, d1, d2, a0)
                      end
              current_term *= ratio
              value += current_term
              j = next_j
              iterations += 1 unless iterations >= budget
              future_index = direction == :up ? j : next_j
              future_ratio = next_j.zero? ? 0.0 : pdf_future_ratio(direction, future_index, mu, x, d1, d2, a0)
              errors[side] = geometric_bound(current_term, future_ratio)
              break if errors[side] <= tolerance / 2.0 && errors[side] <= value.abs * tolerance
            end
          end

          [value, errors.sum, iterations]
        end

        def pdf_future_ratio(direction, j, mu, x, d1, d2, a0)
          return 0.0 if direction == :down && j.zero?
          if direction == :up
            pdf_up_ratio(j, mu, x, d1, d2, a0)
          else
            pdf_down_ratio(j, mu, x, d1, d2, a0)
          end
        end

        def probability_series(kind, x, d1, d2, ncp, tolerance, budget)
          mu = Float(ncp) / 2.0
          a0 = d1 / 2.0
          b0 = d2 / 2.0
          y = d1 * x / (d1 * x + d2)
          mode = poisson_mode(mu)
          value = beta_value(kind, y, a0 + mode, b0)
          partition = 1.0
          iterations = 1
          errors = [0.0, 0.0]

          [[:up, mode], [:down, mode]].each_with_index do |(direction, start), side|
            if direction == :down
              break if start.zero?
              j = start - 1
              weight = start / mu
              value += weight * beta_value(kind, y, a0 + j, b0)
              partition += weight
              iterations += 1 unless iterations >= budget
              future_ratio = j.zero? ? 0.0 : j / mu
              errors[side] = geometric_bound(weight, future_ratio)
              break if errors[side] <= tolerance / 2.0
            else
              j = start
              weight = 1.0
            end
            loop do
              break if iterations >= budget
              next_j = direction == :up ? j + 1 : j - 1
              break if next_j < 0
              ratio = direction == :up ? mu / (j + 1) : j / mu
              weight *= ratio
              value += weight * beta_value(kind, y, a0 + next_j, b0)
              partition += weight
              j = next_j
              iterations += 1 unless iterations >= budget
              future_ratio = direction == :up ? mu / (j + 1) : j.zero? ? 0.0 : j / mu
              errors[side] = geometric_bound(weight, future_ratio)
              break if errors[side] <= tolerance / 2.0
            end
          end

          [value / partition, errors.sum / partition, iterations]
        end

        def geometric_bound(current_value, ratio)
          return 0.0 if current_value.zero? || ratio <= 0.0
          return Float::INFINITY if ratio >= 1.0
          current_value * ratio / (1.0 - ratio)
        end

        def ensure_convergence!(name, value, error, iterations, budget, tolerance)
          return if error.finite? && error <= tolerance
          raise NoncentralFConvergenceError.new(
            "#{name} did not reach tolerance #{tolerance} after #{iterations} terms",
            result: value, error: error, iterations: iterations, budget: budget
          )
        end

        def inverse_cdf(p, d1, d2, ncp, tolerance, budget)
          use_survival = p > 0.5
          target = use_survival ? 1.0 - p : p
          series_tolerance = [tolerance / 4.0, 1e-14].min
          lower = Float::EPSILON
          upper = initial_quantile(d1, d2, ncp)
          iterations = 0

          if use_survival
            while survival(upper, d1, d2, ncp, tolerance: series_tolerance,
                                      max_iterations: budget) > target
              iterations += 1 unless iterations >= budget
              raise_quantile_error!(lower, upper, iterations, budget, tolerance) if iterations >= budget
              lower = upper
              upper *= 2.0
            end
          else
            while cdf(upper, d1, d2, ncp, tolerance: series_tolerance,
                                  max_iterations: budget) < target
              iterations += 1
              raise_quantile_error!(lower, upper, iterations, budget, tolerance) if iterations >= budget
              upper *= 2.0
            end
          end

          loop do
            middle = Math.exp((Math.log(lower) + Math.log(upper)) / 2.0)
            probability = if use_survival
                            survival(middle, d1, d2, ncp, tolerance: series_tolerance,
                                                       max_iterations: budget)
                          else
                            cdf(middle, d1, d2, ncp, tolerance: series_tolerance,
                                                        max_iterations: budget)
                          end
            iterations += 1
            if probability > target
              if use_survival
                lower = middle
              else
                upper = middle
              end
            else
              if use_survival
                upper = middle
              else
                lower = middle
              end
            end
            log_error = Math.log(upper) - Math.log(lower)
            if log_error <= tolerance
              return newton_quantile(upper, target, d1, d2, ncp, use_survival,
                                     series_tolerance, budget, iterations)
            end
            raise_quantile_error!(lower, upper, iterations, budget, tolerance) if iterations >= budget
          end
        end

        def newton_quantile(x, target, d1, d2, ncp, use_survival,
                            series_tolerance, budget, iterations)
          4.times do
            probability = if use_survival
                            survival(x, d1, d2, ncp, tolerance: series_tolerance,
                                               max_iterations: budget)
                          else
                            cdf(x, d1, d2, ncp, tolerance: series_tolerance,
                                                max_iterations: budget)
                          end
            density = pdf(x, d1, d2, ncp, tolerance: series_tolerance,
                              max_iterations: budget)
            return x unless density.positive? && density.finite?
            candidate = x + (probability - target) / density
            return x unless candidate.positive? && candidate.finite?
            x = candidate
            iterations += 1
          end
          x
        end

        def initial_quantile(d1, d2, ncp)
          initial = 1.0 + ncp / d1
          return initial if d2 <= 2.0
          [initial, (d2 / (d2 - 2.0)) * (d1 + ncp) / d1].max
        end

        def raise_quantile_error!(lower, upper, iterations, budget, tolerance)
          raise NoncentralFConvergenceError.new(
            "inverse CDF did not reach tolerance #{tolerance} after #{iterations} iterations",
            result: upper, error: (upper - lower).abs, iterations: iterations, budget: budget
          )
        end
      end
    end
  end
end
