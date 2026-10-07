module Distribution
  module Tweedie
    module Ruby_
      class << self
        include Math

        # Probability density for strictly positive observations.
        #
        # The mixed distribution's probability mass at zero is available
        # from +zero_mass+ or +pmf+.
        def pdf(x, mean, dispersion, power, tolerance: 1e-12)
          validate_parameters(mean, dispersion, power)
          return 0.0 unless x.positive?
          exp(log_pdf(x, mean, dispersion, power, tolerance: tolerance))
        end

        def log_pdf(x, mean, dispersion, power, tolerance: 1e-12)
          validate_parameters(mean, dispersion, power)
          return -Float::INFINITY unless x.positive?

          x = x.to_f
          unless tolerance.is_a?(Numeric) && tolerance > 0 && tolerance < 1
            fail ArgumentError, 'tolerance must be between 0 and 1'
          end

          alpha_shape = (2.0 - power) / (power - 1.0)
          theta_scale = (power - 1.0) * dispersion * mean**(power - 1.0)
          poisson_mean = mean**(2.0 - power) / ((2.0 - power) * dispersion)

          constant = -poisson_mean - x / theta_scale - log(x)
          slope = log(poisson_mean) - alpha_shape * log(theta_scale) +
                  alpha_shape * log(x)

          log_term = lambda do |n|
            constant + n * slope - lgamma(n + 1.0).first -
              lgamma(n * alpha_shape).first
          end

          mode = approximate_mode(poisson_mean, theta_scale, x, alpha_shape)
          mode += 1 while log_term.call(mode + 1) > log_term.call(mode)
          mode -= 1 while mode > 1 && log_term.call(mode - 1) > log_term.call(mode)

          log_series(log_term, mode, tolerance)
        end

        def zero_mass(mean, dispersion, power)
          validate_parameters(mean, dispersion, power)
          poisson_mean = mean**(2.0 - power) / ((2.0 - power) * dispersion)
          exp(-poisson_mean)
        end

        def log_zero_mass(mean, dispersion, power)
          validate_parameters(mean, dispersion, power)
          -mean**(2.0 - power) / ((2.0 - power) * dispersion)
        end

        def pmf(x, mean, dispersion, power)
          x == 0 ? zero_mass(mean, dispersion, power) : 0.0
        end

        alias_method :probability_zero, :zero_mass

        def rng(mean, dispersion, power, seed = nil, random: nil)
          validate_parameters(mean, dispersion, power)
          random ||= Random.new(seed.nil? ? Random.new_seed : seed)
          alpha_shape = (2.0 - power) / (power - 1.0)
          theta_scale = (power - 1.0) * dispersion * mean**(power - 1.0)
          poisson_mean = mean**(2.0 - power) / ((2.0 - power) * dispersion)

          lambda do
            count = poisson_sample(random, poisson_mean)
            next 0.0 if count.zero?

            gamma_shape = count * alpha_shape
            if gamma_shape < 1.0
              gamma_shape += 1.0
              booster = random.rand**(1.0 / (gamma_shape - 1.0))
            else
              booster = 1.0
            end

            normal = standard_normal(random)
            intermediate = gamma_shape - 1.0 / 3.0
            inner = 1.0 / sqrt(9.0 * intermediate)
            gamma_variate = (1.0 + inner * normal)**3

            while gamma_variate <= 0.0 ||
                  log(random.rand) >=
                    0.5 * normal**2 + intermediate -
                    intermediate * gamma_variate +
                    intermediate * log(gamma_variate)
              normal = standard_normal(random)
              gamma_variate = (1.0 + inner * normal)**3
            end

            booster * intermediate * gamma_variate * theta_scale
          end
        end

        private

        def validate_parameters(mean, dispersion, power)
          unless valid_positive?(mean) && valid_positive?(dispersion) &&
                 power.is_a?(Numeric) && power.finite? && power > 1 && power < 2
            fail ArgumentError,
                 'mean and dispersion must be positive and power must be between 1 and 2'
         end
        end

        def valid_positive?(value)
          value.is_a?(Numeric) && value.finite? && value > 0
        end

        def approximate_mode(poisson_mean, theta_scale, x, alpha_shape)
          mode = x * (poisson_mean / theta_scale**alpha_shape)**
                 (1.0 / (1.0 + alpha_shape))
          [1, mode.round.to_i].max
        end

        def log_series(log_term, mode, tolerance)
          maximum = log_term.call(mode)
          total = 1.0

          upper = mode
          lower = mode
          loop do
            upper += 1
            term = exp(log_term.call(upper) - maximum)
            total += term
            ratio = exp(log_term.call(upper) - log_term.call(upper - 1))
            next if ratio >= 1.0

            upper_done = term * ratio / (1.0 - ratio) <= tolerance * total
            lower_done = lower == 1
            unless lower_done
              lower -= 1
              term = exp(log_term.call(lower) - maximum)
              total += term
              ratio = exp(log_term.call(lower) - log_term.call(lower + 1))
              if ratio < 1.0 && term * ratio / (1.0 - ratio) <= tolerance * total
                lower_done = true
              end
              lower = 1 if lower < 1
            end

            break maximum + log(total) if upper_done && lower_done
          end
        end

        def poisson_sample(random, poisson_mean)
          return poisson_knuth(random, poisson_mean) if poisson_mean <= 30.0

          c = 0.767 - 3.36 / poisson_mean
          beta = PI / sqrt(3.0 * poisson_mean)
          alpha = beta * poisson_mean

          loop do
            u = random.rand
            v = random.rand
            x = (alpha - log((1.0 - u) / u)) / beta
            n = (x + 0.5).floor
            next if n < 0

            log_envelope = alpha - beta * x - log(v / (1.0 - v)**2)
            log_target = n * log(poisson_mean) - poisson_mean +
                         lgamma(n + 1.0).first
            log_envelope += log(c) - log(beta)
            if log_envelope <= log_target
              return n
            end
          end
        end

        def poisson_knuth(random, poisson_mean)
          limit = exp(-poisson_mean)
          product = 1.0
          count = 0
          while product > limit
            product *= random.rand
            count += 1
          end
          count - 1
        end

        def standard_normal(random)
          radius = sqrt(-2.0 * log(random.rand))
          angle = 2.0 * PI * random.rand
          radius * cos(angle)
        end
      end
    end
  end
end
