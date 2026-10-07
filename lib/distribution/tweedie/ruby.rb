module Distribution
  module Tweedie
    module Ruby_
      class << self
        include Math

        DEFAULT_TOLERANCE = 1e-12
        def pdf(x, mean, dispersion, power, tolerance = DEFAULT_TOLERANCE)
          x = x.to_f
          validate_common_parameters!(mean.to_f, dispersion.to_f, power.to_f)
          return 0.0 unless x.positive?

          exp(log_pdf(x, mean, dispersion, power, tolerance))
        end

        def log_pdf(x, mean, dispersion, power, tolerance = DEFAULT_TOLERANCE)
          x = x.to_f
          mean = mean.to_f
          dispersion = dispersion.to_f
          power = power.to_f
          validate_common_parameters!(mean, dispersion, power)
          return -Float::INFINITY unless x.positive?
          validate_tolerance!(tolerance)

          parameters = parameters(mean, dispersion, power)
          log_weight_series(x, parameters, tolerance) + linear_log_density(x, parameters)
        end

        def zero_probability(mean, dispersion, power)
          mean = mean.to_f
          dispersion = dispersion.to_f
          power = power.to_f
          validate_common_parameters!(mean, dispersion, power)

          exp(-parameters(mean, dispersion, power)[:poisson_mean])
        end

        def rng(mean, dispersion, power, seed = nil, random: nil)
          random = seed[:random] if seed.is_a?(Hash) && random.nil?
          mean = mean.to_f
          dispersion = dispersion.to_f
          power = power.to_f
          validate_common_parameters!(mean, dispersion, power)

          parameters = parameters(mean, dispersion, power)
          prng = random || random_generator(seed)

          lambda do
            count = poisson_sample(parameters[:poisson_mean], prng)
            next 0.0 if count.zero?

            shape = count * parameters[:gamma_shape]
            gamma_sample(shape, parameters[:gamma_scale], prng)
          end
        end

        private

        def validate_tolerance!(tolerance)
          unless tolerance.is_a?(Numeric) && tolerance > 0 && tolerance < 1
            fail ArgumentError, 'tolerance must be between 0 and 1'
          end
        end

        def validate_common_parameters!(mean, dispersion, power)
          unless mean.positive?
            fail ArgumentError, 'mean must be positive'
          end
          unless dispersion.positive?
            fail ArgumentError, 'dispersion must be positive'
          end
          unless power > 1.0 && power < 2.0
            fail ArgumentError, 'power must satisfy 1 < power < 2'
          end
        end

        def parameters(mean, dispersion, power)
          poisson_mean = mean**(2.0 - power) / (dispersion * (2.0 - power))
          gamma_shape = (2.0 - power) / (power - 1.0)
          gamma_scale = dispersion * (power - 1.0) * mean**(power - 1.0)

          {
            poisson_mean: poisson_mean,
            gamma_shape: gamma_shape,
            gamma_scale: gamma_scale
          }
        end

        def linear_log_density(x, parameters)
          scale = parameters[:gamma_scale]
          poisson_mean = parameters[:poisson_mean]

          -log(x) - poisson_mean - x / scale
        end

        def log_weight_series(x, parameters, tolerance)
          scale = parameters[:gamma_scale]
          shape = parameters[:gamma_shape]
          z = parameters[:poisson_mean] * (x / scale)**shape
          return -Float::INFINITY if z.zero?

          mode = series_mode(z, shape)

          log_weight = lambda do |count|
            count * log(z) - loggamma(count + 1.0) - loggamma(shape * count)
          end

          mode_log_weight = log_weight.call(mode)
          sum = 1.0
          left = mode
          right = mode
          left_done = mode == 1
          right_done = false
          stop = tolerance.fdiv(4.0)

          until left_done && right_done
            unless right_done
              q = tail_ratio(z, shape, right)
              if q < 1.0
                right_bound = q.fdiv(1.0 - q)
                right_done = right_bound <= stop * sum
              end

              unless right_done
                right += 1
                term = exp(log_weight.call(right) - mode_log_weight)
                sum += term
                right_done = term.zero?
              end
            end

            unless left_done
              left -= 1
              if left.zero?
                left_done = true
              else
                term = exp(log_weight.call(left) - mode_log_weight)
                sum += term
                backward_ratio = if left == 1
                                   0.0
                                 else
                                   forward_ratio = tail_ratio(z, shape, left - 1)
                                   forward_ratio.zero? ? 0.0 : 1.0 / forward_ratio
                                 end
                left_done = backward_ratio <= stop || term.zero?
                left_bound = term.zero? ? 0.0 : term * backward_ratio.fdiv(1.0 - backward_ratio)
                left_done = true if left_bound <= stop * sum
              end
            end
          end

          log(sum) + mode_log_weight
        end

        def series_mode(z, shape)
          return 1 if tail_ratio(z, shape, 1) <= 1.0

          estimate = (z**(1.0 / (shape + 1.0))) /
                     shape**(shape / (shape + 1.0))
          upper = [2, Integer(estimate + 2.0)].max
          upper *= 2 while tail_ratio(z, shape, upper) > 1.0

          lower = upper / 2
          lower = 1 if lower < 1 || tail_ratio(z, shape, lower) <= 1.0

          while upper - lower > 1
            middle = (lower + upper) / 2
            if tail_ratio(z, shape, middle) > 1.0
              lower = middle
            else
              upper = middle
            end
          end

          upper
        end

        def tail_ratio(z, shape, count)
          return 0.0 if z.zero?

          next_count = count + 1.0
          log_ratio = log(z) - log(next_count) -
            (loggamma(shape * next_count) - loggamma(shape * count))
          return 0.0 if log_ratio <= -745.0
          return Float::INFINITY if log_ratio >= 709.0

          exp(log_ratio)
        end

        def poisson_sample(mean, prng)
          return poisson_sample_small(mean, prng) if mean <= 30.0

          radius = Integer(10.0 * sqrt(mean) + 6.0)
          mode = Integer(mean)
          lower = [0, mode - radius].max
          upper = mode + radius
          log_weights = (lower..upper).map do |count|
            count * log(mean) - loggamma(count + 1.0) - mean
          end
          max_log_weight = log_weights.max
          weights = log_weights.map { |log_weight| exp(log_weight - max_log_weight) }
          total_weight = weights.sum
          target = prng.rand * total_weight
          cumulative = 0.0
          weights.each_with_index do |_weight, index|
            cumulative += weights[index]
            return lower + index if target <= cumulative
          end

          upper
        end

        def poisson_sample_small(mean, prng)
          limit = exp(-mean)
          product = prng.rand
          count = 0

          while product > limit
            count += 1
            product *= prng.rand
          end

          count
        end

        def gamma_sample(shape, scale, prng)
          boost = false
          if shape < 1.0
            shape += 1.0
            boost = true
          end

          d = shape - 1.0 / 3.0
          c = 1.0 / (3.0 * sqrt(d))

          value = loop do
            normal = standard_normal(prng)
            temporary = (1.0 + c * normal)**3
            uniform = prng.rand

            if temporary > 0.0 &&
               log(uniform) < 0.5 * normal**2 + d * (1.0 - temporary + log(temporary))
              break d * temporary
            end
          end

          value *= prng.rand**(1.0 / (shape - 1.0)) if boost
          value * scale
        end

        def standard_normal(prng)
          sqrt(-2.0 * log(prng.rand)) * cos(2.0 * PI * prng.rand)
        end

        def random_generator(seed)
          return seed if seed.is_a?(Random)

          Random.new(seed.nil? ? Random.new_seed : seed)
        end
      end
    end
  end
end
