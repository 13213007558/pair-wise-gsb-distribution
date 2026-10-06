module Distribution
  module WalleniusNoncentralHypergeometric
    module Ruby_
      class << self
        def pdf(k, m1, m2, n, odds = 1)
          validate_draw_count(k)
          distribution(m1, m2, n, odds) do |probabilities, lower, upper|
            return 0 if k < lower || k > upper
            probabilities[k]
          end
        end

        alias_method :pmf, :pdf

        def cdf(k, m1, m2, n, odds = 1)
          validate_draw_count(k)
          distribution(m1, m2, n, odds) do |probabilities, lower, upper|
            return 0 if k < lower
            return 1 if k >= upper
            return Hypergeometric::Ruby_.cdf(k, m1, n, m1 + m2) if odds == 1
            sum = probabilities[lower..k].inject(0, :+)
            normalize_cdf(sum)
          end
        end

        def mean(m1, m2, n, odds = 1)
          validate_parameters(m1, m2, n)
          validate_odds(odds)
          probabilities = probability_mass(m1, m2, n, odds)
          lower, upper = support(m1, m2, n)
          (lower..upper).inject(0) do |total, color1_count|
            total + color1_count * probabilities[color1_count]
          end
        end

        alias_method :expectation, :mean
        alias_method :expected_value, :mean

        def exact_pdf(k, m1, m2, n, odds = 1)
          validate_draw_count(k)
          validate_parameters(m1, m2, n)
          validate_odds(odds)
          lower, upper = support(m1, m2, n)
          return Rational(0) if k < lower || k > upper

          return Hypergeometric::Ruby_.pdf(k, m1, n, m1 + m2) if odds == 1
          return Rational(1) if (odds.zero? && k == [0, n - m2].max) ||
                                (odds.infinite? && k == [m1, n].min)
          return Rational(0) if odds.zero? || odds.infinite?

          ratio = odds.to_r

          enumerate_paths(m1, m2, ratio, k, n - k)
        end

        alias_method :exact_pmf, :exact_pdf

        def pdf_with_status(k, m1, m2, n, odds = 1, tolerance = 1e-10, max_depth = 20)
          validate_draw_count(k)
          validate_parameters(m1, m2, n)
          validate_odds(odds)
          validate_integration_options(tolerance, max_depth)
          lower, upper = support(m1, m2, n)
          return [0.0, :outside_support] if k < lower || k > upper
          return [pdf(k, m1, m2, n, odds).to_f, :equal_weights] if odds == 1
          return [pdf(k, m1, m2, n, odds).to_f, :degenerate] if odds.zero? || odds.infinite?
          if odds.to_f > 1e100 || odds.to_f < 1e-300
            return [pdf(k, m1, m2, n, odds).to_f, :integral_not_evaluated]
          end

          integrate_pdf(k, m1, m2, n, odds.to_f, tolerance.to_f, max_depth.to_i)
        end

        alias_method :pmf_with_status, :pdf_with_status
        alias_method :integral_pdf, :pdf_with_status

        def rng(m1, m2, n, odds = 1, seed = nil, random: nil)
          prng = random || (seed.respond_to?(:rand) ? seed : Random.new(seed.nil? ? Random.new_seed : seed))
          -> { sample(m1, m2, n, odds, random: prng) }
        end

        def sample(m1, m2, n, odds = 1, random: Random)
          validate_parameters(m1, m2, n)
          validate_odds(odds)

          remaining1 = m1
          remaining2 = m2
          draws = n.times.count do
            probability = draw_probability(remaining1, remaining2, odds)
            draw = prng_draw(random) <= probability
            remaining1 -= 1 if draw
            remaining2 -= 1 unless draw
            draw
          end
          draws
        end

        def sample_draws(m1, m2, n, odds = 1, seed = nil, random: nil)
          prng = random || (seed.respond_to?(:rand) ? seed : Random.new(seed.nil? ? Random.new_seed : seed))
          remaining1 = m1
          remaining2 = m2

          Array.new(n) do
            probability = draw_probability(remaining1, remaining2, odds)
            draw = prng.rand <= probability
            remaining1 -= 1 if draw
            remaining2 -= 1 unless draw
            draw ? 1 : 0
          end
        end

        alias_method :sequential_sample, :sample_draws

        private

        def distribution(m1, m2, n, odds)
          validate_parameters(m1, m2, n)
          validate_odds(odds)
          lower, upper = support(m1, m2, n)
          probabilities = probability_mass(m1, m2, n, odds)
          yield(probabilities, lower, upper)
        end

        def validate_parameters(m1, m2, n)
          [m1, m2, n].each do |value|
            unless value.is_a?(Integer) && value >= 0
              fail(ArgumentError, 'population and draw counts must be non-negative integers')
            end
          end
          fail(ArgumentError, 'cannot draw more items than the population contains') if n > m1 + m2
        end

        def validate_draw_count(k)
          unless k.is_a?(Integer) && k >= 0
            fail(ArgumentError, 'k must be a non-negative integer')
          end
        end

        def validate_integration_options(tolerance, max_depth)
          unless tolerance.is_a?(Numeric) && tolerance.finite? && tolerance > 0 &&
                 max_depth.is_a?(Integer) && max_depth >= 0
            fail(ArgumentError, 'tolerance must be positive and max_depth must be a non-negative integer')
          end
        end

        def validate_odds(odds)
          unless odds.is_a?(Numeric) && odds >= 0 && (odds.infinite? || odds.finite?)
            fail(ArgumentError, 'odds must be a non-negative finite number or Infinity')
          end
        end

        def support(m1, m2, n)
          lower = [0, n - m2].max
          upper = [m1, n].min
          [lower, upper]
        end

        def probability_mass(m1, m2, n, odds)
          return exact_probability_mass(m1, m2, n) if odds == 1
          return degenerate_probability_mass(m1, m2, n, odds) if odds.zero? || odds.infinite?
          return exact_probability_mass_from_rational(m1, m2, n, odds) if odds.is_a?(Rational)

          log_probabilities = Array.new(m1 + 1, -Float::INFINITY)
          log_probabilities[0] = 0.0
          (0...n).each do |draws|
            next_log_probabilities = Array.new(m1 + 1, -Float::INFINITY)
            each_present_log_state(log_probabilities, draws, m1, m2) do |log_probability, type1_drawn|
              remaining1 = m1 - type1_drawn
              remaining2 = m2 - (draws - type1_drawn)
              p = draw_probability(remaining1, remaining2, odds)

              if remaining1 > 0
                add_log_probability(
                  next_log_probabilities,
                  type1_drawn + 1,
                  log_probability + Math.log(p)
                )
              end
              if remaining2 > 0
                add_log_probability(
                  next_log_probabilities,
                  type1_drawn,
                  log_probability + Math.log(1 - p)
                )
              end
            end
            log_probabilities = next_log_probabilities
          end

          normalize_log_mass(log_probabilities, m1, m2, n)
        end

        def exact_probability_mass_from_rational(m1, m2, n, odds)
          lower, upper = support(m1, m2, n)
          probabilities = Array.new(m1 + 1)
          (lower..upper).each do |k|
            probabilities[k] = enumerate_paths(m1, m2, odds, k, n - k)
          end
          probabilities
        end

        def degenerate_probability_mass(m1, m2, n, odds)
          probabilities = Array.new(m1 + 1)
          lower, upper = support(m1, m2, n)
          (lower..upper).each { |k| probabilities[k] = 0.0 }
          forced_k = odds.zero? ? [0, n - m2].max : [m1, n].min
          probabilities[forced_k] = 1.0
          probabilities
        end

        def exact_probability_mass(m1, m2, n)
          lower, upper = support(m1, m2, n)
          probabilities = Array.new(m1 + 1)
          (lower..upper).each do |k|
            probabilities[k] = Hypergeometric::Ruby_.pdf(k, m1, n, m1 + m2)
          end
          probabilities
        end

        def each_present_log_state(log_probabilities, draws, m1, m2)
          lower = [0, draws - m2].max
          upper = [m1, draws].min
          (lower..upper).each do |type1_drawn|
            log_probability = log_probabilities[type1_drawn]
            unless log_probability.nil? || log_probability == -Float::INFINITY
              yield(log_probability, type1_drawn)
            end
          end
        end

        def add_log_probability(log_probabilities, k, value)
          if log_probabilities[k] == -Float::INFINITY
            log_probabilities[k] = value
          else
            log_probabilities[k] = log_add_exp(log_probabilities[k], value)
          end
        end

        def log_add_exp(first, second)
          return second if first == -Float::INFINITY
          return first if second == -Float::INFINITY
          larger, smaller = [first, second].sort.reverse
          larger + Math.log(1 + Math.exp(smaller - larger))
        end

        def normalize_log_mass(log_probabilities, m1, m2, n)
          lower, upper = support(m1, m2, n)
          log_total = (lower..upper).map { |k| log_probabilities[k] || -Float::INFINITY }
                              .inject(-Float::INFINITY) { |total, value| log_add_exp(total, value) }
          (lower..upper).each do |k|
            log_probabilities[k] = 0.0 if log_total == -Float::INFINITY && k == upper
            log_probabilities[k] ||= -Float::INFINITY
            log_probabilities[k] = Math.exp(log_probabilities[k] - log_total) if log_total != -Float::INFINITY
          end
          log_probabilities
        end

        def normalize_cdf(value)
          return 0 if value.negative?
          return 1 if value > 1
          value
        end

        def draw_probability(remaining1, remaining2, odds)
          return 1.0 if remaining2.zero?
          return 0.0 if remaining1.zero?
          return remaining1.quo(remaining1 + remaining2).to_f if odds == 1

          if odds.zero?
            0.0
          elsif odds.infinite?
            1.0
          elsif odds.is_a?(Rational)
            Rational(remaining1) * odds / (Rational(remaining1) * odds + Rational(remaining2))
          else
            log_odds = Math.log(odds.to_f)
            log_relative_weight = log_odds + Math.log(remaining1.to_f) - Math.log(remaining2.to_f)
            return 1.0 if log_relative_weight > 709
            return 0.0 if log_relative_weight < -709
            if log_relative_weight >= 0
              1 / (1 + Math.exp(-log_relative_weight))
            else
              exponential = Math.exp(log_relative_weight)
              exponential / (1 + exponential)
            end
          end
        end

        def prng_draw(random)
          random.respond_to?(:rand) ? random.rand : Random.rand
        end

        def enumerate_paths(m1, m2, odds, type1_drawn, type2_drawn)
          return Rational(1) if type1_drawn.zero? && type2_drawn.zero?

          remaining1 = m1 - type1_drawn
          remaining2 = m2 - type2_drawn
          probability = Rational(0)

          if type1_drawn.positive?
            before1 = remaining1 + 1
            step = (Rational(before1) * odds).quo(Rational(before1) * odds + Rational(remaining2))
            probability += step * enumerate_paths(m1, m2, odds, type1_drawn - 1, type2_drawn)
          end

          if type2_drawn.positive?
            before2 = remaining2 + 1
            step = Rational(before2).quo(Rational(remaining1) * odds + Rational(before2))
            probability += step * enumerate_paths(m1, m2, odds, type1_drawn, type2_drawn - 1)
          end

          probability
        end

        def integrate_pdf(k, m1, m2, n, odds, tolerance, max_depth)
          ways = Math.binomial_coefficient(m1, k) * Math.binomial_coefficient(m2, n - k)
          remaining1 = m1 - k
          remaining2 = m2 - n + k

          if k.positive?
            first, first_status = adaptive_integral(
              0, 1,
              ways * odds * k,
              odds * (remaining1 + 1) + remaining2 - 1,
              k - 1,
              n - k,
              odds,
              tolerance,
              max_depth
            )
          else
            first = 0.0
            first_status = :ok
          end

          second_draws = n - k
          if second_draws.positive?
            second, second_status = adaptive_integral(
              0, 1,
              ways * (n - k),
              odds * remaining1 + remaining2,
              k,
              n - k - 1,
              odds,
              tolerance,
              max_depth
            )
          else
            second = 0.0
            second_status = :ok
          end

          value = first + second
          status = [first_status, second_status].include?(:did_not_converge) ? :did_not_converge : :ok
          [[value, 0].max, status]
        end

        def adaptive_integral(lower, upper, constant, exponent_x, exponent_y1, exponent_y2, odds, tolerance, max_depth)
          if exponent_x >= 0
            f = method(:integrand).to_proc.curry[
              constant, exponent_x, exponent_y1, exponent_y2, odds
            ]
          else
            scale = exponent_x + 1
            f = ->(y) do
              x = y**(1 / scale)
              (constant / scale) * (1 - x**odds)**exponent_y1 * (1 - x)**exponent_y2
            end
          end
          middle = (lower + upper) / 2.0
          whole = simpson(f, lower, middle, upper)
          recursive_integral(f, lower, middle, upper, whole, tolerance, max_depth)
        end

        def recursive_integral(f, lower, middle, upper, whole, tolerance, depth)
          left_middle = (lower + middle) / 2.0
          right_middle = (middle + upper) / 2.0
          left = simpson(f, lower, left_middle, middle)
          right = simpson(f, middle, right_middle, upper)
          refined = left + right
          error = (refined - whole).abs

          if depth.zero? || (!error.nan? && error <= 15 * tolerance)
            [refined, depth.zero? && error > tolerance ? :did_not_converge : :ok]
          else
            left_value, left_status = recursive_integral(
              f, lower, left_middle, middle, left, tolerance / 2, depth - 1
            )
            right_value, right_status = recursive_integral(
              f, middle, right_middle, upper, right, tolerance / 2, depth - 1
            )
            status = [left_status, right_status].include?(:did_not_converge) ? :did_not_converge : :ok
            [left_value + right_value, status]
          end
        end

        def simpson(f, lower, middle, upper)
          (upper - lower) / 6.0 * (f.call(lower) + 4 * f.call(middle) + f.call(upper))
        end

        def integrand(constant, exponent_x, exponent_y1, exponent_y2, odds, x)
          constant * x**exponent_x * (1 - x**odds)**exponent_y1 * (1 - x)**exponent_y2
        end
      end
    end
  end
end
