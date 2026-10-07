module Distribution
  class TruncatedNormal
    SQRT_2 = Math.sqrt(2.0)
    SQRT_2PI = Math.sqrt(2.0 * Math::PI)
    LOG_2PI = Math.log(2.0 * Math::PI)

    HAZARD_COEFFICIENTS = [1.0, 1.0, -2.0, 10.0, -74.0, 706.0,
                           -8162.0, 110_410.0, -1_708_394.0,
                           29_752_066.0].freeze
    TAIL_VARIANCE_COEFFICIENTS = [0.0, 1.0, -6.0, 50.0, -518.0,
                                  6354.0, -89_782.0, 1_435_330.0,
                                  -25_625_910.0].freeze

    attr_reader :mu, :sigma, :lower, :upper, :mean, :variance, :sd,
                :random, :attempts

    def initialize(mu: 0.0, sigma: 1.0, lower: -Float::INFINITY,
                   upper: Float::INFINITY, random: Random)
      @mu = mu.to_f
      @sigma = sigma.to_f
      @lower = lower.nil? ? -Float::INFINITY : lower.to_f
      @upper = upper.nil? ? Float::INFINITY : upper.to_f
      @random = random
      @attempts = 0

      raise ArgumentError, "sigma must be positive" unless @sigma.positive?
      raise ArgumentError, "lower must be less than upper" unless @lower < @upper

      if !@lower.finite? && !@upper.finite?
        @log_z = 0.0
        @mean = @mu
        @variance = @sigma * @sigma
      elsif @lower.finite? && @upper.finite?
        initialize_bounded
      elsif @lower.finite?
        initialize_lower
      else
        initialize_upper
      end
      @sd = Math.sqrt(@variance)
    end

    def pdf(x)
      return 0.0 if x < @lower || x > @upper

      z = (x - @mu) / @sigma
      Math.exp(log_phi(z) - @log_z) / @sigma
    end

    def cdf(x)
      return 0.0 if x <= @lower
      return 1.0 if x >= @upper

      z = (x - @mu) / @sigma
      if @lower.finite? && @upper.finite?
        bounded_cdf(z)
      elsif @lower.finite?
        1.0 - Math.exp(log_sf(z) - log_sf(@lower_z))
      elsif @upper.finite?
        Math.exp(log_cdf(z) - log_cdf(@upper_z))
      else
        normal_cdf(z)
      end
    end

    def quantile(probability)
      check_probability!(probability)
      return @lower if probability.zero?
      return @upper if probability == 1.0

      z =
        if @lower.finite? && @upper.finite?
          bounded_quantile(probability)
        elsif @lower.finite?
          lower_quantile(probability)
        elsif @upper.finite?
          upper_quantile(probability)
        else
          normal_quantile(probability)
        end
      @mu + z * @sigma
    end

    alias_method :p_value, :quantile

    def sample(probability = random_probability)
      quantile(probability)
    end

    alias_method :rand, :sample

    def call(*arguments)
      sample(*arguments)
    end

    def reset_attempts
      @attempts = 0
    end

    private

    def initialize_lower
      @lower_z = (@lower - @mu) / @sigma
      @log_z = log_sf(@lower_z)
      @mean = @mu + @sigma * hazard_rate(@lower_z)
      @variance = @sigma * @sigma * lower_tail_variance(@lower_z)
    end

    def initialize_upper
      @upper_z = (@upper - @mu) / @sigma
      @log_z = log_cdf(@upper_z)
      @mean = @mu - @sigma * hazard_rate(@upper_z.abs)
      @variance = @sigma * @sigma * lower_tail_variance(-@upper_z)
    end

    def initialize_bounded
      a = (@lower - @mu) / @sigma
      b = (@upper - @mu) / @sigma
      @lower_z = a
      @upper_z = b
      @log_z = bounded_log_mass(a, b)

      if width_is_narrow?(a, b)
        @mean, @variance = narrow_bounded_moments(a, b)
      elsif a >= 0.0
        @mean, @variance = positive_bounded_moments(a, b)
      elsif b <= 0.0
        mean_z, variance_z = positive_bounded_moments(-b, -a)
        @mean = -mean_z
        @variance = variance_z
      else
        @mean, @variance = crossing_bounded_moments(a, b)
      end
      @mean = @mu + @sigma * @mean
      @variance = @sigma * @sigma * @variance
    end

    def width_is_narrow?(a, b)
      width = b - a
      width <= 1.0e-4 && width * (1.0 + a.abs) <= 10.0
    end

    def narrow_bounded_moments(a, b)
      width = b - a
      midpoint = (a + b) / 2.0
      linear = -midpoint * width
      quadratic = -(width * width) / 4.0
      coefficients = [1.0, linear]

      2.upto(120) do |n|
        coefficients << (
          linear * coefficients[n - 1] +
          quadratic * 2.0 * (n - 1) * coefficients[n - 2]
        ) / n.to_f
      end

      integrals = Array.new(3, 0.0)
      3.times do |moment|
        previous = nil
        coefficients.each_with_index do |coefficient, n|
          term = coefficient / (n + moment + 1).to_f
          break if previous && term.abs >= previous
          integrals[moment] += term
          previous = term.abs
        end
      end

      normalized_mean = integrals[1] / integrals[0]
      mean_z = a + width * normalized_mean
      normalized_raw = integrals[2] / integrals[0]
      variance_z = width * width *
                   (normalized_raw - normalized_mean * normalized_mean)
      [mean_z, variance_z]
    end

    def positive_bounded_moments(a, b)
      tail_a = mills_ratio(a)
      tail_b = mills_ratio(b)
      r = bounded_mass_ratio(a, b, tail_a, tail_b)
      hazard_a = hazard_rate(a, tail_a)
      hazard_b = hazard_rate(b, tail_b)
      one_minus = 1.0 - r
      mean_z = (hazard_a - r * hazard_b) / one_minus
      raw_z = (1.0 + a * hazard_a - r * (1.0 + b * hazard_b)) / one_minus
      [mean_z, raw_z - mean_z * mean_z]
    end

    def crossing_bounded_moments(a, b)
      left_mass = normal_cdf(a)
      right_mass = normal_sf(b)
      total = Math.exp(@log_z)
      left_weight = left_mass / total
      right_weight = right_mass / total
      left_hazard = hazard_rate(-a)
      right_hazard = hazard_rate(b)
      mean_z = right_weight * right_hazard - left_weight * left_hazard
      raw_z = (normal_cdf(b) - b * normal_pdf(b) -
               normal_cdf(a) + a * normal_pdf(a)) / total
      [mean_z, raw_z - mean_z * mean_z]
    end

    def bounded_log_mass(a, b)
      if a >= 0.0
        log_sf(a) + log1p(-bounded_mass_ratio(a, b))
      elsif b <= 0.0
        log_cdf(b) + log1p(-bounded_mass_ratio(-b, -a))
      else
        Math.log(normal_cdf(b) - normal_cdf(a))
      end
    end

    def bounded_mass_ratio(a, b, tail_a = mills_ratio(a),
                           tail_b = mills_ratio(b))
      Math.exp(log_phi(b) - log_phi(a) + Math.log(tail_b) -
               Math.log(tail_a))
    end

    def bounded_cdf(z)
      a = @lower_z
      b = @upper_z
      if a >= 0.0
        r = bounded_mass_ratio(a, b)
        Math.exp(log1p(-bounded_mass_ratio(a, z)) -
                 log1p(-r))
      elsif b <= 0.0
        reflected_a = -b
        reflected_b = -a
        reflected_x = -z
        r = bounded_mass_ratio(reflected_a, reflected_b)
        Math.exp(log1p(-bounded_mass_ratio(reflected_a, reflected_x)) -
                 log1p(-r))
      else
        (normal_cdf(z) - normal_cdf(a)) / Math.exp(@log_z)
      end
    end

    def bounded_quantile(probability)
      a = @lower_z
      b = @upper_z
      if a >= 0.0
        r = bounded_mass_ratio(a, b)
        target = log_sf(a) + Math.log(1.0 - probability * (1.0 - r))
        inverse_log_sf(target)
      elsif b <= 0.0
        r = bounded_mass_ratio(-b, -a)
        target = log_sf(-b) + Math.log(1.0 - probability * (1.0 - r))
        -inverse_log_sf(target)
      elsif probability <= 0.5
        normal_quantile(normal_cdf(a) + probability * Math.exp(@log_z))
      else
        target = normal_sf(b) + (1.0 - probability) * Math.exp(@log_z)
        if target < 0.005
          inverse_log_sf(Math.log(target))
        else
          normal_quantile(1.0 - target)
        end
      end
    end

    def lower_quantile(probability)
      a = @lower_z
      inverse_log_sf(log_sf(a) + Math.log(1.0 - probability))
    end

    def upper_quantile(probability)
      b = @upper_z
      inverse_log_cdf(log_cdf(b) + Math.log(probability))
    end

    def check_probability!(probability)
      unless probability.is_a?(Numeric) && probability >= 0.0 && probability <= 1.0
        raise ArgumentError, "probability must be between 0 and 1"
      end
    end

    def random_probability
      @attempts += 1
      if @random.respond_to?(:call) && !@random.respond_to?(:rand)
        @random.call
      else
        @random.rand
      end
    end

    def normal_pdf(z)
      Math.exp(log_phi(z))
    end

    def normal_cdf(z)
      if z < 0.0
        0.5 * Math.erfc(-z / SQRT_2)
      else
        1.0 - 0.5 * Math.erfc(z / SQRT_2)
      end
    end

    def normal_sf(z)
      normal_cdf(-z)
    end

    def log_phi(z)
      -0.5 * z * z - 0.5 * LOG_2PI
    end

    def log_cdf(z)
      return log_phi(z) - Math.log(hazard_rate(z)) if z >= 0.0
      return log_phi(z) + Math.log(mills_ratio(-z)) if z < -26.0

      Math.log(0.5) + Math.log(Math.erfc(-z / SQRT_2))
    end

    def log_sf(z)
      if z <= 26.0
        Math.log(0.5) + Math.log(Math.erfc(z / SQRT_2))
      else
        log_phi(z) + Math.log(mills_ratio(z))
      end
    end

    def mills_ratio(z)
      z = z.to_f
      return normal_sf(z) / normal_pdf(z) if z < 8.0

      sum = 1.0 / z
      term = sum
      k = 1
      loop do
        next_term = term * (2.0 * k - 1.0) / (z * z)
        break if next_term.abs >= term.abs
        k.odd? ? sum -= next_term : sum += next_term
        term = next_term
        k += 1
      end
      sum
    end

    def hazard_rate(z, tail = mills_ratio(z))
      return 1.0 / tail if z < 10.0

      sum = 0.0
      term = 1.0
      previous = nil
      HAZARD_COEFFICIENTS.each do |coefficient|
        value = coefficient * term
        break if previous && value.abs >= previous
        sum += value
        previous = value.abs unless value.zero?
        term /= z * z
      end
      z * sum
    end

    def lower_tail_variance(z)
      if z < 10.0
        hazard = hazard_rate(z)
        return 1.0 + z * hazard - hazard * hazard
      end

      variance = 0.0
      previous = nil
      TAIL_VARIANCE_COEFFICIENTS.each_with_index do |coefficient, n|
        term = coefficient / z**(2 * n)
        break if n > 1 && previous && term.abs >= previous
        variance += term
        previous = term.abs unless coefficient.zero?
      end
      variance
    end

    def log1p(value)
      if value.abs < 1.0e-4
        value - value * value / 2.0 + value**3 / 3.0
      else
        Math.log(1.0 + value)
      end
    end

    def inverse_log_sf(target)
      x = if target < -2.0
            Math.sqrt(-2.0 * target)
          else
            1.0
          end
      8.times do
        tail = mills_ratio(x)
        x -= (log_phi(x) + Math.log(tail) - target) * -tail
      end
      x
    end

    def inverse_log_cdf(target)
      -inverse_log_sf(target)
    end

    def normal_quantile(probability)
      return -Float::INFINITY if probability.zero?
      return Float::INFINITY if probability == 1.0

      tail_probability = [probability, 1.0 - probability].min
      t = Math.sqrt(-2.0 * Math.log(tail_probability))
      x = t - (2.30753 + 0.27061 * t) /
          (1.0 + 0.99229 * t + 0.04481 * t * t)
      x = -x if probability < 0.5

      2.times do
        error = normal_cdf(x) - probability
        x -= error / normal_pdf(x)
      end
      x
    end
  end
end
