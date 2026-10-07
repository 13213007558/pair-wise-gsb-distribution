module Distribution
  class TruncatedNormal
    attr_reader :lower, :upper, :normal_mean, :normal_stddev, :rng, :attempts

    def initialize(lower = -Float::INFINITY, upper = Float::INFINITY,
                   mean: 0.0, stddev: 1.0, rng: nil, seed: nil)
      @lower = lower.to_f
      @upper = upper.to_f
      @normal_mean = mean.to_f
      @normal_stddev = stddev.to_f

      unless @normal_stddev.finite? && @normal_stddev > 0
        raise ArgumentError, 'stddev must be finite and positive'
      end
      unless @normal_mean.finite?
        raise ArgumentError, 'mean must be finite'
      end
      unless valid_bound?(@lower, true) && valid_bound?(@upper, false)
        raise ArgumentError, 'invalid truncation bounds'
      end
      unless @lower < @upper
        raise ArgumentError, 'lower bound must be less than upper bound'
      end

      @rng = setup_rng(rng, seed)
      @attempts = 0

      @a = standardize(@lower)
      @b = standardize(@upper)
      @log_mass = log_interval_mass(@a, @b)
      unless @log_mass.finite? && @log_mass > -744.0
        raise RangeError, 'truncation interval has zero probability in Float'
      end

      if finite?(@a) && finite?(@b)
        @standard_mean, @standard_variance = finite_standard_moments(@a, @b)
      elsif finite?(@a)
        hazard = Math.exp(log_standard_pdf(@a) - log_standard_survival(@a))
        @standard_mean = hazard
        @standard_variance = 1 + @a * hazard - hazard**2
      elsif finite?(@b)
        inverse_mills = Math.exp(log_standard_pdf(@b) - log_standard_cdf(@b))
        @standard_mean = -inverse_mills
        @standard_variance = 1 - @b * @standard_mean - @standard_mean**2
      else
        @standard_mean = 0.0
        @standard_variance = 1.0
      end

      @mean = @normal_mean + @standard_mean * @normal_stddev
      @stddev = @normal_stddev * Math.sqrt(@standard_variance)
    end

    attr_reader :mean, :stddev

    def variance
      @standard_variance * @normal_stddev**2
    end

    def pdf(x)
      Math.exp(log_pdf(x))
    end

    def log_pdf(x)
      value = x.to_f
      return -Float::INFINITY if value < @lower || value > @upper

      z = standardize(value)
      log_standard_pdf(z) - @log_mass - Math.log(@normal_stddev)
    end

    def cdf(x)
      value = x.to_f
      return 0.0 if value <= @lower
      return 1.0 if value >= @upper

      z = standardize(value)
      if finite?(@a) && finite?(@b)
        cdf_with_finite_bounds(z)
      elsif finite?(@a)
        -expm1(log_standard_survival(z) - log_standard_survival(@a))
      elsif finite?(@b)
        if z >= 0.0
          1.0 - Math.exp(log_standard_survival(z) - log_standard_cdf(@b))
        else
          Math.exp(log_standard_cdf(z) - log_standard_cdf(@b))
        end
      else
        Math.exp(log_standard_cdf(z))
      end
    end

    def quantile(probability)
      p = probability.to_f
      raise RangeError, 'probability must be in [0, 1]' unless p >= 0.0 && p <= 1.0
      return @lower if p.zero?
      return @upper if p == 1.0

      z =
        if finite?(@a) && finite?(@b)
          quantile_with_finite_bounds(p)
        elsif finite?(@a)
          inverse_standard_survival_log(Math.log(1.0 - p) + log_standard_survival(@a))
        elsif finite?(@b)
          inverse_standard_cdf_log(Math.log(p) + log_standard_cdf(@b))
        else
          inverse_standard_cdf_log(Math.log(p))
        end

      z = @a if finite?(@a) && z < @a
      z = @b if finite?(@b) && z > @b
      @normal_mean + z * @normal_stddev
    end

    def sample
      @attempts += 1
      probability = @rng.call
      probability = 0.0 if probability < 0.0
      probability = 1.0 if probability > 1.0
      quantile(probability)
    end
    alias rand sample

    def reset_attempts
      @attempts = 0
    end

    private

    def setup_rng(custom_rng, seed)
      if custom_rng.respond_to?(:call)
        custom_rng
      elsif custom_rng.respond_to?(:rand)
        -> { custom_rng.rand }
      elsif !custom_rng.nil?
        raise ArgumentError, 'rng must respond to call or rand'
      else
        Random.new(seed || Random.new_seed).method(:rand).to_proc
      end
    end

    def valid_bound?(value, is_lower)
      value.finite? || (is_lower && value == -Float::INFINITY) ||
        (!is_lower && value == Float::INFINITY)
    end

    def finite?(value)
      value.finite?
    end

    def standardize(value)
      (value - @normal_mean) / @normal_stddev
    end

    def finite_standard_moments(a, b)
      width = b - a
      center = (a + b) / 2.0

      if width * (center.abs + width) < 0.1
        truncated_standard_moments(a, b)
      elsif a >= 0.0
        positive_standard_moments(a, b)
      elsif b <= 0.0
        reflected_mean, reflected_variance = positive_standard_moments(-b, -a)
        [-reflected_mean, reflected_variance]
      else
        density_a = Math.exp(log_standard_pdf(a))
        density_b = Math.exp(log_standard_pdf(b))
        mass = Math.exp(@log_mass)
        standard_mean = (density_a - density_b) / mass
        second_moment = 1 + (a * density_a - b * density_b) / mass
        [standard_mean, second_moment - standard_mean**2]
      end
    end

    def positive_standard_moments(a, b)
      width = b - a
      density_ratio_a = Math.exp(log_standard_pdf(a) - @log_mass)
      density_ratio_b = Math.exp(log_standard_pdf(b) - @log_mass)
      standard_mean = density_ratio_a - density_ratio_b
      offset_mean = standard_mean - a
      variance = 1 - a * offset_mean - offset_mean**2 - width * density_ratio_b
      [standard_mean, variance]
    end

    def split_probability(a, b)
      return 0.0 unless a <= 0.0 && b >= 0.0

      Math.exp(logaddexp(log_standard_cdf(a), Math.log(0.5) + @log_mass))
    end

    def quantile_with_finite_bounds(p)
      if @a >= 0.0
        inverse_standard_survival_log(
          logaddexp(log_standard_survival(@b), Math.log(1.0 - p) + @log_mass)
        )
      elsif @b <= 0.0
        inverse_standard_cdf_log(
          logaddexp(log_standard_cdf(@a), Math.log(p) + @log_mass)
        )
      elsif p <= split_probability(@a, @b)
        inverse_standard_cdf_log(
          logaddexp(log_standard_cdf(@a), Math.log(p) + @log_mass)
        )
      else
        inverse_standard_survival_log(
          logaddexp(log_standard_survival(@b), Math.log(1.0 - p) + @log_mass)
        )
      end
    end

    def log_interval_mass(a, b)
      if !finite?(a)
        log_standard_cdf(b)
      elsif !finite?(b)
        log_standard_survival(a)
      elsif a >= 0.0
        log_sub_exp(log_standard_survival(a), log_standard_survival(b))
      elsif b <= 0.0
        log_sub_exp(log_standard_cdf(b), log_standard_cdf(a))
      else
        Math.log(0.5 * (Math.erf(b / Math.sqrt(2)) - Math.erf(a / Math.sqrt(2))))
      end
    end

    def cdf_with_finite_bounds(z)
      if @a >= 0.0
        positive_finite_cdf(z)
      elsif z <= 0.0
        negative_finite_cdf(z)
      else
        crossing_finite_cdf(z)
      end
    end

    def negative_finite_cdf(z)
      if z - @a <= @b - z
        local = -expm1(log_standard_cdf(@a) - log_standard_cdf(z))
        local * Math.exp(log_standard_cdf(z) - @log_mass)
      else
        local = -expm1(log_standard_cdf(z) - log_standard_cdf(@b))
        1.0 - local * Math.exp(log_standard_cdf(@b) - @log_mass)
      end
    end

    def positive_finite_cdf(z)
      if z - @a <= @b - z
        local = -expm1(log_standard_survival(z) - log_standard_survival(@a))
        local * Math.exp(log_standard_survival(@a) - @log_mass)
      else
        right_mass = log_sub_exp(log_standard_survival(z), log_standard_survival(@b))
        1.0 - Math.exp(right_mass - @log_mass)
      end
    end

    def crossing_finite_cdf(z)
      left_log = Math.log(0.5 * (Math.erf(z / Math.sqrt(2)) - Math.erf(@a / Math.sqrt(2))))
      right_log = Math.log(0.5 * (Math.erf(@b / Math.sqrt(2)) - Math.erf(z / Math.sqrt(2))))
      return 1.0 - Math.exp(right_log - @log_mass) if left_log > right_log

      Math.exp(left_log - @log_mass)
    end

    def truncated_standard_moments(a, b)
      center = (a + b) / 2.0
      half_width = (b - a) / 2.0
      reference = a >= 0.0 ? a : (b <= 0.0 ? b : 0.0)
      reference_log_pdf = log_standard_pdf(reference)

      total = 0.0
      first_central_moment = 0.0
      second_central_moment = 0.0

      GAUSS_LEGENDRE_32.each do |node, weight|
        z = center + half_width * node
        relative_weight = half_width * weight *
                          Math.exp(log_standard_pdf(z) - reference_log_pdf)
        offset = z - center
        total += relative_weight
        first_central_moment += offset * relative_weight
        second_central_moment += offset**2 * relative_weight
      end

      centered_mean = first_central_moment / total
      [
        center + centered_mean,
        second_central_moment / total - centered_mean**2
      ]
    end

    def log_standard_cdf(z)
      return Math.log(0.5 * Math.erfc(-z / Math.sqrt(2))) if z < 6.0

      Math.log(-expm1(log_standard_survival(z)))
    end

    def log_standard_survival(z)
      return Math.log(0.5 * Math.erfc(z / Math.sqrt(2))) if z > -6.0

      Math.log(-expm1(log_standard_cdf(z)))
    end

    def log_standard_pdf(z)
      -(Math.log(2 * Math::PI) / 2.0 + z * z / 2.0)
    end

    def log_sub_exp(first, second)
      return first if first == second
      return second if second > first

      first + Math.log(1.0 - Math.exp(second - first))
    end

    def expm1(value)
      value.abs < 1e-4 ? value + value**2 / 2.0 + value**3 / 6.0 : Math.exp(value) - 1.0
    end

    def logaddexp(first, second)
      larger = [first, second].max
      larger + Math.log(Math.exp(first - larger) + Math.exp(second - larger))
    end

    def inverse_standard_cdf_log(log_probability)
      return -Float::INFINITY if log_probability == -Float::INFINITY
      return Float::INFINITY if log_probability >= 0.0

      if log_probability > -700.0
        probability = Math.exp(log_probability)
        x = inverse_standard_cdf(probability)
      end
      refine_inverse_standard_cdf_log(x, log_probability, :log_standard_cdf, 1.0)
    end

    def inverse_standard_survival_log(log_probability)
      return Float::INFINITY if log_probability == -Float::INFINITY
      return -Float::INFINITY if log_probability >= 0.0

      if log_probability > -700.0
        probability = Math.exp(log_probability)
        x = -inverse_standard_cdf(probability)
      else
        x = asymptotic_inverse_log(log_probability)
      end
      refine_inverse_standard_cdf_log(x, log_probability, :log_standard_survival, -1.0)
    end

    def inverse_standard_cdf(probability)
      p = probability
      return -Float::INFINITY if p.zero?
      return Float::INFINITY if p == 1.0

      lower_a = [-3.969683028665376e+01, 2.209460984245205e+02,
                 -2.759285104469687e+02, 1.383577518672690e+02,
                 -3.066479806614716e+01, 2.506628277459239e+00]
      lower_b = [-5.447609879822406e+01, 1.615858368580409e+02,
                 -1.556989798598866e+02, 6.680131188771972e+01,
                 -1.328068155288572e+01]
      upper_c = [-7.784894002430293e-03, -3.223964580411365e-01,
                 -2.400758277161838e+00, -2.549732539343734e+00,
                 4.374664141464968e+00, 2.938163982698783e+00]
      upper_d = [7.784695709041462e-03, 3.224671290700398e-01,
                 2.445134137142996e+00, 3.754408661907416e+00]

      p_low = 0.02425
      p_high = 1.0 - p_low

      if p < p_low
        q = Math.sqrt(-2.0 * Math.log(p))
        x = polynomial(lower_a, q) / (q * polynomial(lower_b, q) + 1.0)
      elsif p <= p_high
        q = p - 0.5
        r = q * q
        x = q * polynomial(upper_c, r) / (polynomial(upper_d, r) + 1.0)
      else
        q = Math.sqrt(-2.0 * Math.log(1.0 - p))
        x = -polynomial(lower_a, q) / (q * polynomial(lower_b, q) + 1.0)
      end

      refine_inverse_standard_cdf(x, p)
    end

    def refine_inverse_standard_cdf(x, target)
      4.times do
        error = Math.exp(log_standard_cdf(x)) - target
        delta = error / Math.exp(log_standard_pdf(x))
        next_x = x - delta
        break x = next_x if (next_x - x).abs <= Float::EPSILON * x.abs

        x = next_x
      end
      [[x, -37.0].max, 37.0].min
    end

    def asymptotic_inverse_log(log_probability)
      y = -2.0 * log_probability
      root_y = Math.sqrt(y)
      log_y = Math.log(y)
      x = root_y - (2.30753 + 0.27061 / root_y) /
                  (1.0 + 0.99229 / root_y + 0.04481 * log_y / y)
      [[x, -37.0].max, 37.0].min
    end

    def refine_inverse_standard_cdf_log(x, target_log, log_function, derivative_sign)
      6.times do
        current_log = send(log_function, x)
        density_log = log_standard_pdf(x)
        break if current_log.nan?

        next_x = x - derivative_sign * (current_log - target_log) *
                 Math.exp(current_log - density_log)
        break x if next_x.nan?

        next_x = [[next_x, -37.0].max, 37.0].min
        break x = next_x if (next_x - x).abs <= Float::EPSILON * x.abs

        x = next_x
      end
      [[x, -37.0].max, 37.0].min
    end

    def polynomial(coefficients, variable)
      coefficients.inject { |result, coefficient| result * variable + coefficient }
    end

    def self.build_gauss_legendre(order)
      nodes = []
      weights = []

      (0...((order + 1) / 2)).each do |index|
        x = Math.cos(Math::PI * (index + 0.75) / (order + 0.5))

        loop do
          p_minus_1 = 0.0
          p = 1.0

          (1..order).each do |degree|
            p_minus_2 = p_minus_1
            p_minus_1 = p
            p = ((2.0 * degree - 1.0) * x * p_minus_1 -
                 (degree - 1.0) * p_minus_2) / degree
          end

          derivative = order * (x * p - p_minus_1) / (x**2 - 1.0)
          previous = x
          x -= p / derivative
          break if (x - previous).abs < 1e-15
        end

        p_minus_1 = 0.0
        p = 1.0
        (1..order).each do |degree|
          p_minus_2 = p_minus_1
          p_minus_1 = p
          p = ((2.0 * degree - 1.0) * x * p_minus_1 -
               (degree - 1.0) * p_minus_2) / degree
        end
        derivative = order * (x * p - p_minus_1) / (x**2 - 1.0)
        weight = 2.0 / ((1.0 - x**2) * derivative**2)

        nodes.concat([-x, x])
        weights.concat([weight, weight])
      end

      nodes.zip(weights).sort_by { |node, _| node }
    end
    private_class_method :build_gauss_legendre

    GAUSS_LEGENDRE_32 = build_gauss_legendre(32)
  end
end
