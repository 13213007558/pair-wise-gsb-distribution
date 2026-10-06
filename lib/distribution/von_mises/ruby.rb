module Distribution
  module VonMises
    module Ruby_
      TWO_PI = 2.0 * Math::PI
      ASYMPTOTIC_KAPPA = 100.0
      SCALED_TAIL_LIMIT = 40.0
      QUANTILE_ITERATIONS = 64

      class << self
        # Probability density at +theta+ for a von Mises distribution.
        # +mu+ is the mean direction in radians and may be outside [0, 2*pi].
        def pdf(theta, mu = 0, kappa = 1)
          theta = finite_float(theta, 'theta')
          mu = finite_float(mu, 'mu')
          kappa = concentration(kappa)
          delta = wrap_pi(theta - mu)

          Math.exp(log_pdf(delta, kappa))
        end

        alias_method :density, :pdf

        # CDF on the one-turn interval [+cut+, +cut+2*pi].
        def cdf(x, mu = 0, kappa = 1, cut = 0)
          x = finite_float(x, 'x')
          mu = finite_float(mu, 'mu')
          kappa = concentration(kappa)
          cut = finite_float(cut, 'cut')
          upper_bound = cut + TWO_PI

          return 1.0 if x == upper_bound

          if x < cut || x > upper_bound
            raise RangeError,
                  "x must be in [cut, cut + 2*pi] (#{cut}..#{upper_bound}), got #{x}"
          end

          length = wrap_two_pi(x - cut)
          probability = forward_arc_probability(cut, length, mu, kappa)
          clamp(probability, 0.0, 1.0)
        end

        # Probability of the directed arc starting at +start_angle+ and having
        # signed +arc_length+. A negative length runs clockwise. Arcs with zero
        # length have probability zero; arcs covering a turn have probability one.
        def arc_probability(start_angle, arc_length, mu = 0, kappa = 1)
          start_angle = finite_float(start_angle, 'start_angle')
          arc_length = finite_float(arc_length, 'arc_length')
          mu = finite_float(mu, 'mu')
          kappa = concentration(kappa)

          return 0.0 if arc_length == 0.0
          return 1.0 if arc_length.abs >= TWO_PI

          if arc_length > 0.0
            forward_arc_probability(start_angle, arc_length, mu, kappa)
          else
            forward_arc_probability(start_angle + arc_length, -arc_length, mu, kappa)
          end
        end

        def quantile(probability, mu = 0, kappa = 1, cut = 0)
          probability = finite_float(probability, 'probability')
          mu = finite_float(mu, 'mu')
          kappa = concentration(kappa)
          cut = finite_float(cut, 'cut')

          unless probability >= 0.0 && probability <= 1.0
            raise RangeError, "probability must be in [0, 1], got #{probability}"
          end

          return cut.to_f if probability == 0.0
          return cut + TWO_PI if probability == 1.0
          return cut + TWO_PI * probability if kappa == 0.0

          lower = cut.to_f
          upper = cut + TWO_PI

          QUANTILE_ITERATIONS.times do
            middle = (lower + upper) / 2.0
            if cdf(middle, mu, kappa, cut) < probability
              lower = middle
            else
              upper = middle
            end
          end

          (lower + upper) / 2.0
        end

        alias_method :p_value, :quantile

        private

        def forward_arc_probability(start_angle, arc_length, mu, kappa)
          from = start_angle - mu
          to = from + arc_length
          cumulative_on_line(to, kappa) - cumulative_on_line(from, kappa)
        end

        def cumulative_on_line(value, kappa)
          cycles = (value / TWO_PI).floor
          remainder = value - cycles * TWO_PI
          base = cycles.to_f

          return base + half_circle_probability(remainder, kappa) if remainder <= Math::PI

          base + 1.0 - half_circle_probability(TWO_PI - remainder, kappa)
        end

        def half_circle_probability(distance, kappa)
          return 0.0 if distance == 0.0
          return 0.5 if distance == Math::PI
          return distance / TWO_PI if kappa == 0.0

          scale = Math.sqrt(kappa)
          upper = [distance * scale, SCALED_TAIL_LIMIT].min
          clamp(integrate_scaled(upper, kappa) / scale, 0.0, 0.5)
        end

        def integrate_scaled(upper, kappa)
          scale = Math.sqrt(kappa)

          density = lambda do |u|
            angle = u / scale
            Math.exp(cosine_exponent(angle, kappa) + log_mode_density(kappa))
          end

          left = 0.0
          right = upper
          middle = upper / 2.0
          fl = density.call(left)
          fm = density.call(middle)
          fr = density.call(right)
          whole = simpson(left, right, fl, fm, fr)
          intervals = [[left, right, fl, fm, fr, whole]]
          tolerance = 1e-12
          result = 0.0

          while (interval = intervals.pop)
            a, b, fa, fmid, fb, total = interval
            m = (a + b) / 2.0
            lm = (a + m) / 2.0
            rm = (m + b) / 2.0
            flm = density.call(lm)
            frm = density.call(rm)
            left_total = simpson(a, m, fa, flm, fmid)
            right_total = simpson(m, b, fmid, frm, fb)
            refined = left_total + right_total
            error = (refined - total).abs / 15.0

            if intervals.size > 64 ||
              lm == a || rm == b || error <= tolerance
              result += refined + (refined - total) / 15.0
              next
            end

            intervals.push([a, m, fa, flm, fmid, left_total],
                           [m, b, fmid, frm, fb, right_total])
          end

          result
        end

        def simpson(left, right, fl, fm, fr)
          (right - left) * (fl + 4.0 * fm + fr) / 6.0
        end

        def cosine_exponent(delta, kappa)
          -2.0 * kappa * Math.sin(delta / 2.0)**2
        end

        def log_pdf(delta, kappa)
          cosine_exponent(delta, kappa) + log_mode_density(kappa)
        end

        def log_mode_density(kappa)
          return -Math.log(TWO_PI) if kappa == 0.0

          if kappa <= ASYMPTOTIC_KAPPA
            kappa - Math.log(TWO_PI) - Math.log(i0_series(kappa))
          else
            0.5 * Math.log(kappa / TWO_PI) - Math.log(i0_asymptotic_factor(kappa))
          end
        end

        def i0_series(kappa)
          sum = 1.0
          term = 1.0
          n = 1

          loop do
            term *= (kappa**2).fdiv(4 * n * n)
            sum += term
            break if term <= Float::EPSILON * sum
            n += 1
          end

          sum
        end

        def i0_asymptotic_factor(kappa)
          term = 1.0
          sum = 1.0
          best_abs = 1.0

          1.upto(1000) do |m|
            term *= ((2 * m - 1)**2).to_f / (8.0 * m * kappa)
            break if term.abs >= best_abs

            best_abs = term.abs
            sum += term
            break if term.abs <= Float::EPSILON
          end

          sum
        end

        def wrap_pi(value)
          value = value % TWO_PI
          value -= TWO_PI if value > Math::PI
          value
        end

        def wrap_two_pi(value)
          value % TWO_PI
        end

        def concentration(value)
          kappa = finite_float(value, 'kappa')
          if kappa < 0.0
            raise RangeError, "kappa must be non-negative, got #{kappa}"
          end

          kappa
        end

        def finite_float(value, name)
          number = Float(value)
          unless number.finite?
            raise ArgumentError, "#{name} must be finite, got #{value.inspect}"
          end

          number
        end

        def clamp(value, lower, upper)
          return lower if value < lower
          return upper if value > upper
          value
        end
      end
    end
  end
end
