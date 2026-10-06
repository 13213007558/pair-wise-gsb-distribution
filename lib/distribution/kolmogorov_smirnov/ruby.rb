require 'bigdecimal'

module Distribution
  module KolmogorovSmirnov
    module Ruby_
      # Exact finite-n distribution of the one-sample two-sided
      # Kolmogorov-Smirnov statistic D_n under a continuous null.
      #
      # The exact cdf is computed with Pomeranz's recursion
      # (Pomeranz, 1974): if a_j = max(0, j/n - d) and
      # b_j = min(1, (j-1)/n + d), then P(D_n < d) equals the
      # probability that a_j < U_(j) < b_j for all order statistics
      # U_(j) of n i.i.d. Uniform(0,1) variables. Writing w_j for the
      # density of U_(j) restricted to the event and
      # M_j(t) = \int_{a_j}^t w_j(u) (1-u)^{-(n-j)} du, the Markov
      # property of uniform order statistics gives
      #   M_1(t)   = n (t - a_1),
      #   M_{j+1}(t) = (n-j) \int_{a_{j+1}}^t M_j^*(v) dv,
      # where M_j^* extends M_j constantly beyond b_j, and
      # P(D_n < d) = M_n(b_n). Every M_j is piecewise polynomial with
      # knots at the a_j and b_j, and all coefficients stay
      # non-negative, so the recursion involves no cancellation and
      # can be evaluated in Float, BigDecimal or exact Rational
      # arithmetic.
      class << self
        # Exact P(D_n <= d) for the two-sided one-sample statistic.
        def cdf(d, n)
          check_n(n)
          d = d.to_f
          return 0.0 if d <= 0.5 / n
          return 1.0 if d >= 1.0
          return small_d_cdf(d, n) if d <= 1.0 / n

          pomeranz_cdf_float(d, n)
        end

        # Exact P(D_n > d) for the two-sided one-sample statistic.
        #
        # This is NOT simply 1 - cdf everywhere: in the right tail
        # (cdf > 0.9) the cdf is recomputed with BigDecimal arithmetic
        # at a precision adapted to the tail magnitude, so very small
        # survival probabilities keep significant digits.
        def sf(d, n)
          check_n(n)
          d = d.to_f
          return 1.0 if d <= 0.5 / n
          return 0.0 if d >= 1.0
          return 1.0 - small_d_cdf(d, n) if d <= 1.0 / n

          c = pomeranz_cdf_float(d, n)
          return 1.0 - c if c <= 0.9

          upper_tail_sf(d, n)
        end

        # Exact quantile (inverse cdf) of D_n.
        #
        # Returns the smallest d with cdf(d) >= pr; the previous
        # representable Float satisfies cdf(d) < pr. Comparisons in
        # the extreme right tail (pr > 0.9) are made with the
        # accurately computed survival function, i.e. the returned d
        # satisfies sf(d) <= 1 - pr and sf(d.prev_float) > 1 - pr.
        def p_value(pr, n)
          check_n(n)
          raise ArgumentError, 'pr must lie in [0, 1]' if pr < 0 || pr > 1
          return 0.5 / n if pr == 0
          return 1.0 if pr == 1

          lo = 0.5 / n
          hi = 1.0
          if pr <= 0.9
            50.times do
              mid = lo + (hi - lo) / 2.0
              cdf(mid, n) < pr ? lo = mid : hi = mid
            end
            d = hi
            d = d.next_float while cdf(d, n) < pr
            d = d.prev_float while cdf(d.prev_float, n) >= pr
            d
          else
            q = 1.0 - pr
            50.times do
              mid = lo + (hi - lo) / 2.0
              sf(mid, n) > q ? lo = mid : hi = mid
            end
            d = hi
            d = d.next_float while sf(d, n) > q
            d = d.prev_float while sf(d.prev_float, n) <= q
            d
          end
        end

        alias_method :quantile, :p_value

        # APPROXIMATE P(D_n <= d): Kolmogorov's limiting (n -> infinity)
        # distribution evaluated at sqrt(n) * d. This is an asymptotic
        # approximation, not the exact finite-n cdf; use #cdf for the
        # exact value.
        def cdf_asymptotic(d, n)
          check_n(n)
          return 0.0 if d <= 0

          z2 = n * d * d
          if z2 < 1.18
            kolmogorov_cdf_series(Math.sqrt(z2))
          else
            1.0 - kolmogorov_sf_series(z2)
          end
        end

        # APPROXIMATE P(D_n > d): upper tail of Kolmogorov's limiting
        # (n -> infinity) distribution evaluated at sqrt(n) * d. This is
        # an asymptotic approximation, not the exact finite-n survival
        # function; use #sf for the exact value.
        def sf_asymptotic(d, n)
          check_n(n)
          return 1.0 if d <= 0

          z2 = n * d * d
          if z2 < 1.18
            1.0 - kolmogorov_cdf_series(Math.sqrt(z2))
          else
            kolmogorov_sf_series(z2)
          end
        end

        private

        def check_n(n)
          raise ArgumentError, 'n must be a positive integer' unless n.is_a?(Integer) && n >= 1
        end

        # Exact P(D_n <= d) = n! (2d - 1/n)^n, valid for
        # 1/(2n) <= d <= 1/n (the band intervals are then disjoint).
        def small_d_cdf(d, n)
          Math.exp(Math.lgamma(n + 1).first + n * Math.log(2.0 * d - 1.0 / n))
        end

        def pomeranz_cdf_float(d, n)
          a = Array.new(n) do |j|
            t = (j + 1).fdiv(n) - d
            t < 0.0 ? 0.0 : t
          end
          b = Array.new(n) do |j|
            t = j.fdiv(n) + d
            t > 1.0 ? 1.0 : t
          end
          pomeranz_core(a, b, n)
        end

        def pomeranz_cdf_bigdecimal(d, n, precision)
          previous_limit = BigDecimal.limit(precision)
          begin
            dd = BigDecimal(d.to_s)
            nn = BigDecimal(n)
            a = Array.new(n) do |j|
              t = BigDecimal(j + 1) / nn - dd
              t < 0 ? BigDecimal(0) : t
            end
            b = Array.new(n) do |j|
              t = BigDecimal(j) / nn + dd
              t > 1 ? BigDecimal(1) : t
            end
            pomeranz_core(a, b, n)
          ensure
            BigDecimal.limit(previous_limit)
          end
        end

        # Exact upper tail P(D_n > d) with adaptive precision: the
        # working precision is sized from Kolmogorov's leading
        # asymptotic term and increased until at least ~15 significant
        # decimal digits of the survival probability are resolved.
        def upper_tail_sf(d, n)
          estimate = 2.0 * Math.exp(-2.0 * n * d * d)
          estimate = 1e-300 if estimate <= 0.0
          needed = (-Math.log10(estimate)).ceil
          precision = [30, needed + 20].max
          10.times do
            sf = 1 - pomeranz_cdf_bigdecimal(d, n, precision)
            return sf.to_f if sf > BigDecimal(10)**(15 - precision)

            precision += 30
          end
          raise 'KolmogorovSmirnov.sf: unable to reach required precision'
        end

        # Pomeranz recursion core. a[j], b[j] (0-based) are the clipped
        # band limits for U_(j+1); works with any numeric type
        # supporting +, -, *, / and comparison (Float, BigDecimal,
        # Rational). Requires 1/(2n) < d < 1.
        def pomeranz_core(a, b, n)
          knots = (a + b).uniq.sort
          # M_1(t) = n (t - a_1) on [a_1, b_1]
          pieces = []
          knots.each_cons(2) do |left, right|
            next if left < a[0] || right > b[0]

            pieces << [left, right, [n * (left - a[0]), n]]
          end
          (1...n).each do |j|
            last_left, _last_right, last_coeffs = pieces.last
            mass_at_b = horner(last_coeffs, b[j - 1] - last_left)
            extended = pieces.dup
            extended << [b[j - 1], 1, [mass_at_b]] if b[j - 1] < 1
            factor = n - j
            new_pieces = []
            running = 0
            extended.each do |left, right, coeffs|
              break if left >= b[j]
              next if right <= a[j]

              if left < a[j]
                # Integration starts strictly inside this piece
                # (only possible when d <= 1/n): re-expand the
                # polynomial around a[j].
                shift = a[j] - left
                coeffs = rebase(coeffs, shift)
                left = a[j]
              end
              hi = right > b[j] ? b[j] : right
              new_coeffs = [factor * running]
              coeffs.each_with_index do |c, k|
                new_coeffs << factor * c / (k + 1)
              end
              new_pieces << [left, hi, new_coeffs]
              width = right - left
              integral = 0
              coeffs.each_with_index do |c, k|
                integral += c * width**(k + 1) / (k + 1)
              end
              running += integral
            end
            pieces = new_pieces
          end
          last_left, _last_right, last_coeffs = pieces.last
          horner(last_coeffs, b[n - 1] - last_left)
        end

        # Coefficients of p(x + s) given the coefficients of p.
        def rebase(coeffs, s)
          degree = coeffs.size - 1
          out = Array.new(degree + 1, 0)
          (0..degree).each do |k|
            ck = coeffs[k]
            next if ck == 0

            binomial = 1
            (0..k).each do |q|
              out[q] += ck * binomial * s**(k - q)
              binomial = binomial * (k - q) / (q + 1)
            end
          end
          out
        end

        def horner(coeffs, x)
          coeffs.reverse.inject(0) { |sum, c| sum * x + c }
        end

        # Kolmogorov's limiting upper tail 2 sum_j (-1)^(j-1)
        # exp(-2 j^2 z2); converges fast for z2 >= 1.18.
        def kolmogorov_sf_series(z2)
          sum = 0.0
          j = 1
          loop do
            term = 2.0 * Math.exp(-2.0 * j * j * z2)
            sum = j.odd? ? sum + term : sum - term
            break if term < 1e-17

            j += 1
          end
          sum
        end

        # Kolmogorov's limiting cdf in the theta-function form
        # sqrt(2 pi)/z sum_j exp(-(2j+1)^2 pi^2 / (8 z^2));
        # converges fast for z < sqrt(1.18).
        def kolmogorov_cdf_series(z)
          sum = 0.0
          j = 0
          loop do
            term = Math.exp(-((2 * j + 1)**2) * Math::PI**2 / (8.0 * z * z))
            sum += term
            break if term < 1e-17

            j += 1
          end
          sum * Math.sqrt(2.0 * Math::PI) / z
        end
      end
    end
  end
end
