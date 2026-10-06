module Distribution
  module Vonmises
    module Ruby_
      TWO_PI = 2.0 * Math::PI

      # Module contain pure ruby implementation of the
      # von Mises distribution (rng, pdf, cdf, quantile and
      # circular arc probability functions).
      #
      # All angles are expressed in radians. The mean direction +mu+
      # is not required to lie inside [0, 2pi): any real value is
      # accepted, because angles are always reduced modulo one period
      # before evaluation. The concentration +kappa+ may range from
      # 0 (circular uniform) to arbitrarily large values (narrow peak);
      # all normalizing constants are computed in scaled form so no
      # Bessel function evaluation can overflow.
      class << self
        # Returns a lambda that emits a sequence of random angles
        # from a von Mises distribution, using the rejection
        # algorithm of Best & Fisher (1979).
        #
        # == Arguments
        #   * +mu+    - Mean direction, in radians (any real value)
        #   * +kappa+ - Concentration, must be >= 0
        #   * +seed+  - Seed to set the initial state, randomized if ommited
        #
        def rng(mu = 0.0, kappa = 1.0, seed = nil)
          seed = Random.new_seed if seed.nil?
          prng = Random.new(seed)
          mu, kappa = normalize_parameters(mu, kappa)
          return -> { mu + TWO_PI * (prng.rand - 0.5) } if kappa == 0.0

          a = 1.0 + Math.sqrt(1.0 + 4.0 * kappa * kappa)
          b = (a - Math.sqrt(2.0 * a)) / (2.0 * kappa)
          r = (1.0 + b * b) / (2.0 * b)
          lambda do
            loop do
              z = Math.cos(Math::PI * prng.rand)
              f = (1.0 + r * z) / (r + z)
              c = kappa * (r - f)
              u = prng.rand
              if (c * (2.0 - c) > u) || (c > 0.0 && Math.log(c / u) + 1.0 - c >= 0.0)
                return mu + (prng.rand < 0.5 ? -1.0 : 1.0) * Math.acos(f)
              end
            end
          end
        end

        # Von Mises probability density function.
        #
        # == Arguments
        #   * +x+     - Angle at which the density is evaluated, in radians
        #   * +mu+    - Mean direction, in radians (any real value)
        #   * +kappa+ - Concentration, must be >= 0
        #
        # The density is computed as
        #   exp(kappa * (cos(x - mu) - 1)) / (2 * pi * I0(kappa) * exp(-kappa))
        # so the Bessel normalization constant I0(kappa) * exp(-kappa)
        # stays below 1 and never overflows, however large kappa is.
        #
        def pdf(x, mu = 0.0, kappa = 1.0)
          mu, kappa = normalize_parameters(mu, kappa)
          return 1.0 / TWO_PI if kappa == 0.0

          t = reduce_angle(x - mu)
          Math.exp(kappa * (Math.cos(t) - 1.0)) / (TWO_PI * bessel_i0_scaled(kappa))
        end

        # Von Mises cumulative distribution function anchored at +cut+.
        #
        # == Arguments
        #   * +x+     - Angle at which the cdf is evaluated, in radians
        #   * +mu+    - Mean direction, in radians (any real value)
        #   * +kappa+ - Concentration, must be >= 0
        #   * +cut+   - Cut point (zero of the cdf), in radians.
        #               Defaults to +mu - PI+, opposite the mean direction.
        #
       # F(x) is the probability of the counterclockwise arc (cut, x].
        # The evaluation angle is taken modulo the caller-chosen
        # period, so F is strictly increasing from 0 to 1 on
        # [cut, cut + 2*PI) and angles that differ by whole periods
        # share the same value. Probability mass located across the
        # cut is fully accounted for: F(x) -> 1 as x approaches
        # cut + 2*PI from below.
        #
        def cdf(x, mu = 0.0, kappa = 1.0, cut = nil)
          mu, kappa = normalize_parameters(mu, kappa)
          cut = mu - Math::PI if cut.nil?
          d = (x - cut) % TWO_PI
          arc_integral(d, reduce_angle(cut - mu), series_coefficients(kappa))
        end

        # Von Mises inverse cdf / P-value function.
        #
        # == Arguments
        #   * +qn+    - Integral value, must be in [0, 1]
        #   * +mu+    - Mean direction, in radians (any real value)
        #   * +kappa+ - Concentration, must be >= 0
        #   * +cut+   - Cut point, in radians. Defaults to +mu - PI+.
        #
        # Returns the unique angle in [cut, cut + 2*PI] whose cdf equals
        # +qn+, so the quantile is monotone over the caller-specified
        # period.
        #
        def quantile(qn, mu = 0.0, kappa = 1.0, cut = nil)
          fail RangeError, 'cdf value(qn) must be from [0, 1]. '\
          "Cannot find quantile for qn=#{qn}" if qn > 1 || qn < 0

          mu, kappa = normalize_parameters(mu, kappa)
          cut = mu - Math::PI if cut.nil?
          return cut if qn == 0
          return cut + TWO_PI if qn == 1
          return cut + TWO_PI * qn if kappa == 0.0

          rhos = series_coefficients(kappa)
          t0 = reduce_angle(cut - mu)
          lo = 0.0
          hi = TWO_PI
          55.times do
            mid = 0.5 * (lo + hi)
            if arc_integral(mid, t0, rhos) < qn
              lo = mid
            else
              hi = mid
            end
          end
          cut + 0.5 * (lo + hi)
        end

        alias p_value quantile

        # Probability of a directed circular arc, given by its starting
        # angle and its signed length (both in radians).
        #
        # == Arguments
        #   * +start+  - Angle where the arc starts, in radians
        #   * +length+ - Signed arc length: positive sweeps
        #                counterclockwise (increasing angle), negative
        #                sweeps clockwise
        #   * +mu+     - Mean direction, in radians (any real value)
        #   * +kappa+  - Concentration, must be >= 0
        #
        # Arcs that cross 2*PI are followed through the angle wrap, so
        # a short arc from 359 deg to 1 deg is simply start=359 deg with
        # length=2 deg. A length of 0 is an empty arc (probability 0),
        # while a length of +/-2*PI or more is the full circle
        # (probability 1): the two are never confused.
        #
        def arc_probability(start, length, mu = 0.0, kappa = 1.0)
          mu, kappa = normalize_parameters(mu, kappa)
          return 0.0 if length == 0.0
          return 1.0 if length.abs >= TWO_PI

          if length < 0.0
            start += length
            length = -length
          end
          arc_integral(length, reduce_angle(start - mu), series_coefficients(kappa))
        end

        private

        # A negative concentration is equivalent to shifting the mean
        # direction by half a period.
        def normalize_parameters(mu, kappa)
          kappa < 0 ? [mu + Math::PI, -kappa] : [mu, kappa]
        end

        # Reduces an angle to [-PI, PI] so that trigonometric
        # evaluations keep full precision even for angles that are
        # many periods away from the origin.
        def reduce_angle(t)
          t %= TWO_PI
          t > Math::PI ? t - TWO_PI : t
        end

        # Scaled modified Bessel function I0(x) * exp(-x), which lies
        # in (0, 1] for x >= 0 and therefore cannot overflow.
        #
        # For x < 30 the absolutely convergent Taylor series of I0 is
        # summed and scaled afterwards; for larger x the asymptotic
        # expansion of the scaled function is used directly.
        def bessel_i0_scaled(x)
          if x < 30.0
            y = 0.25 * x * x
            sum = 1.0
            term = 1.0
            k = 0
            loop do
              k += 1
              term *= y / (k * k)
              sum += term
              break if term <= sum * 1e-17
            end
            sum * Math.exp(-x)
          else
            sum = 1.0
            term = 1.0
            k = 0
            loop do
              k += 1
              term *= ((2 * k - 1) * (2 * k - 1)) / (8.0 * k * x)
              break if term <= 1e-17
              sum += term
            end
            sum / Math.sqrt(TWO_PI * x)
          end
        end

        # Ratios rho_j = I_j(kappa) / I_0(kappa) for j = 1, 2, ...
        # needed by the Fourier series of the cdf.
        #
        # The consecutive ratios r_j = I_{j+1} / I_j satisfy the
        # backward recurrence r_{j-1} = 1 / (2j/kappa + r_j), evaluated
        # here as a continued fraction starting from r_N = 0 deep in
        # the tail (Miller's algorithm). Every r_j lies in (0, 1), so
        # the computation is stable and overflow-free for any kappa.
        # The list is truncated once rho_j is negligible.
        def series_coefficients(kappa)
          return [] if kappa == 0.0

          root = Math.sqrt(kappa)
          depth = ([kappa, 10.0 * root].min + 10.0 * root + 100.0).ceil
          ratios = []
          r = 0.0
          depth.downto(1) do |j|
            r = 1.0 / (2.0 * j / kappa + r)
            ratios[j - 1] = r
          end
          rhos = []
          rho = 1.0
          j = 0
          while j < depth
            rho *= ratios[j]
            break if rho < 1e-19

            rhos << rho
            j += 1
          end
          rhos
        end

        # Integral of the density over the counterclockwise arc of
        # length +d+ (0 < d < 2*PI) starting at relative angle +t0+
        # (already reduced to [-PI, PI] and measured from +mu+).
        #
        # From the Fourier series of the density,
        #   pdf(t) = 1/(2 pi) + (1/pi) * sum rho_j cos(j t),
        # the arc integral is
        #   d/(2 pi) + (2/pi) * sum rho_j cos(j (t0 + d/2)) sin(j d/2) / j.
        # The sine/cosine product form keeps short arcs free of
        # cancellation and makes the mass on both sides of any cut
        # add up to exactly one full period.
        def arc_integral(d, t0, rhos)
          mid = reduce_angle(t0 + 0.5 * d)
          half = 0.5 * d
          sum = 0.0
          rhos.each_with_index do |rho, j|
            n = j + 1
            sum += rho * Math.cos(n * mid) * Math.sin(n * half) / n
          end
          d / TWO_PI + 2.0 * sum / Math::PI
        end
      end
    end
  end
end
