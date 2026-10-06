module Distribution
  module WalleniusHypergeometric
    module Ruby_
      # Result of a computation that involves numerical integration.
      # +status+ is :ok when the adaptive integrator converged and
      # :not_converged when it did not (the returned value is then the
      # best estimate available). +evaluations+ counts integrand calls.
      Result = Struct.new(:value, :status, :evaluations) do
        def success?
          status == :ok
        end
      end

      DEFAULT_TOLERANCE = 1e-10
      DEFAULT_MAX_DEPTH = 24

      class << self
        # Probability mass function of Wallenius' noncentral
        # hypergeometric distribution (two colors).
        #
        # * +k+: number of drawn balls of the interesting class
        # * +m+: number of balls of the interesting class in the urn
        # * +n+: sample size (drawn without replacement)
        # * +total+: total number of balls in the urn
        # * +odds+: weight (bias) of the interesting class relative to
        #   the other class, i.e. w1/w2. +odds+ == 1 reduces to the
        #   ordinary (central) hypergeometric distribution.
        #
        # Returns only the probability. Use +pdf_with_status+ if you
        # need to know whether the numerical integration converged.
        def pdf(k, m, n, total, odds = 1.0, opts = {})
          pdf_with_status(k, m, n, total, odds, opts).value
        end

        # Same as +pdf+ but returns a Result carrying the integration
        # status (:ok / :not_converged).
        def pdf_with_status(k, m, n, total, odds = 1.0, opts = {})
          validate!(m, n, total, odds)
          tol = opts.fetch(:tol, DEFAULT_TOLERANCE)
          max_depth = opts.fetch(:max_depth, DEFAULT_MAX_DEPTH)

          degenerate = degenerate_pdf(k, m, n, total, odds)
          return Result.new(degenerate, :ok, 0) unless degenerate.nil?

          lo, hi = support(m, n, total)
          return Result.new(0.0, :ok, 0) if k < lo || k > hi

          x = k
          y = n - k
          w1 = odds.to_f
          w2 = 1.0
          d = (w1 * (m - x)) + (w2 * ((total - m) - y))
          a1 = w1 / d
          a2 = w2 / d

          integral, ok, evals = wallenius_integral(a1, a2, x, y, tol, max_depth)
          coef = Math.binomial_coefficient(m, x) * Math.binomial_coefficient(total - m, y)
          value = (coef * integral).to_f
          value = 0.0 if value.negative? && value > -1e-12
          Result.new(value, ok ? :ok : :not_converged, evals)
        end

        # Cumulative distribution function: P(X <= k).
        def cdf(k, m, n, total, odds = 1.0, opts = {})
          cdf_with_status(k, m, n, total, odds, opts).value
        end

        # Same as +cdf+ but returns a Result with an aggregated
        # integration status over all summed point probabilities.
        def cdf_with_status(k, m, n, total, odds = 1.0, opts = {})
          validate!(m, n, total, odds)
          lo, hi = support(m, n, total)
          return Result.new(0.0, :ok, 0) if k < lo
          return Result.new(1.0, :ok, 0) if k >= hi

          total_evals = 0
          ok = true
          sum = 0.0
          (lo..k).each do |i|
            r = pdf_with_status(i, m, n, total, odds, opts)
            sum += r.value
            total_evals += r.evaluations
            ok &&= r.success?
          end
          sum = 1.0 if sum > 1.0 && sum < 1.0 + 1e-9
          Result.new(sum, ok ? :ok : :not_converged, total_evals)
        end

        # Expected value of the distribution, computed by summing
        # x * P(X = x) over the whole support.
        def mean(m, n, total, odds = 1.0, opts = {})
          mean_with_status(m, n, total, odds, opts).value
        end

        # Same as +mean+ but returns a Result with an aggregated
        # integration status.
        def mean_with_status(m, n, total, odds = 1.0, opts = {})
          validate!(m, n, total, odds)
          lo, hi = support(m, n, total)

          total_evals = 0
          ok = true
          sum = 0.0
          (lo..hi).each do |i|
            r = pdf_with_status(i, m, n, total, odds, opts)
            sum += i * r.value
            total_evals += r.evaluations
            ok &&= r.success?
          end
          Result.new(sum, ok ? :ok : :not_converged, total_evals)
        end

        # Draw a sample by simulating the sequential (ball-by-ball)
        # Wallenius process: at each of the +n+ draws, a class is picked
        # with probability proportional to (remaining balls of the
        # class) * (weight of the class). A class that is exhausted
        # simply gets probability zero. Returns the number of drawn
        # balls of the interesting class.
        def rng(m, n, total, odds = 1.0, seed = nil)
          validate!(m, n, total, odds)
          seed = Random.new_seed if seed.nil?
          prng = Random.new(seed)

          remaining1 = m
          remaining2 = total - m
          w1 = odds.to_f
          drawn = 0
          n.times do
            break if (remaining1 + remaining2).zero?

            p1 = pick_probability(w1, remaining1, remaining2)
            if prng.rand < p1
              drawn += 1
              remaining1 -= 1
            else
              remaining2 -= 1
            end
          end
          drawn
        end

        private

        def validate!(m, n, total, odds)
          fail ArgumentError, 'total must be non-negative' if total.negative?
          fail ArgumentError, 'm must be between 0 and total' if m.negative? || m > total
          fail ArgumentError, 'n must be between 0 and total' if n.negative? || n > total
          fail ArgumentError, 'odds must be non-negative' if odds.to_f.nan? || odds.negative?
        end

        def support(m, n, total)
          [[0, n - (total - m)].max, [n, m].min]
        end

        # Handles every case that does not need numerical integration.
        # Returns the probability, or nil when the general integral
        # form must be used.
        def degenerate_pdf(k, m, n, total, odds)
          lo, hi = support(m, n, total)
          return nil if k < lo || k > hi
          return (k.zero? ? 1.0 : 0.0) if n.zero? || m.zero?
          return (k == n ? 1.0 : 0.0) if m == total
          return (k == m ? 1.0 : 0.0) if n == total
          return (k == lo ? 1.0 : 0.0) if odds.to_f.zero?
          return (k == hi ? 1.0 : 0.0) if odds.to_f.infinite?
          if odds == 1
            return Distribution::Hypergeometric::Ruby_.pdf(k, m, n, total).to_f
          end
          nil
        end

        # Probability of picking the interesting class in one
        # sequential draw, given the remaining balls of each class.
        def pick_probability(w1, remaining1, remaining2)
          return 0.0 if remaining1.zero? || w1.zero?
          return 1.0 if remaining2.zero? || w1.infinite?
          a = w1 * remaining1
          b = remaining2.to_f
          a / (a + b)
        end

        # Wallenius' pmf integral, computed after the substitution
        # tau = -log(t), which turns the integrand into
        #   (1 - e^(-a1*tau))^x * (1 - e^(-a2*tau))^y * e^(-tau)
        # on [0, TAU_MAX]. This form is smooth, vanishes at both ends
        # and stays well behaved even for extreme odds. The interval
        # is split at the transition points 1/a1 and 1/a2 so the
        # adaptive integrator never straddles a sharp transition.
        TAU_MAX = 30.0

        def wallenius_integral(a1, a2, x, y, tol, max_depth)
          points = [0.0]
          [a1, a2].each do |a|
            bp = 1.0 / a
            points << bp if bp.positive? && bp < TAU_MAX
          end
          points << TAU_MAX
          points = points.uniq.sort

          total = 0.0
          ok_all = true
          evaluations = 0
          tol_each = tol / (points.size - 1)
          points.each_cons(2) do |lo_t, hi_t|
            integrand = lambda do |tau|
              wallenius_integrand(tau, a1, a2, x, y)
            end
            value, ok, evals = adaptive_simpson(integrand, lo_t, hi_t, tol_each, max_depth)
            total += value
            ok_all &&= ok
            evaluations += evals
          end
          [total, ok_all, evaluations]
        end

        # Integrand in tau space, evaluated in the log domain for
        # numerical stability.
        def wallenius_integrand(tau, a1, a2, x, y)
          u1 = -expm1(-a1 * tau)
          u2 = -expm1(-a2 * tau)
          s = -tau
          s += x * Math.log(u1) if x.positive?
          s += y * Math.log(u2) if y.positive?
          s < -745.0 ? 0.0 : Math.exp(s)
        end

        # exp(x) - 1, accurate for small |x|.
        def expm1(x)
          if x.abs < 1e-5
            x * (1.0 + x * (0.5 + x * (1.0 / 6.0 + x / 24.0)))
          else
            Math.exp(x) - 1.0
          end
        end

        # Adaptive Simpson integration of +func+ over [a, b].
        # Returns [value, converged, evaluations].
        def adaptive_simpson(func, a, b, tol, max_depth)
          evaluations = 0
          f = lambda do |t|
            evaluations += 1
            func.call(t)
          end
          fa = f.call(a)
          fb = f.call(b)
          fm = f.call(0.5 * (a + b))
          whole = (b - a) * (fa + 4.0 * fm + fb) / 6.0
          value, ok = simpson_rec(f, a, b, fa, fm, fb, whole, tol, max_depth)
          ok = false if value.nan?
          [value, ok, evaluations]
        end

        def simpson_rec(f, a, b, fa, fm, fb, whole, tol, depth)
          m = 0.5 * (a + b)
          lm = 0.5 * (a + m)
          rm = 0.5 * (m + b)
          flm = f.call(lm)
          frm = f.call(rm)
          left = (m - a) * (fa + 4.0 * flm + fm) / 6.0
          right = (b - m) * (fm + 4.0 * frm + fb) / 6.0
          delta = left + right - whole
          return [left + right + delta / 15.0, true] if delta.abs <= 15.0 * tol
          return [left + right, false] if depth <= 0

          lv, lok = simpson_rec(f, a, m, fa, flm, fm, left, tol * 0.5, depth - 1)
          rv, rok = simpson_rec(f, m, b, fm, frm, fb, right, tol * 0.5, depth - 1)
          [lv + rv, lok && rok]
        end
      end
    end
  end
end
