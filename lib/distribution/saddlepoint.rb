module Distribution
  # Saddlepoint approximation of tail probabilities (and densities) for
  # distributions known only through their cumulant generating function
  # (CGF) and its derivatives, when no closed-form pdf/cdf is available.
  #
  # The user supplies callbacks for K(t), K'(t) and K''(t) (K'''(t) is
  # optional and improves accuracy near the mean), declares whether the
  # distribution is continuous or lattice-valued (with its span), and may
  # restrict the domain on which the CGF is defined.
  #
  # Tail probabilities use the Lugannani-Rice formula; lattice
  # distributions use a continuity correction of half the span.
  # Near the mean (where the formula degenerates to 0/0) a series
  # expansion in the saddlepoint is used instead. If the saddlepoint
  # equation has no root inside the CGF domain, the result is flagged
  # as not converged instead of forcing a root.
  #
  # All results are *approximations*: the returned Result object carries
  # +approximate? == true+ plus solver diagnostics (saddlepoint, number
  # of iterations, convergence flag, messages).
  #
  # == Example (Gamma(shape, scale), continuous)
  #   spa = Distribution::Saddlepoint.new(
  #     cgf:        ->(t) { -3.0 * Math.log(1 - 2.0 * t) },
  #     cgf_prime:  ->(t) { 6.0 / (1 - 2.0 * t) },
  #     cgf_second: ->(t) { 12.0 / (1 - 2.0 * t)**2 },
  #     domain: -Float::INFINITY...0.5
  #   )
  #   spa.cdf(10.0).value # => approximation of P(X <= 10)
  #
  # == Example (Binomial(n, p), lattice with span 1)
  #   spa = Distribution::Saddlepoint.new(
  #     cgf:        ->(t) { 20 * Math.log(0.7 + 0.3 * Math.exp(t)) },
  #     cgf_prime:  ->(t) { 20 * 0.3 * Math.exp(t) / (0.7 + 0.3 * Math.exp(t)) },
  #     cgf_second: ->(t) { ... },
  #     lattice: 1
  #   )
  #   spa.cdf(8).value # => approximation of P(X <= 8)
  #
  # == References
  # * Lugannani, R. & Rice, S. (1980). Saddle point approximation for the
  #   distribution of the sum of independent random variables.
  # * Butler, R. W. (2007). Saddlepoint Approximations with Applications.
  class Saddlepoint
    # Default tolerance for the saddlepoint equation K'(t) = x.
    TOLERANCE = 1e-10
    # Maximum number of safeguarded Newton iterations.
    MAX_ITERATIONS = 100
    # Saddlepoints with |t| below this are treated with series expansions
    # to avoid the 0/0 cancellation in the Lugannani-Rice formula.
    NEAR_MEAN_THRESHOLD = 1e-4

    # Wraps an approximate value together with solver diagnostics.
    class Result
      # Approximated value (nil when the approximation could not be computed).
      attr_reader :value
      # Saddlepoint t solving K'(t) = x (nil when not found).
      attr_reader :saddlepoint
      # Number of solver iterations performed.
      attr_reader :iterations
      # Lugannani-Rice quantities (nil when not applicable).
      attr_reader :w, :u
      # Diagnostic messages (e.g. why the solver failed).
      attr_reader :messages

      def initialize(value:, saddlepoint:, iterations:, converged:, w: nil, u: nil, messages: [])
        @value = value
        @saddlepoint = saddlepoint
        @iterations = iterations
        @converged = converged
        @w = w
        @u = u
        @messages = messages
      end

      # Always true: every value produced by Saddlepoint is an approximation.
      def approximate?
        true
      end

      # Whether the saddlepoint equation was solved inside the CGF domain.
      def converged?
        @converged
      end

      def to_f
        @value && @value.to_f
      end

      # Solver diagnostics as a hash.
      def diagnostics
        {
          approximate: true,
          converged: @converged,
          saddlepoint: @saddlepoint,
          iterations: @iterations,
          w: @w,
          u: @u,
          messages: @messages.dup
        }
      end
    end

    # +cgf+::        callable K(t) (required)
    # +cgf_prime+::  callable K'(t) (required)
    # +cgf_second+:: callable K''(t) (required)
    # +cgf_third+::  callable K'''(t) (optional; estimated numerically
    #                from K'' when omitted)
    # +lattice+::    +nil+/+false+ for a continuous distribution, or the
    #                positive span of a lattice distribution (e.g. 1 for
    #                integer-valued ones)
    # +domain+::     interval (Range or [min, max]) on which the CGF is
    #                defined; defaults to the whole real line
    def initialize(cgf:, cgf_prime:, cgf_second:, cgf_third: nil,
                   lattice: nil, domain: nil,
                   tolerance: TOLERANCE, max_iterations: MAX_ITERATIONS)
      @k = cgf
      @k1 = cgf_prime
      @k2 = cgf_second
      @k3 = cgf_third
      @lattice = lattice && lattice.to_f
      raise ArgumentError, 'lattice span must be positive' if @lattice && @lattice <= 0
      @t_min, @t_max = normalize_domain(domain)
      raise ArgumentError, 'domain must contain 0 (K(0) = 0)' unless @t_min < 0 && @t_max > 0
      @tolerance = tolerance
      @max_iterations = max_iterations
    end

    # True for lattice (discrete) distributions, false for continuous ones.
    def lattice?
      !@lattice.nil?
    end

    # Mean of the distribution, K'(0).
    def mean
      @k1.call(0.0)
    end

    # Approximate P(X <= x). For lattice distributions a continuity
    # correction of half the span is applied.
    def cdf(x)
      x = continuity_correct(x)
      tail_probability(x, :cdf)
    end

    # Approximate P(X > x) (survival function), computed tail-stably.
    def sf(x)
      x = continuity_correct(x)
      tail_probability(x, :sf)
    end

    # Saddlepoint density approximation (continuous distributions).
    def pdf(x)
      raise 'pdf is only defined for continuous distributions; use pmf' if lattice?
      density(x)
    end

    # Saddlepoint point-mass approximation (lattice distributions).
    def pmf(x)
      raise 'pmf is only defined for lattice distributions; use pdf' unless lattice?
      density(x, mass: true)
    end

    private

    def normalize_domain(domain)
      case domain
      when nil
        [-Float::INFINITY, Float::INFINITY]
      when Range
        [domain.begin.to_f, domain.end.to_f]
      when Array
        [domain[0].to_f, domain[1].to_f]
      else
        raise ArgumentError, 'domain must be a Range, a [min, max] pair or nil'
      end
    end

    def continuity_correct(x)
      return x.to_f unless lattice?
      # P(X <= x) and P(X > x) both use x + span/2 for lattice spans.
      x.to_f + @lattice / 2.0
    end

    # Third derivative at 0, from the callback or a central difference
    # of K''. Used only by the near-mean series expansions.
    def third_cumulant
      return @k3.call(0.0) if @k3
      h = 1e-3
      h = [h, (@t_max - @t_min) / 4.0].min if (@t_max - @t_min).finite?
      (@k2.call(h) - @k2.call(-h)) / (2.0 * h)
    end

    # Safeguarded Newton solve of K'(t) = x on the CGF domain.
    # Returns [t, iterations, messages]; t is nil when no root exists
    # inside the domain (e.g. x outside the attainable mean range).
    def solve_saddlepoint(x)
      mu = @k1.call(0.0)
      return [0.0, 0, []] if x == mu

      messages = []
      lo = @t_min
      hi = @t_max
      t = initial_guess(x, mu)
      iterations = 0

      while iterations < @max_iterations
        iterations += 1
        g = @k1.call(t) - x
        unless g.finite?
          messages << "K'(t) not finite at t=#{t}; x=#{x} may lie outside " +                      'the attainable mean range of the CGF'
          return [nil, iterations, messages]
        end
        return [t, iterations, messages] if g.abs <= @tolerance * (1.0 + x.abs)

        if g > 0
          hi = t
        else
          lo = t
        end

        curvature = @k2.call(t)
        step = curvature > 0 ? g / curvature : Float::INFINITY
        t_new = t - step
        # Safeguard: fall back to bisection (or a bounded expansion when
        # the bracket end is infinite) whenever Newton leaves the bracket.
        t_new = next_bracket_point(lo, hi) unless t_new > lo && t_new < hi
        unless t_new.finite? && t_new != t
          # The iteration stagnated without meeting the residual
          # tolerance: there is no root of K'(t) = x inside the domain
          # (e.g. x outside the support of the distribution).
          messages << "no root of K'(t) = #{x} inside the CGF domain " +                      '(x may lie outside the support of the distribution)'
          return [nil, iterations, messages]
        end
        t = t_new
      end

      messages << 'saddlepoint solver did not converge within '                   "#{@max_iterations} iterations"
      [nil, iterations, messages]
    rescue StandardError => e
      [nil, iterations, messages + ["saddlepoint solver failed: #{e.message}"]]
    end

    def initial_guess(x, mu)
      # Small step from 0 in the direction of the root, kept inside domain.
      direction = x > mu ? 1.0 : -1.0
      bound = direction > 0 ? @t_max : @t_min
      guess = direction * 0.1
      guess = bound / 2.0 if bound.finite? && guess.abs > bound.abs / 2.0
      guess
    end

    def next_bracket_point(lo, hi)
      if lo.finite? && hi.finite?
        (lo + hi) / 2.0
      elsif hi.finite?
        hi - 1.0
      elsif lo.finite?
        lo + 1.0
      else
        0.0
      end
    end

    # Signed likelihood-ratio quantity w and Lugannani-Rice u.
    # Near t = 0, t*x - K(t) is computed by series to avoid cancellation,
    # and the 1/w - 1/u correction is replaced by its limit.
    def compute_w_u(t, x)
      k2 = @k2.call(t)
      if t.abs < NEAR_MEAN_THRESHOLD
        k3 = third_cumulant
        a = 0.5 * k2 * t * t + (k3 / 3.0) * t * t * t
        root = Math.sqrt([2.0 * a, 0.0].max)
        w = t >= 0 ? root : -root
        u = t * Math.sqrt(k2)
        # Series limit of (1/w - 1/u), avoiding 0/0.
        correction = k3 / (6.0 * k2**1.5)
        return [w, u, correction]
      end
      a = t * x - @k.call(t)
      root = Math.sqrt([2.0 * a, 0.0].max)
      w = t >= 0 ? root : -root
      u = t * Math.sqrt(k2)
      [w, u, nil]
    end

    def tail_probability(x, kind)
      t, iterations, messages = solve_saddlepoint(x)
      unless t
        return Result.new(value: nil, saddlepoint: nil, iterations: iterations,
                          converged: false,
                          messages: messages + ["no saddlepoint for x=#{x} inside the CGF domain"])
      end

      w, u, correction_limit = compute_w_u(t, x)
      phi_w = Math.exp(-0.5 * w * w) / SQ2PI
      correction = correction_limit || (1.0 / w - 1.0 / u)

      value =
        if kind == :cdf
          normal_cdf(w) + phi_w * correction
        else
          normal_cdf(-w) - phi_w * correction
        end
      value = 0.0 if value < 0 && value > -1e-10
      value = 1.0 if value > 1 && value < 1 + 1e-10

      Result.new(value: value, saddlepoint: t, iterations: iterations,
                 converged: true, w: w, u: u, messages: messages)
    end

    def density(x, mass: false)
      x = x.to_f
      t, iterations, messages = solve_saddlepoint(x)
      unless t
        return Result.new(value: nil, saddlepoint: nil, iterations: iterations,
                          converged: false,
                          messages: messages + ["no saddlepoint for x=#{x} inside the CGF domain"])
      end
      k2 = @k2.call(t)
      log_f = @k.call(t) - t * x - 0.5 * Math.log(2.0 * Math::PI * k2)
      value = Math.exp(log_f)
      value *= @lattice if mass
      Result.new(value: value, saddlepoint: t, iterations: iterations,
                 converged: true, messages: messages)
    end

    def normal_cdf(z)
      0.5 * (1.0 + Math.erf(z / Math.sqrt(2.0)))
    end
  end
end
