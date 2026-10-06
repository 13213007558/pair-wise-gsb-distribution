module Distribution
  # Adaptive Rejection Sampling (ARS) for univariate log-concave
  # distributions, after Gilks & Wild (1992, Applied Statistics 41:337-348).
  #
  # The sampler only needs an unnormalized log-density and its derivative;
  # the normalizing constant is never required. Tangent lines at a set of
  # abscissae form a piecewise-exponential upper hull, and chords between
  # adjacent abscissae form a squeeze (lower hull). Every time the
  # log-density has to be evaluated at a rejected (or hull-evaluated)
  # point, that point is added to the abscissa set and the hull is refined.
  #
  # The upper hull is checked for integrability on every update, and the
  # log-concavity assumptions (non-increasing derivative, tangents above
  # the density, chords below it) are verified both at initialization and
  # whenever a new point is evaluated. Any violation raises
  # NotLogConcaveError / NonIntegrableHullError immediately instead of
  # silently returning samples from a wrong distribution.
  #
  # Bounded, one-sided unbounded and two-sided unbounded supports are
  # handled through +-Float::INFINITY bounds.
  #
  # == Usage
  #   sampler = Distribution::AdaptiveRejectionSampler.new(
  #     log_pdf:  ->(x) { -0.5 * x * x }, # unnormalized N(0,1) log-density
  #     dlog_pdf: ->(x) { -x },
  #     lower: -Float::INFINITY,
  #     upper: Float::INFINITY,
  #     points: [-1.0, 0.0, 1.0],
  #     rng: Random.new(42),
  #     max_attempts: 1000
  #   )
  #   sampler.sample    # => one Float
  #   sampler.sample(5) # => Array of 5 Floats
  class AdaptiveRejectionSampler
    # Base class for every error raised by the sampler.
    class Error < StandardError; end

    # Raised when the supplied log-density or its derivative violates the
    # log-concavity assumptions (non-increasing derivative, tangent lines
    # above the density, chords below it).
    class NotLogConcaveError < Error; end

    # Raised when the tangent upper hull is not integrable over the
    # support (e.g. a non-positive tangent slope on a left-unbounded
    # interval), so no valid rejection envelope exists.
    class NonIntegrableHullError < Error; end

    # Raised when a single call to #sample exceeds the caller-supplied
    # maximum number of proposals.
    class MaxAttemptsExceededError < Error; end

    # Relative tolerance for derivative monotonicity checks.
    DERIV_TOL = 1e-8
    # Relative tolerance for hull dominance checks.
    HULL_TOL = 1e-8

    attr_reader :lower, :upper, :max_attempts

    # Number of abscissae currently defining the hulls.
    attr_reader :size

    # +log_pdf+:   callable returning the (unnormalized) log-density.
    # +dlog_pdf+:  callable returning the derivative of +log_pdf+.
    # +lower+:     lower bound of the support (-Float::INFINITY allowed).
    # +upper+:     upper bound of the support (Float::INFINITY allowed).
    # +points+:    initial abscissae, strictly inside the support. If the
    #              support is left-unbounded, the derivative at the
    #              leftmost point must be positive; if right-unbounded,
    #              the derivative at the rightmost point must be negative.
    # +rng+:       random source, must respond to +rand+ like ::Random.
    # +max_attempts+: maximum proposals per #sample call before giving up.
    def initialize(log_pdf:, dlog_pdf:, lower:, upper:, points:, rng: Random.new, max_attempts: 1000)
      raise ArgumentError, 'log_pdf must respond to #call' unless log_pdf.respond_to?(:call)
      raise ArgumentError, 'dlog_pdf must respond to #call' unless dlog_pdf.respond_to?(:call)
      raise ArgumentError, 'rng must respond to #rand' unless rng.respond_to?(:rand)

      lower = Float(lower)
      upper = Float(upper)
      raise ArgumentError, 'bounds must not be NaN' if lower.nan? || upper.nan?
      raise ArgumentError, 'lower bound must be smaller than upper bound' unless lower < upper

      points = points.map { |p| Float(p) }.sort
      raise ArgumentError, 'at least one starting point is required' if points.empty?
      if points.any?(&:nan?) || points.any? { |p| p <= lower || p >= upper }
        raise ArgumentError, 'starting points must be finite and strictly inside the support'
      end
      raise ArgumentError, 'starting points must be distinct' if points.each_cons(2).any? { |a, b| a == b }

      max_attempts = Integer(max_attempts)
      raise ArgumentError, 'max_attempts must be positive' if max_attempts < 1

      @log_pdf = log_pdf
      @dlog_pdf = dlog_pdf
      @lower = lower
      @upper = upper
      @rng = rng
      @max_attempts = max_attempts
      @failure = nil

      @xs = points
      @hs = points.map { |x| checked_log_pdf(x) }
      @dhs = points.map { |x| checked_dlog_pdf(x) }
      rebuild_hull
      @size = @xs.size
    end

    # Draws one sample, or an Array of +n+ samples when +n+ is given.
    def sample(n = nil)
      return draw if n.nil?

      n = Integer(n)
      raise ArgumentError, 'sample size must not be negative' if n.negative?

      Array.new(n) { draw }
    end

    private

    def draw
      raise @failure if @failure

      attempts = 0
      loop do
        attempts += 1
        if attempts > @max_attempts
          fail_with(MaxAttemptsExceededError,
                    "no accepted proposal after #{@max_attempts} attempts")
        end

        segment = select_segment(rand_unit)
        x = sample_segment(segment, rand_unit)
        ux = @hs[segment] + @dhs[segment] * (x - @xs[segment])
        lx = lower_hull_at(x)
        log_w = Math.log(rand_unit)

        # Squeeze test: cheap acceptance that never touches the density
        # and therefore cannot distort the target distribution.
        return x if log_w <= lx - ux

        hx = @log_pdf.call(x)
        if hx.nan? || hx == Float::INFINITY
          fail_with(NotLogConcaveError, "log_pdf returned #{hx} at x=#{x}")
        end

        unless hx == -Float::INFINITY
          # The tangent hull must dominate the density and the chord
          # squeeze must lie below it; otherwise the density is not
          # log-concave and sampling must stop.
          if hx > ux + HULL_TOL * (1.0 + ux.abs)
            fail_with(NotLogConcaveError,
                      "tangent hull violated at x=#{x}: log_pdf=#{hx} above hull=#{ux}")
          end
          if lx > -Float::INFINITY && hx < lx - HULL_TOL * (1.0 + lx.abs)
            fail_with(NotLogConcaveError,
                      "squeeze hull violated at x=#{x}: log_pdf=#{hx} below chord=#{lx}")
          end
          insert_point(x, hx, checked_dlog_pdf(x))
        end

        return x if log_w <= hx - ux
      end
    rescue Error => e
      @failure ||= e
      raise
    end

    def fail_with(klass, message = nil)
      raise klass, message || klass.name
    end

    # Uniform draw on (0, 1), clamped away from the exact endpoints so
    # that logarithms and inverse-CDF transforms stay finite.
    def rand_unit
      v = @rng.rand.to_f
      v = Float::MIN if v <= 0.0
      v = 1.0 - Float::EPSILON if v >= 1.0
      v
    end

    def checked_log_pdf(x)
      h = @log_pdf.call(x)
      unless h.is_a?(Numeric) && h.finite?
        raise ArgumentError, "log_pdf must be finite at every starting point (got #{h} at x=#{x})"
      end

      h.to_f
    end

    def checked_dlog_pdf(x)
      dh = @dlog_pdf.call(x)
      if dh.nil? || dh != dh || dh == Float::INFINITY || dh == -Float::INFINITY
        fail_with(NotLogConcaveError, "dlog_pdf is not finite at x=#{x}")
      end

      dh.to_f
    end

    def insert_point(x, hx, dhx)
      i = @xs.bsearch_index { |xi| xi >= x } || @xs.size
      return if i < @xs.size && @xs[i] == x

      @xs.insert(i, x)
      @hs.insert(i, hx)
      @dhs.insert(i, dhx)
      rebuild_hull
      @size = @xs.size
    end

    # Recomputes tangent intersections, per-segment areas and the
    # cumulative segment distribution. Validates the log-concavity and
    # integrability assumptions on every call.
    def rebuild_hull
      k = @xs.size
      @dhs.each_cons(2).with_index do |(d1, d2), i|
        next if d2 <= d1 + DERIV_TOL * (1.0 + d1.abs + d2.abs)

        fail_with(NotLogConcaveError,
                  "dlog_pdf is not non-increasing between x=#{@xs[i]} and x=#{@xs[i + 1]}")
      end

      @zs = Array.new(k + 1)
      @zs[0] = @lower
      @zs[k] = @upper
      (0...k - 1).each { |i| @zs[i + 1] = tangent_intersection(i) }

      log_areas = Array.new(k) { |i| segment_log_area(i) }
      max = log_areas.max
      if max.nil? || max.nan? || max == Float::INFINITY || max == -Float::INFINITY
        fail_with(NonIntegrableHullError, 'tangent upper hull has zero or infinite total area')
      end

      log_total = max + Math.log(log_areas.sum { |la| la == -Float::INFINITY ? 0.0 : Math.exp(la - max) })
      cumulative = 0.0
      @cum = log_areas.map do |la|
        cumulative += la == -Float::INFINITY ? 0.0 : Math.exp(la - log_total)
        cumulative
      end
      @cum[-1] = 1.0
    end

    # Intersection of the tangents at abscissae i and i+1. For a
    # log-concave density it always lies between the two abscissae.
    def tangent_intersection(i)
      x1, x2 = @xs[i], @xs[i + 1]
      h1, h2 = @hs[i], @hs[i + 1]
      d1, d2 = @dhs[i], @dhs[i + 1]

      if (d1 - d2).abs <= DERIV_TOL * (1.0 + d1.abs + d2.abs)
        # Parallel tangents are only consistent with a linear stretch of
        # the log-density; otherwise the function is not concave.
        unless (h2 - h1 - d1 * (x2 - x1)).abs <= HULL_TOL * (1.0 + h1.abs + h2.abs)
          fail_with(NotLogConcaveError,
                    "inconsistent parallel tangents between x=#{x1} and x=#{x2}")
        end
        return 0.5 * (x1 + x2)
      end

      z = (h2 - h1 - x2 * d2 + x1 * d1) / (d1 - d2)
      tol = HULL_TOL * (1.0 + x1.abs + x2.abs)
      if z < x1 - tol || z > x2 + tol
        fail_with(NotLogConcaveError,
                  "tangent intersection #{z} outside [#{x1}, #{x2}]")
      end
      z.clamp(x1, x2)
    end

    # Log of the (unnormalized) area under one exponential hull segment.
    def segment_log_area(i)
      a = @zs[i]
      b = @zs[i + 1]
      m = @dhs[i]
      c = @hs[i] - m * @xs[i]

      if m.zero?
        if a == -Float::INFINITY || b == Float::INFINITY
          fail_with(NonIntegrableHullError, 'flat tangent on an unbounded interval is not integrable')
        end
        return c + Math.log(b - a)
      end
      if a == -Float::INFINITY
        if m <= 0.0
          fail_with(NonIntegrableHullError,
                    'leftmost tangent must have positive slope on a left-unbounded support')
        end
        return c + m * b - Math.log(m)
      end
      if b == Float::INFINITY
        if m >= 0.0
          fail_with(NonIntegrableHullError,
                    'rightmost tangent must have negative slope on a right-unbounded support')
        end
        return c + m * a - Math.log(-m)
      end

      if m.positive?
        c + m * b + Math.log(-Math.expm1(m * (a - b))) - Math.log(m)
      else
        c + m * a + Math.log(-Math.expm1(m * (b - a))) - Math.log(-m)
      end
    end

    def select_segment(u)
      @cum.bsearch_index { |c| c >= u } || @cum.size - 1
    end

    # Inverse-CDF draw from one truncated exponential hull segment.
    def sample_segment(i, u)
      a = @zs[i]
      b = @zs[i + 1]
      m = @dhs[i]

      return a + u * (b - a) if m.zero?

      if m.positive?
        return b + Math.log(u) / m if a == -Float::INFINITY

        b + Math.log(u + (1.0 - u) * Math.exp(m * (a - b))) / m
      else
        return a + Math.log(1.0 - u) / m if b == Float::INFINITY

        a + Math.log((1.0 - u) + u * Math.exp(m * (b - a))) / m
      end
    end

    # Piecewise-linear squeeze (chords between adjacent abscissae);
    # -infinity outside the outermost abscissae.
    def lower_hull_at(x)
      return -Float::INFINITY if @xs.size < 2 || x < @xs.first || x > @xs.last

      i = @xs.bsearch_index { |xi| xi >= x } || (@xs.size - 1)
      i -= 1 if i.positive? && @xs[i] > x
      i = @xs.size - 2 if i >= @xs.size - 1
      x1, x2 = @xs[i], @xs[i + 1]
      h1, h2 = @hs[i], @hs[i + 1]
      h1 + (h2 - h1) * (x - x1) / (x2 - x1)
    end
  end
end
