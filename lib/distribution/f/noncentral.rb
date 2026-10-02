require 'distribution/f/noncentral_mixture'

module Distribution
  module F
    # Non-central Fisher-Snedecor (non-central F) distribution with
    # numerator degrees of freedom +d1+, denominator degrees of freedom
    # +d2+ and non-centrality parameter +lambda+.
    #
    # The CDF is evaluated as a Poisson mixture of regularized incomplete
    # beta functions (see NoncentralMixture), never by shifting a central
    # F.  All four core functions are provided: #pdf, #cdf, #survival and
    # #quantile (inverse CDF).
    #
    # Every call accepts a +tolerance+ (absolute accuracy for pdf, cdf and
    # survival; root accuracy for quantile) and a +max_iterations+ budget
    # (the maximum number of Poisson mixture terms).  Results are returned
    # as Distribution::F::NoncentralResult, a Float that also carries
    # +error_bound+, +iterations+ and +converged?+.  A truncated result is
    # reported explicitly through +converged?+ rather than presented as if
    # the requested tolerance had been reached.
    module Noncentral
      DEFAULT_TOLERANCE = 1e-12
      DEFAULT_MAX_ITERATIONS = 200_000
      # exp(UNDERFLOW_LOG) is safely zero in ordinary floating point.
      UNDERFLOW_LOG = -740.0

      class << self
        # Probability density function at +x+.
        def pdf(x, d1, d2, lambda = 0.0,
                tolerance: DEFAULT_TOLERANCE,
                max_iterations: DEFAULT_MAX_ITERATIONS)
          x = to_float!(x, 'x')
          check_parameters(d1, d2, lambda)
          return exact_result(0.0) if x < 0.0 || x.infinite?
          return central(:pdf, x, d1, d2) if lambda.zero?
          return pdf_at_zero(d1, lambda) if x.zero?

          mu = lambda / 2.0
          mode = mu.floor.to_i
          b = d2 / 2.0
          ln_ratio = Math.log((d1 * x) / d2)
          ln_one_plus = Math.log1p((d1 * x) / d2)

          term = lambda do |j|
            a = d1 / 2.0 + j
            Math.exp(log_density(a, b, d1, d2, ln_ratio, ln_one_plus))
          end
          up_ratio = lambda do |j|
            a = d1 / 2.0 + j
            # p_{j+1}/p_j.  The F-density ratio between shapes a and a+1
            # is (a+b)/a * z, z = d1*x/(d1*x+d2).
            (mu / (j + 1).to_f) * ((a + b) / a) *
              ((d1 * x) / (d1 * x + d2))
          end
          down_ratio = lambda do |j|
            a = d1 / 2.0 + j
            # down_ratio.call(k) = p_{k-1}/p_k.  The F-density ratio from
            # shape a-1 up to a is (a+b-1)/(a-1)*z, so its reciprocal (with
            # the Poisson weight j/mu) is (a-1)/(a+b-1) * (1/z).
            (j.to_f / mu) * ((a - 1.0) / (a + b - 1.0)) *
              ((d1 * x + d2) / (d1 * x))
          end

          value, bound, iterations, converged =
            NoncentralMixture.sum(mu, mode, max_iterations, tolerance,
                                  term, up_ratio, down_ratio)
          F::NoncentralResult.new(value, bound, iterations, converged)
        end

        # Cumulative distribution function P(X <= x).
        def cdf(x, d1, d2, lambda = 0.0,
                tolerance: DEFAULT_TOLERANCE,
                max_iterations: DEFAULT_MAX_ITERATIONS)
          x = to_float!(x, 'x')
          check_parameters(d1, d2, lambda)
          return exact_result(0.0) if x < 0.0 || x.zero?
          return exact_result(1.0) if x.infinite?
          return central(:cdf, x, d1, d2) if lambda.zero?

          mixture_probability(d1, d2, lambda, tolerance, max_iterations, x) do |a, b, z|
            regularized_beta(a, b, z)
          end
        end

        # Survival function P(X > x) = 1 - CDF(x), evaluated directly in
        # the complementary beta form so probabilities close to one keep
        # their accuracy instead of being obtained as 1 - a rounded CDF.
        def survival(x, d1, d2, lambda = 0.0,
                     tolerance: DEFAULT_TOLERANCE,
                     max_iterations: DEFAULT_MAX_ITERATIONS)
          x = to_float!(x, 'x')
          check_parameters(d1, d2, lambda)
          return exact_result(1.0) if x < 0.0 || x.zero?
          return exact_result(0.0) if x.infinite?
          return central(:q, x, d1, d2) if lambda.zero?

          mixture_probability(d1, d2, lambda, tolerance, max_iterations, x) do |a, b, z|
            # I_{1-z}(b, a) = 1 - I_z(a, b), so the result sums to 1 - CDF.
            regularized_beta(b, a, 1.0 - z)
          end
        end

        alias_method :q, :survival

        # Inverse CDF (quantile): return x such that CDF(x) = +p+.
        def quantile(p, d1, d2, lambda = 0.0,
                     tolerance: 1e-10,
                     max_iterations: DEFAULT_MAX_ITERATIONS)
          p = to_float!(p, 'p')
          check_parameters(d1, d2, lambda)
          fail ArgumentError, 'p must be within [0, 1]' unless (0.0..1.0).cover?(p)
          return exact_result(0.0) if p.zero?
          return F::NoncentralResult.new(Float::INFINITY, 0.0, 0, true) if p == 1.0
          return central(:quantile, p, d1, d2) if lambda.zero?

          # Invert survival (not CDF) for p near 1, since 1 - a rounded CDF
          # loses all accuracy in the upper tail.  Target q = min(p, 1-p).
          upper = p > 0.5
          target = upper ? 1.0 - p : p

          tail_at = lambda do |point|
            if upper
              survival(point, d1, d2, lambda,
                       tolerance: [tolerance / 4.0, 1e-14].max,
                       max_iterations: max_iterations)
            else
              cdf(point, d1, d2, lambda,
                  tolerance: [tolerance / 4.0, 1e-14].max,
                  max_iterations: max_iterations)
            end
          end
          pdf_at = lambda do |point|
            pdf(point, d1, d2, lambda,
                tolerance: [tolerance / 4.0, 1e-12].max,
                max_iterations: max_iterations)
          end

          lo = 0.0
          central = central_quantile(p, d1, d2)
          hi = central * (1.0 + lambda / d1.to_f)
          hi = Float::EPSILON unless hi > 0.0
          hi_result = tail_at.call(hi)
          all_converged = hi_result.converged?
          total_iterations = hi_result.iterations

          # Bracketing: for the lower tail the upper bound must have
          # tail >= target; for the upper tail it must have tail <= target.
          256.times do
            enough = upper ? hi_result.to_f <= target : hi_result.to_f >= target
            break if enough
            hi *= 2.0
            hi_result = tail_at.call(hi)
            break unless hi_result.to_f.finite?
            all_converged &&= hi_result.converged?
            total_iterations += hi_result.iterations
          end
          enough = upper ? hi_result.to_f <= target : hi_result.to_f >= target
          unless enough
            return F::NoncentralResult.new(hi, Float::INFINITY,
                                         total_iterations, false)
          end

          x = central.clamp(hi / 2.0, hi)
          gap = Float::INFINITY

          100.times do
            fx = tail_at.call(x)
            density = pdf_at.call(x)
            total_iterations += fx.iterations + density.iterations
            all_converged &&= fx.converged? && density.converged?
            # Lower tail: gap = CDF - target (root when zero).
            # Upper tail: gap = survival - target, decreasing in x; the
            # Newton step is x - gap/density in both cases.
            gap = fx.to_f - target

            if upper ? gap >= 0.0 : gap <= 0.0
              lo = x
            else
              hi = x
            end

            width_error = (hi - lo) / 2.0
            if width_error <= tolerance * (1.0 + x) && all_converged
              return F::NoncentralResult.new(x.clamp(lo, hi), width_error,
                                             total_iterations, true)
            end

            candidate =
              if density > 0.0 && density.finite? && gap.finite?
                x - gap / density.to_f
              else
                nil
              end
            if candidate.nil? || !(candidate > lo) || !(candidate < hi)
              candidate = (lo + hi) / 2.0
            end
            x = candidate
          end

          width_error = (hi - lo) / 2.0
          F::NoncentralResult.new((lo + hi) / 2.0,
                                   [width_error, gap.abs].max,
                                   total_iterations, false)
        end

        alias_method :p_value, :quantile
        alias_method :inverse_cdf, :quantile

        private

        def check_parameters(d1, d2, lambda)
          unless d1.is_a?(Numeric) && d1.to_f.finite? && d1 > 0
            fail ArgumentError, 'numerator degrees of freedom d1 must be positive'
          end
          unless d2.is_a?(Numeric) && d2.to_f.finite? && d2 > 0
            fail ArgumentError, 'denominator degrees of freedom d2 must be positive'
          end
          unless lambda.is_a?(Numeric) && lambda.to_f.finite? && lambda >= 0
            fail ArgumentError, 'non-centrality lambda must be non-negative'
          end
        end

        def to_float!(value, name)
          value = value.to_f
          fail ArgumentError, "#{name} must be a real number" if value.nan?
          value
        end

        def exact_result(value)
          F::NoncentralResult.new(value, 0.0, 0, true)
        end

        def central(method, *args)
          value = Distribution::F::Ruby_.send(method, *args).to_f
          bound = method == :quantile ? 8.88e-16 * (1.0 + value.abs) : 4.0 * Float::EPSILON
          F::NoncentralResult.new(value, bound, 0, true)
        end

        def central_quantile(p, d1, d2)
          Distribution::F::Ruby_.quantile(p, d1, d2).to_f
        end

        # Density at x = 0: the limiting value.  It is infinite when
        # d1 < 2, exp(-lambda/2) when d1 = 2 and zero when d1 > 2.
        def pdf_at_zero(d1, lambda)
          value =
            if d1 < 2.0
              Float::INFINITY
            elsif d1 == 2.0
              Math.exp(-lambda / 2.0)
            else
              0.0
            end
          exact_result(value)
        end

        def log_density(a, b, d1, d2, ln_ratio, ln_one_plus)
          Math.log(d1.to_f / d2) - Math.logbeta(a, b) +
            (a - 1.0) * ln_ratio - (a + b) * ln_one_plus
        end

        def mixture_probability(d1, d2, lambda, tolerance, max_iterations, x)
          mu = lambda / 2.0
          mode = mu.floor.to_i
          b = d2 / 2.0
          z = (d1 * x) / (d1 * x + d2)

          term = lambda do |j|
            a = d1 / 2.0 + j
            yield(a, b, z)
          end

          value, bound, iterations, converged =
            NoncentralMixture.sum(mu, mode, max_iterations, tolerance, term)
          F::NoncentralResult.new(clamp_probability(value, bound),
                               bound, iterations, converged)
        end

        def clamp_probability(value, bound)
          if value > 1.0 && value <= 1.0 + bound
            1.0
          elsif value < 0.0 && value >= -bound
            0.0
          else
            value
          end
        end

        # Regularized incomplete beta I_z(a, b), with an underflow short
        # circuit that avoids the continued-fraction machinery when the
        # integrand prefactor already underflows.
        def regularized_beta(a, b, z)
          return 0.0 if z <= 0.0
          return 1.0 if z >= 1.0

          ln_pre = -Math.logbeta(a, b) +
                   a * Math.log(z) + b * Math::Log.log1p(-z)
          if ln_pre < UNDERFLOW_LOG
            return z < (a + 1.0) / (a + b + 2.0) ? 0.0 : 1.0
          end

          MathExtension::IncompleteBeta.evaluate(a, b, z)
        end

      end
    end

    # A numeric result that also reports the numerical error bound, the
    # mixture budget consumed and whether the requested tolerance was met.
    # It compares and coerces like a Float, so it works directly in
    # arithmetic and numeric RSpec matchers; #to_f gives the bare value.
    class NoncentralResult < Numeric
      attr_reader :value, :error_bound, :iterations

      def initialize(value, error_bound, iterations, converged)
        @value = value.to_f
        @error_bound = error_bound.to_f
        @iterations = iterations.to_i
        @converged = converged
      end

      def converged?
        @converged
      end

      def <=>(other)
        @value <=> other.to_f
      end

      def to_f
        @value
      end

      def to_i
        @value.to_i
      end

      def finite?
        @value.finite?
      end

      def infinite?
        @value.infinite?
      end

      def nan?
        @value.nan?
      end

      def abs
        @value.abs
      end

      # Delegate arithmetic and other Float behaviour to the wrapped value,
      # while comparisons remain exact through <tt><=></tt>.
      def method_missing(name, *args, &block)
        if @value.respond_to?(name)
          @value.public_send(name, *args, &block)
        else
          super
        end
      end

      def respond_to_missing?(name, include_private = false)
        @value.respond_to?(name, include_private) || super
      end

      def coerce(other)
        [other.to_f, @value]
      end

      def to_s
        format('%g (error_bound=%g, iterations=%d, converged=%s)',
               @value, @error_bound, @iterations, @converged)
      end
      alias_method :inspect, :to_s
    end
  end
end
