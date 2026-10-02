module Distribution
  module PoissonBinomial
    # Pure Ruby engine for the Poisson-binomial distribution.
    #
    # == Numerical notes
    #
    # The pmf is computed by iterative convolution of the Bernoulli
    # probabilities (O(n^2) time, O(n) memory). Each convolution step
    # introduces a relative error of a few machine epsilons, so the
    # accumulated error grows slowly with n (roughly O(n * 2^-53)).
    # The resulting pmf is renormalized so its total mass is exactly 1,
    # preserving probability mass normalization.
    #
    # Probabilities below ~1e-308 underflow to 0.0 in the far tails;
    # this only affects individual masses that are already below double
    # precision resolution.
    #
    # The right tail (sf) is computed by summing pmf terms from the
    # largest k downwards, so tiny tail probabilities (down to ~1e-308)
    # stay accurate. Computing it as 1 - cdf would lose all precision
    # below ~1e-16 because of floating point cancellation.
    #
    # == Input limits
    #
    # * +probs+ must be an Array of real numbers in [0, 1].
    #   0 and 1 are treated as deterministic events (they never / always
    #   fail) and are handled exactly.
    # * An empty Array is allowed and represents a point mass at 0.
    # * Practical size is n up to ~10^4 (O(n^2) running time).
    module Ruby_
      class << self
        # Full pmf as an Array of Floats: entry +i+ is P(X = i).
        # The array is normalized so that its sum is 1.
        def pmf_array(probs)
          ps = validate_probs(probs)
          deterministic_ones = ps.count(1.0)
          pmf = [1.0]
          ps.each do |p|
            next if p == 0.0 || p == 1.0
            q = 1.0 - p
            nxt = Array.new(pmf.size + 1, 0.0)
            pmf.each_with_index do |mass, i|
              nxt[i] += mass * q
              nxt[i + 1] += mass * p
            end
            pmf = nxt
          end
          total = pmf.inject(0.0) { |ac, m| ac + m }
          pmf = pmf.map { |m| m / total } if total > 0.0 && total != 1.0
          Array.new(deterministic_ones, 0.0).concat(pmf)
        end

        # P(X = k)
        def pdf(k, probs)
          return 0.0 if k != k.to_i
          pmf = pmf_array(probs)
          i = k.to_i
          return 0.0 if i < 0 || i >= pmf.size
          pmf[i]
        end

        # P(X <= k)
        def cdf(k, probs)
          pmf = pmf_array(probs)
          i = k.floor
          return 0.0 if i < 0
          return 1.0 if i >= pmf.size - 1
          ac = 0.0
          (0..i).each { |j| ac += pmf[j] }
          ac > 1.0 ? 1.0 : ac
        end

        # Right tail: P(X >= k).
        # Computed by summing the pmf downwards from the largest k,
        # so very small tail probabilities remain accurate (a plain
        # 1 - cdf subtraction would underflow to 0 below ~1e-16).
        def sf(k, probs)
          pmf = pmf_array(probs)
          i = k.ceil
          return 1.0 if i <= 0
          return 0.0 if i >= pmf.size
          ac = 0.0
          j = pmf.size - 1
          while j >= i
            ac += pmf[j]
            j -= 1
          end
          ac > 1.0 ? 1.0 : ac
        end

        # Smallest k with P(X <= k) >= qn
        def p_value(qn, probs)
          unless qn.is_a?(Numeric) && qn >= 0 && qn <= 1
            fail ArgumentError, "quantile must be a real number in [0, 1], got #{qn.inspect}"
          end
          pmf = pmf_array(probs)
          ac = 0.0
          pmf.each_with_index do |m, i|
            ac += m
            return i if qn <= ac
          end
          pmf.size - 1
        end

        alias_method :quantile, :p_value

        private

        def validate_probs(probs)
          unless probs.is_a?(Array)
            fail ArgumentError, "expected an Array of Bernoulli probabilities, got #{probs.inspect}"
          end
          probs.map do |p|
            unless p.is_a?(Numeric) && !p.is_a?(Complex) && p >= 0 && p <= 1
              fail ArgumentError, "invalid probability #{p.inspect}: must be a real number in [0, 1]"
            end
            Float(p)
          end
        end
      end
    end
  end
end
