module Distribution
  module FisherNoncentralHypergeometric
    module Ruby_
      class << self
        # Lower and upper bounds of the support, fixed by the margins:
        # k in [max(0, n - m2), min(n, m1)]
        def support(m1, m2, n)
          [[0, n - m2].max, [n, m1].min]
        end

        # Fisher's noncentral hypergeometric probability mass function.
        #
        # Probability of drawing +k+ type-1 items in a sample of +n+ drawn
        # without replacement from a population of +m1+ type-1 and +m2+
        # type-2 items, where each type-1 item is +odds+ times as likely
        # to be sampled as a type-2 item (Fisher's definition: the pmf is
        # proportional to C(m1,k) * C(m2,n-k) * odds**k).
        def pdf(k, m1, m2, n, odds)
          validate_parameters(m1, m2, n, odds)
          lo, hi = support(m1, m2, n)
          return 0.0 if k < lo || k > hi
          probabilities(m1, m2, n, odds)[k - lo]
        end

        alias_method :exact_pdf, :pdf

        # Cumulative distribution function: P(X <= k)
        def cdf(k, m1, m2, n, odds)
          validate_parameters(m1, m2, n, odds)
          lo, hi = support(m1, m2, n)
          return 0.0 if k < lo
          return 1.0 if k >= hi
          probs = probabilities(m1, m2, n, odds)
          probs[0..(k - lo)].inject(0.0, :+)
        end

        alias_method :exact_cdf, :cdf

        # Survival function: P(X > k)
        def sf(k, m1, m2, n, odds)
          validate_parameters(m1, m2, n, odds)
          lo, hi = support(m1, m2, n)
          return 1.0 if k < lo
          return 0.0 if k >= hi
          probs = probabilities(m1, m2, n, odds)
          probs[(k - lo + 1)..-1].inject(0.0, :+)
        end

        # Smallest k such that cdf(k) >= pr
        def quantile(pr, m1, m2, n, odds)
          validate_parameters(m1, m2, n, odds)
          lo, hi = support(m1, m2, n)
          probs = probabilities(m1, m2, n, odds)
          ac = 0.0
          (lo..hi).each do |k|
            ac += probs[k - lo]
            return k if ac >= pr
          end
          hi
        end

        alias_method :p_value, :quantile
        alias_method :exact_p_value, :quantile

        def mean(m1, m2, n, odds)
          lo, = support(m1, m2, n)
          probs = probabilities(m1, m2, n, odds)
          probs.each_with_index.inject(0.0) { |s, (p, i)| s + (lo + i) * p }
        end

        def variance(m1, m2, n, odds)
          lo, = support(m1, m2, n)
          probs = probabilities(m1, m2, n, odds)
          mu = probs.each_with_index.inject(0.0) { |s, (p, i)| s + (lo + i) * p }
          probs.each_with_index.inject(0.0) do |s, (p, i)|
            d = (lo + i) - mu
            s + d * d * p
          end
        end

        # Full probability vector over the support, computed from the mode
        # outwards with the recurrence
        #   w(k+1)/w(k) = odds * (m1-k)(n-k) / ((k+1)(m2-n+k+1))
        # so no binomial coefficient is ever evaluated as a Float and
        # relative weights only shrink away from the mode, keeping tails
        # accurate even for extremely large or small odds.
        def probabilities(m1, m2, n, odds)
          validate_parameters(m1, m2, n, odds)
          lo, hi = support(m1, m2, n)
          return [1.0] if lo == hi

          mode = find_mode(lo, hi, m1, m2, n, odds)
          probs = Array.new(hi - lo + 1, 0.0)
          probs[mode - lo] = 1.0

          k = mode
          while k < hi
            probs[k - lo + 1] = probs[k - lo] * recurrence_ratio(k, m1, m2, n, odds)
            k += 1
          end

          k = mode
          while k > lo
            probs[k - lo - 1] = probs[k - lo] / recurrence_ratio(k - 1, m1, m2, n, odds)
            k -= 1
          end

          sum = probs.inject(0.0, :+)
          probs.map { |p| p / sum }
        end

        private

        def validate_parameters(m1, m2, n, odds)
          fail(ArgumentError, 'm1 must be a non-negative integer') unless m1.is_a?(Integer) && m1 >= 0
          fail(ArgumentError, 'm2 must be a non-negative integer') unless m2.is_a?(Integer) && m2 >= 0
          fail(ArgumentError, 'n must be an integer with 0 <= n <= m1+m2') unless n.is_a?(Integer) && n >= 0 && n <= m1 + m2
          fail(ArgumentError, 'odds must be non-negative') unless odds.is_a?(Numeric) && odds >= 0
        end

        # w(k+1)/w(k); strictly decreasing in k over the support
        def recurrence_ratio(k, m1, m2, n, odds)
          odds * (m1 - k) * (n - k) / ((k + 1).to_f * (m2 - n + k + 1))
        end

        # Any mode of the distribution: the smallest k with ratio(k) <= 1,
        # found by binary search since the ratio decreases monotonically.
        def find_mode(lo, hi, m1, m2, n, odds)
          return lo if recurrence_ratio(lo, m1, m2, n, odds) <= 1
          return hi if recurrence_ratio(hi - 1, m1, m2, n, odds) >= 1
          a = lo
          b = hi - 1
          while a < b
            mid = (a + b) / 2
            if recurrence_ratio(mid, m1, m2, n, odds) <= 1
              b = mid
            else
              a = mid + 1
            end
          end
          a
        end
      end
    end
  end
end
