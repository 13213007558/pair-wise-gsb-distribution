module Distribution
  # Finite mixture of continuous component distributions and discrete
  # point masses (atoms).
  #
  # A mixture is built from
  # * continuous components: any object responding to +cdf(x)+ and
  #   +quantile(p)+ (plus +pdf(x)+ or +density(x)+ if Mixture#density
  #   is needed), each with a non-negative weight, and
  # * atoms: explicit point masses at given locations.
  #
  # Weights are normalized to sum to 1 by default; pass
  # <tt>normalize: false</tt> to keep them untouched (they should then
  # already sum to 1). Negative weights are rejected outright. Atoms
  # given at the same location are merged into a single atom before
  # anything else is computed.
  #
  # The quantile function follows the generalized inverse CDF
  # definition: Q(p) = inf { x : F(x) >= p }. Jumps caused by atoms and
  # plateaus of the CDF (gaps between disjoint supports) are resolved
  # exactly according to that definition instead of being smoothed over
  # by continuous root-finding.
  #
  # == Example
  #   mix = Distribution::Mixture.new(
  #     components: [[0.25, dist_a], [0.25, dist_b]],
  #     atoms: [[0.5, 0.0]]
  #   )
  #   mix.cdf(0)        # => 0.5
  #   mix.mass_at(0)    # => 0.5
  #   mix.quantile(0.2) # => 0.0
  class Mixture
    # Continuous components as normalized [weight, component] pairs.
    attr_reader :components
    # Merged atoms as [location, mass] pairs, sorted by location.
    attr_reader :atoms
    # Sum of all (normalized) weights.
    attr_reader :total_weight

    # Number of bisection steps used when a quantile falls inside a
    # region where several continuous components overlap.
    BISECTION_STEPS = 100
    # Maximum doublings used to bracket an unbounded quantile search.
    BRACKET_EXPANSION_STEPS = 64

    def initialize(components: [], atoms: [], normalize: true)
      component_pairs = components.map do |weight, component|
        Mixture.validate_weight!(weight)
        unless component.respond_to?(:cdf) && component.respond_to?(:quantile)
          fail ArgumentError, 'Components must respond to #cdf and #quantile'
        end
        [weight.to_f, component]
      end

      atom_entries = atoms.is_a?(Hash) ? atoms.map { |location, weight| [weight, location] } : atoms
      merged_atoms = Hash.new(0.0)
      atom_entries.each do |weight, location|
        Mixture.validate_weight!(weight)
        merged_atoms[location.to_f] += weight.to_f
      end

      if normalize
        total = component_pairs.sum(&:first) + merged_atoms.values.sum
        fail ArgumentError, 'Weights must sum to a positive total' unless total > 0
        component_pairs.each { |pair| pair[0] /= total }
        merged_atoms.each_key { |location| merged_atoms[location] /= total }
      end

      @components = component_pairs.reject { |weight, _| weight == 0 }
      @atoms = merged_atoms.reject { |_, mass| mass == 0 }
                           .map { |location, mass| [location, mass] }
                           .sort_by(&:first)
      @total_weight = @components.sum(&:first) + @atoms.sum { |_, mass| mass }
    end

    # Mixture probability density at x. Atoms carry no density; only the
    # continuous components contribute.
    def density(x)
      @components.sum do |weight, component|
        weight * Mixture.component_density(component, x)
      end
    end
    alias pdf density

    # Point mass located exactly at x (0.0 when x is not an atom).
    def mass_at(x)
      pair = @atoms.find { |location, _| location == x }
      pair ? pair[1] : 0.0
    end

    # Mixture cumulative distribution function at x: the weighted sum of
    # the component CDFs plus the mass of all atoms at or below x.
    def cdf(x)
      total = 0.0
      @components.each { |weight, component| total += weight * component.cdf(x) }
      @atoms.each { |location, mass| total += mass if location <= x }
      total
    end

    # Generalized inverse CDF: Q(p) = inf { x : F(x) >= p }.
    #
    # Atom jumps are resolved exactly (Q(p) is the atom location whenever
    # p falls inside the jump). Plateaus of F return the left edge of the
    # plateau. When p falls inside a region where a single continuous
    # component varies, the component's own quantile function is used;
    # only genuinely overlapping components are bracketed by bisection
    # that preserves the generalized-inverse semantics.
    def quantile(p)
      if p < 0 || p > 1
        fail RangeError, 'cdf value(p) must be from [0, 1]. '+        "Cannot find quantile for p=#{p}"
      end
      return -Float::INFINITY if p <= 0
      return right_endpoint if p >= @total_weight

      left = -Float::INFINITY
      @atoms.each do |location, mass|
        f_at = cdf(location)
        f_minus = f_at - mass
        return continuous_quantile(left, location, p) if p <= f_minus
        return location if p <= f_at
        left = location
      end
      continuous_quantile(left, Float::INFINITY, p)
    end
    alias p_value quantile

    # :nodoc:
    def self.validate_weight!(weight)
      unless weight.is_a?(Numeric)
        fail ArgumentError, "Weights must be numeric (got #{weight.inspect})"
      end
      if weight < 0
        fail ArgumentError, "Negative weights are not allowed (got #{weight})"
      end
    end

    # :nodoc:
    def self.component_density(component, x)
      if component.respond_to?(:pdf)
        component.pdf(x)
      elsif component.respond_to?(:density)
        component.density(x)
      else
        fail ArgumentError, 'Components must respond to #pdf or #density '+        'to compute the mixture density'
      end
    end

    private

    # Smallest x with F(x) >= total_weight: the largest atom location or
    # the largest component right endpoint, whichever is greater.
    def right_endpoint
      endpoints = @atoms.map(&:first)
      endpoints += @components.map { |_, component| component.quantile(1) }
      endpoints.max || -Float::INFINITY
    end

    # Component CDF extended to +-Infinity without calling user code.
    def component_cdf(component, x)
      case x
      when -Float::INFINITY then 0.0
      when Float::INFINITY then 1.0
      else component.cdf(x)
      end
    end

    # Quantile inside the open interval (left, right), where F is
    # continuous, F(left+) < p and F(right-) >= p.
    def continuous_quantile(left, right, p)
      active = @components.select do |_, component|
        component_cdf(component, left) < component_cdf(component, right)
      end
      return left if active.empty?

      if active.size == 1
        weight, component = active.first
        f_left = left == -Float::INFINITY ? 0.0 : cdf(left)
        offset = f_left - weight * component_cdf(component, left)
        q = (p - offset) / weight
        q = 0.0 if q < 0.0
        q = 1.0 if q > 1.0
        component.quantile(q)
      else
        bisect_quantile(left, right, p)
      end
    end

    # Bisection for inf { x : F(x) >= p } on a continuous stretch of F.
    # The invariant F(lo) < p <= F(hi) is kept, so the result converges
    # to the left edge of any plateau instead of an arbitrary point.
    def bisect_quantile(left, right, p)
      lo = left
      hi = right
      if lo == -Float::INFINITY
        lo = -1.0
        steps = 0
        while cdf(lo) >= p && steps < BRACKET_EXPANSION_STEPS
          lo *= 2
          steps += 1
        end
      end
      if hi == Float::INFINITY
        hi = 1.0
        steps = 0
        while cdf(hi) < p && steps < BRACKET_EXPANSION_STEPS
          hi *= 2
          steps += 1
        end
      end
      BISECTION_STEPS.times do
        mid = 0.5 * (lo + hi)
        break if mid == lo || mid == hi
        if cdf(mid) < p
          lo = mid
        else
          hi = mid
        end
      end
      hi
    end
  end
end
