module Distribution
  # Dirichlet-multinomial distribution: a multinomial distribution over
  # +k+ categories whose event probabilities follow a Dirichlet
  # distribution with concentration vector +alphas+. It models
  # overdispersed counts for several categories observed in a single
  # batch of +n+ trials.
  #
  # == Usage
  #   dm = Distribution::DirichletMultinomial.new(10, [1.0, 2.0, 3.0])
  #   dm.log_pmf([2, 3, 5])        # joint log pmf, counts sum to n
  #   dm.marginal(0)               # beta-binomial margin of category 0
  #   dm.conditional(0 => 4)       # distribution of the remaining counts
  #   dm.merge_categories([0, 1], [2])
  #   dm.sample                    # vector of k counts summing exactly to n
  class DirichletMultinomial
    attr_reader :n, :alphas

    # +n+: total number of trials (non-negative Integer).
    # +alphas+: concentration vector, at least two positive numbers.
    def initialize(n, alphas)
      unless n.is_a?(Integer) && n >= 0
        fail ArgumentError, 'n must be a non-negative Integer'
      end
      unless alphas.is_a?(Array) && alphas.size >= 2
        fail ArgumentError, 'alphas must be an Array with at least two categories'
      end
      unless alphas.all? { |a| a.is_a?(Numeric) && a > 0 }
        fail ArgumentError, 'concentrations must be positive numbers'
      end

      @n = n
      @alphas = alphas.map(&:to_f).freeze
    end

    # Number of categories.
    def k
      @alphas.size
    end

    # Total concentration (sum of the concentration vector).
    def total_concentration
      @total_concentration ||= @alphas.inject(:+)
    end

    # Joint log pmf of the count vector +counts+. Returns
    # <tt>-Float::INFINITY</tt> for counts outside the support
    # (negative entries or entries not summing to +n+).
    def log_pmf(counts)
      counts = validate_counts(counts)
      return -Float::INFINITY if counts.nil?

      log_normalizer + Math.lgamma(@n + 1).first +
        counts.each_with_index.inject(0.0) do |ac, (x, i)|
          ac + Math.lgamma(x + @alphas[i]).first -
            Math.lgamma(@alphas[i]).first - Math.lgamma(x + 1).first
        end
    end

    # Joint pmf of the count vector +counts+.
    def pmf(counts)
      Math.exp(log_pmf(counts))
    end

    # Marginal distribution of the categories at +indices+, pooling the
    # remaining categories into a single "other" category. With a single
    # index this is the beta-binomial margin of that category.
    def marginal(*indices)
      indices = normalize_indices(indices.flatten)
      fail ArgumentError, 'at least one category must be left out' if indices.size >= k

      pooled = total_concentration - indices.inject(0.0) { |ac, i| ac + @alphas[i] }
      self.class.new(@n, indices.map { |i| @alphas[i] } + [pooled])
    end

    # Conditional distribution of the unobserved categories, given the
    # counts of some categories. +known+ is a Hash of
    # <tt>index => count</tt> or an Array with +nil+ for unknown
    # categories. Returns a new DirichletMultinomial over the remaining
    # categories, with its own reduced total and concentration vector.
    def conditional(known)
      known = normalize_known(known)
      remaining = (0...k).to_a - known.keys
      fail ArgumentError, 'cannot condition on every category' if remaining.empty?

      used = known.values.inject(0, :+)
      fail ArgumentError, 'known counts exceed n' if used > @n

      self.class.new(@n - used, remaining.map { |i| @alphas[i] })
    end

    # Distribution of counts pooled into +groups+ (arrays of category
    # indices partitioning <tt>0...k</tt>). Merging categories of a
    # Dirichlet-multinomial yields the same law as sampling first and
    # pooling the counts afterwards.
    def merge_categories(*groups)
      groups = groups.first if groups.size == 1 && groups.first.is_a?(Array) && groups.first.first.is_a?(Array)
      flat = groups.flatten
      unless flat.sort == (0...k).to_a
        fail ArgumentError, 'groups must partition the category indices'
      end

      self.class.new(@n, groups.map { |g| g.inject(0.0) { |ac, i| ac + @alphas[i] } })
    end

    # Draws a vector of +k+ counts summing exactly to +n+, using the
    # Pólya urn scheme (exact, no multinomial rounding). Pass a seeded
    # Random for reproducible draws.
    def sample(rng = Random)
      counts = Array.new(k, 0)
      weights = @alphas.dup
      total = total_concentration
      @n.times do
        r = rng.rand * total
        pick = k - 1
        cumulative = 0.0
        weights.each_with_index do |w, i|
          cumulative += w
          if r < cumulative
            pick = i
            break
          end
        end
        counts[pick] += 1
        weights[pick] += 1.0
        total += 1.0
      end
      counts
    end

    private

    def log_normalizer
      @log_normalizer ||=
        Math.lgamma(total_concentration).first -
        Math.lgamma(@n + total_concentration).first
    end

    def validate_counts(counts)
      fail ArgumentError, "expected #{k} counts" unless counts.is_a?(Array) && counts.size == k
      return nil unless counts.all? { |x| x.is_a?(Integer) && x >= 0 }
      return nil unless counts.inject(0, :+) == @n

      counts
    end

    def normalize_indices(indices)
      indices.each do |i|
        unless i.is_a?(Integer) && i >= 0 && i < k
          fail ArgumentError, "invalid category index #{i.inspect}"
        end
      end
      fail ArgumentError, 'duplicated category indices' if indices.uniq!
      indices
    end

    def normalize_known(known)
      pairs =
        case known
        when Hash then known.to_a
        when Array
          known.each_with_index.reject { |c, _| c.nil? }.map { |c, i| [i, c] }
        else
          fail ArgumentError, 'known must be a Hash or an Array with nils'
        end
      pairs.each do |i, c|
        normalize_indices([i])
        unless c.is_a?(Integer) && c >= 0
          fail ArgumentError, "known count for category #{i} must be a non-negative Integer"
        end
      end
      pairs.to_h
    end
  end
end
