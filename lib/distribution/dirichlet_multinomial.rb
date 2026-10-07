module Distribution
  # Dirichlet-multinomial (multivariate Polya) distribution.
  #
  # Counts are represented as an Integer array whose elements sum exactly to
  # the distribution's total count. Concentration values are copied so that
  # later changes to the caller's array cannot alter the distribution.
  class DirichletMultinomial
    attr_reader :n, :alpha, :groups

    def initialize(n, alpha, options = {})
      @n = validate_non_negative_integer(n, :n)
      @alpha = normalize_concentration(alpha)

      groups = options.fetch(:groups, nil) || (0...@alpha.length).to_a
      unless groups.is_a?(Array) && groups.length == @alpha.length
        fail ArgumentError, 'groups must have one label per category'
      end
      @groups = groups.map { |label| label.is_a?(Array) ? label.dup.freeze : label }.freeze

      random = options.fetch(:random, nil)
      @random = random.is_a?(Integer) ? Random.new(random) : random
    end

    def self.log_pmf(counts, n, alpha)
      new(n, alpha).log_pmf(counts)
    end

    def self.sample(n, alpha, rng = nil)
      new(n, alpha).sample(rng)
    end

    def total_count
      @n
    end

    def concentration
      @alpha
    end

    def length
      @alpha.length
    end

    alias_method :categories, :length

    def log_pmf(counts)
      return -Float::INFINITY unless valid_counts?(counts)

      total_alpha = @alpha.inject(0.0, :+)
      result = log_factorial(@n) + Math.lgamma(total_alpha).first -
               Math.lgamma(@n + total_alpha).first

      counts.each_with_index do |count, index|
        result += Math.lgamma(count + @alpha[index]).first -
                  Math.lgamma(@alpha[index]).first -
                  log_factorial(count)
      end
      result
    end

    alias_method :joint_log_pmf, :log_pmf
    alias_method :log_prob, :log_pmf

    def pmf(counts)
      probability = Math.exp(log_pmf(counts))
      probability.nan? ? 0.0 : probability
    end

    alias_method :prob, :pmf
    alias_method :pdf, :pmf
    alias_method :joint_pmf, :pmf

    def sample(rng = nil)
      generator = generator_for(rng)
      counts = Array.new(@alpha.length, 0)
      weights = @alpha.dup

      @n.times do
        index = categorical_sample(weights, generator)
        counts[index] += 1
        weights[index] += 1
      end

      counts
    end

    # Enumerates every count vector in the finite support. This is intended for
    # exact checks on small distributions rather than large Monte Carlo tests.
    def each_count(&block)
      unless block_given?
        return Enumerator.new do |yielder|
          enumerate_counts([], @n, @alpha.length, yielder)
        end
      end

      enumerate_counts([], @n, @alpha.length, block)
      nil
    end

    # Marginal obtained by aggregating categories into one or more groups.
    # Categories omitted from the grouping are aggregated into an :other group.
    #
    #   marginal(0)                 # BetaBinomial for category 0 versus rest
    #   marginal([0, 2])            # categories 0 and 2 merged versus rest
    #   marginal([0], [1, 2])       # two groups
    #   marginal(a: [0], b: [1, 2]) # labeled groups
    def marginal(*groups)
      labels, groups = extract_groups(groups)
      fail ArgumentError, 'at least one marginal group is required' if groups.empty?

      if groups.any?(&:empty?)
        fail ArgumentError, 'marginal groups must not be empty'
      end

      selected = groups.map { |group| group.map { |index| normalized_index(index) } }
      flattened = selected.flatten
      unless flattened.length == flattened.uniq.length
        fail ArgumentError, 'a category cannot appear in more than one group'
      end

      covered = flattened.sort
      unless covered.all? { |index| (0...@alpha.length).include?(index) }
        fail ArgumentError, 'marginal category is out of range'
      end

      remaining = (0...@alpha.length).to_a - covered
      unless remaining.empty?
        selected << remaining
        labels << :other
      end

      concentrations = selected.map { |group| group.inject(0.0) { |sum, index| sum + @alpha[index] } }

      if selected.length == 1
        DirichletMultinomial::BetaBinomial.new(
          @n, concentrations.first, 0.0, [labels.first], random: @random
        )
      elsif selected.length == 2 && remaining == selected.last
        DirichletMultinomial::BetaBinomial.new(
          @n, concentrations.first, concentrations.last, labels, random: @random
        )
      else
        self.class.new(@n, concentrations, groups: labels, random: @random)
      end
    end

    alias_method :marginal_distribution, :marginal

    def merge(*groups)
      marginal(*groups)
    end

    alias_method :merge_categories, :merge

    # Conditional distribution of the remaining categories after fixing some
    # observed counts. Fixed counts may be a Hash or an Array containing nil
    # for categories that remain random.
    def conditional(observed)
      fixed = normalize_observed(observed)
      observed_total = fixed.values.inject(0, :+)
      if observed_total > @n
        fail ArgumentError, 'observed counts exceed the total count'
      end

      remaining_indices = (0...@alpha.length).to_a - fixed.keys
      if remaining_indices.empty?
        fail ArgumentError, 'at least one category must remain random'
      end

      remaining_alpha = remaining_indices.map { |index| @alpha[index] }
      remaining_groups = remaining_indices.map { |index| @groups[index] }

      self.class.new(
        @n - observed_total,
        remaining_alpha,
        groups: remaining_groups,
        random: @random
      )
    end

    alias_method :conditional_distribution, :conditional
    alias_method :given, :conditional

    private

    def normalize_concentration(values)
      unless values.respond_to?(:to_ary)
        fail ArgumentError, 'alpha must be an array of positive concentrations'
      end

      converted = values.to_ary.map do |value|
        begin
          concentration = Float(value)
        rescue ArgumentError, TypeError
          raise ArgumentError, 'alpha must contain numeric values'
        end
        unless concentration.finite? && concentration > 0
          fail ArgumentError, 'alpha concentrations must be positive and finite'
        end
        concentration
      end

      fail ArgumentError, 'alpha must not be empty' if converted.empty?
      converted.freeze
    end

    def validate_non_negative_integer(value, name)
      unless value.is_a?(Integer) && value >= 0
        fail ArgumentError, "#{name} must be a non-negative Integer"
      end
      value
    end

    def valid_counts?(counts)
      return false unless counts.respond_to?(:to_ary)
      counts = counts.to_ary
      return false unless counts.length == @alpha.length
      return false unless counts.all? { |count| count.is_a?(Integer) && count >= 0 }
      counts.inject(0, :+) == @n
    end

    def log_factorial(value)
      Math.lgamma(value + 1).first
    end

    def generator_for(rng)
      generator = rng || @random || Random
      generator.is_a?(Integer) ? Random.new(generator) : generator
    end

    def categorical_sample(weights, generator)
      threshold = generator.rand(weights.inject(0.0, :+))
      cumulative = 0.0

      weights.each_with_index do |weight, index|
        cumulative += weight
        return index if threshold < cumulative
      end

      weights.length - 1
    end

    def enumerate_counts(prefix, remaining, categories, yielder)
      if categories == 1
        emit(yielder, prefix + [remaining])
        return
      end

      (0..remaining).each do |count|
        enumerate_counts(prefix + [count], remaining - count, categories - 1, yielder)
      end
    end

    def emit(yielder, counts)
      if yielder.is_a?(Proc)
        yielder.call(counts)
      else
        yielder << counts
      end
    end

    def extract_groups(groups)
      if groups.length == 1 && groups.first.is_a?(Hash)
        mapping = groups.first
        [mapping.keys.to_a, mapping.values.map { |group| Array(group) }]
      else
        labels = []
        mapped = groups.map do |group|
          values = Array(group)
          labels << (values.length == 1 ? values.first : values.freeze)
          values
        end
        [labels, mapped]
      end
    end

    def normalized_index(index)
      Integer(index)
    rescue ArgumentError, TypeError
      raise ArgumentError, 'category indices must be integers'
    end

    def normalize_observed(observed)
      case observed
      when Hash
        observed.each_with_object({}) do |(index, count), fixed|
          normalized = normalized_index(index)
          unless (0...@alpha.length).include?(normalized)
            fail ArgumentError, 'observed category is out of range'
          end
          fixed[normalized] = validate_observed_count(count)
        end
      when Array
        unless observed.length == @alpha.length
          fail ArgumentError, 'observed count array must have one entry per category'
        end

        observed.each_with_index.each_with_object({}) do |(count, index), fixed|
          fixed[index] = validate_observed_count(count) unless count.nil?
        end
      else
        fail ArgumentError, 'observed counts must be a Hash or Array'
      end
    end

    def validate_observed_count(count)
      unless count.is_a?(Integer) && count >= 0
        fail ArgumentError, 'observed counts must be non-negative Integers'
      end
      count
    end

    # Beta-binomial marginal of one Dirichlet-multinomial category aggregated
    # against all remaining categories. Degenerate zero concentration is used
    # for the omitted side when every category is selected.
    class BetaBinomial
      attr_reader :n, :alpha, :beta, :groups

      def initialize(n, alpha, beta, groups = nil, options = {})
        @n = validate_non_negative_integer(n, :n)
        @alpha = normalize_parameter(alpha, :alpha)
        @beta = normalize_parameter(beta, :beta)
        if @alpha.zero? && @beta.zero?
          fail ArgumentError, 'at least one concentration must be positive'
        end
        @groups = (groups.nil? || groups.empty? ? [0, :other] : groups).dup.freeze
        @random = options.fetch(:random, nil)
      end

      def total_count
        @n
      end

      def concentration
        [@alpha, @beta]
      end

      def log_pmf(count)
        return -Float::INFINITY unless count.is_a?(Integer) && count >= 0 && count <= @n
        return count.zero? ? 0.0 : -Float::INFINITY if @alpha.zero?
        return count == @n ? 0.0 : -Float::INFINITY if @beta.zero?

        failures = @n - count
        total = @alpha + @beta

        log_factorial(@n) -
          log_factorial(count) -
          log_factorial(failures) -
          Math.lgamma(@alpha).first -
          Math.lgamma(@beta).first +
          Math.lgamma(total).first +
          Math.lgamma(count + @alpha).first +
          Math.lgamma(failures + @beta).first -
          Math.lgamma(@n + total).first
      end

      alias_method :log_prob, :log_pmf

      def pmf(count)
        probability = Math.exp(log_pmf(count))
        probability.nan? ? 0.0 : probability
      end

      alias_method :prob, :pmf
      alias_method :pdf, :pmf

      def cdf(value)
        return 0.0 if value < 0
        return 1.0 if value >= @n
        upper = value.floor
        (0..upper).inject(0.0) { |sum, count| sum + pmf(count) }
      end

      def sample(rng = nil)
        return 0 if @alpha.zero?
        return @n if @beta.zero?

        generator = generator_for(rng)
        successes = 0
        success_weight = @alpha
        failure_weight = @beta

        @n.times do
          threshold = generator.rand(success_weight + failure_weight)
          if threshold < success_weight
            successes += 1
            success_weight += 1
          else
            failure_weight += 1
          end
        end

        successes
      end

      private

      def validate_non_negative_integer(value, name)
        unless value.is_a?(Integer) && value >= 0
          fail ArgumentError, "#{name} must be a non-negative Integer"
        end
        value
      end

      def normalize_parameter(value, name)
        parameter = Float(value)
        unless parameter.finite? && parameter >= 0
          fail ArgumentError, "#{name} must be a non-negative finite number"
        end
        parameter
      rescue ArgumentError, TypeError
        raise ArgumentError, "#{name} must be numeric"
      end

      def log_factorial(value)
        Math.lgamma(value + 1).first
      end

      def generator_for(rng)
        generator = rng || @random || Random
        generator.is_a?(Integer) ? Random.new(generator) : generator
      end
    end
  end
end
