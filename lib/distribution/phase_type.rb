module Distribution
  # A phase-type (PH) distribution is the distribution of the absorption
  # time of a finite Markov jump process with a transient generator
  # +T+ (subintensity matrix) and an initial probability row vector
  # +alpha+ over the transient states.
  #
  # The defect of +alpha+, <tt>1 - sum(alpha)</tt>, is the probability
  # mass at zero (immediate absorption) and is reported separately as
  # +zero_mass+; it is never silently redistributed over the states.
  #
  # The generator is validated on construction: it must be a square
  # subintensity matrix (negative diagonal, non-negative off-diagonal,
  # row sums <= 0) and it must be nonsingular, i.e. absorption happens
  # with probability one. A model that cannot absorb raises
  # NotAbsorbingError instead of being renormalized into a proper
  # lifetime distribution.
  #
  # Pure Ruby implementation; no external matrix library is required.
  # Matrix exponentials are evaluated by uniformization of the
  # generator, and moments by an LU factorization of -T.
  #
  # == Usage
  #   # Erlang-2 with rate 2
  #   erlang = Distribution::PhaseType.new([1, 0], [[-2, 2], [0, -2]])
  #   erlang.survival(0.5)   #=> P(T > 0.5) = 2 * exp(-1)
  #   erlang.pdf(0.5)
  #   erlang.moment(1)       #=> 1.0
  #   erlang.residual_lifetime(0.5) #=> PH distribution of T - 0.5 given T > 0.5
  class PhaseType
    # Raised when alpha or the generator do not form a valid
    # transient phase-type model.
    class InvalidModelError < ArgumentError; end

    # Raised when the generator is singular, i.e. the process has a
    # closed class of transient states and never absorbs with
    # probability one. Such a model does not define a proper lifetime
    # and is refused instead of being normalized.
    class NotAbsorbingError < InvalidModelError; end

    TOLERANCE = 1e-10

    # Initial probability row vector over the transient states.
    attr_reader :alpha
    # Transient generator (subintensity matrix).
    attr_reader :generator
    # Probability mass at zero: 1 - sum(alpha).
    attr_reader :zero_mass
    # Number of transient states.
    attr_reader :dimension

    def initialize(alpha, generator)
      @generator = validate_generator(generator)
      @alpha = validate_alpha(alpha, @generator.size)
      @dimension = @generator.size
      @zero_mass = 1.0 - @alpha.inject(0.0, :+)
      @exit_vector = @generator.map { |row| -row.inject(0.0, :+) }
      @uniformization_rate = (0...@dimension).map { |i| @generator[i][i].abs }.max
      @lu, @lu_permutation = lu_factor(@generator.map { |row| row.map { |v| -v } })
    end

    # P(T > t); survival(0) is 1 - zero_mass.
    def survival(t)
      raise ArgumentError, 't must be numeric' unless t.is_a? Numeric
      return 1.0 if t < 0
      return 1.0 - @zero_mass if t.zero?
      dot(alpha_exp(t), ones)
    end

    # Density of the continuous part at t. Integrates to
    # <tt>1 - zero_mass</tt>, not to one.
    def pdf(t)
      raise ArgumentError, 't must be numeric' unless t.is_a? Numeric
      return 0.0 if t <= 0
      dot(alpha_exp(t), @exit_vector)
    end

    # P(T <= t); cdf(0) equals zero_mass.
    def cdf(t)
      raise ArgumentError, 't must be numeric' unless t.is_a? Numeric
      return 0.0 if t < 0
      return @zero_mass if t.zero?
      1.0 - survival(t)
    end

    # Raw k-th moment, E[T^k] = k! * alpha * (-T)^(-k) * 1.
    # The mass at zero contributes nothing for k >= 1.
    def moment(k)
      raise ArgumentError, 'k must be a non-negative integer' unless k.is_a?(Integer) && k >= 0
      return 1.0 if k.zero?
      vector = ones
      k.times { vector = lu_solve(vector) }
      factorial(k) * dot(@alpha, vector)
    end

    def mean
      moment(1)
    end

    def variance
      m = mean
      moment(2) - m * m
    end

    # Distribution of the residual lifetime T - t given T > t:
    # a phase-type distribution with the same generator and the
    # renormalized state probabilities alpha * exp(T t) / P(T > t).
    def residual_lifetime(t)
      raise ArgumentError, 't must be non-negative' unless t.is_a?(Numeric) && t >= 0
      survived = survival(t)
      raise InvalidModelError, "survival past #{t} has probability zero" if survived <= 0
      conditional = alpha_exp(t).map { |p| p / survived }
      PhaseType.new(conditional, @generator)
    end

    private

    def ones
      Array.new(@dimension, 1.0)
    end

    def dot(u, v)
      u.each_with_index.inject(0.0) { |sum, (x, i)| sum + x * v[i] }
    end

    def factorial(k)
      (1..k).inject(1, :*)
    end

    # alpha * exp(T t) by uniformization:
    # exp(T t) = sum_k exp(-q t) (q t)^k / k! * (I + T/q)^k,
    # with q = max |t_ii|. All iterates stay non-negative, so the
    # series is numerically stable for any subintensity matrix.
    def alpha_exp(t)
      return @alpha.dup if t.zero?
      q = @uniformization_rate
      x = q * t
      transition = (0...@dimension).map do |i|
        (0...@dimension).map do |j|
          @generator[i][j] / q + (i == j ? 1.0 : 0.0)
        end
      end
      vector = @alpha.dup
      result = Array.new(@dimension, 0.0)
      weight = Math.exp(-x)
      k = 0
      loop do
        @dimension.times { |j| result[j] += weight * vector[j] }
        k += 1
        vector = vector_matrix_product(vector, transition)
        weight *= x / k
        break if k > x && weight < 1e-17
      end
      result
    end

    def vector_matrix_product(vector, matrix)
      n = @dimension
      (0...n).map do |j|
        sum = 0.0
        n.times { |i| sum += vector[i] * matrix[i][j] }
        sum
      end
    end

    # LU factorization with partial pivoting of the matrix +a+.
    # Returns [lu, permutation]; a zero pivot means the transient
    # generator is singular and absorption is not certain.
    def lu_factor(a)
      n = a.size
      lu = a.map(&:dup)
      permutation = (0...n).to_a
      scale = lu.flatten.map(&:abs).max
      (0...n).each do |k|
        pivot = (k...n).max_by { |i| lu[i][k].abs }
        if lu[pivot][k].abs <= TOLERANCE * scale
          raise NotAbsorbingError,
                'generator is singular: the process never absorbs with ' +                'probability one and does not define a proper lifetime'
        end
        if pivot != k
          lu[k], lu[pivot] = lu[pivot], lu[k]
          permutation[k], permutation[pivot] = permutation[pivot], permutation[k]
        end
        (k + 1...n).each do |i|
          lu[i][k] /= lu[k][k]
          (k + 1...n).each { |j| lu[i][j] -= lu[i][k] * lu[k][j] }
        end
      end
      [lu, permutation]
    end

    # Solve (-T) x = b using the cached LU factorization.
    def lu_solve(b)
      n = @dimension
      x = (0...n).map { |i| b[@lu_permutation[i]].to_f }
      (0...n).each do |i|
        (0...i).each { |j| x[i] -= @lu[i][j] * x[j] }
      end
      (n - 1).downto(0) do |i|
        (i + 1...n).each { |j| x[i] -= @lu[i][j] * x[j] }
        x[i] /= @lu[i][i]
      end
      x
    end

    def validate_generator(generator)
      unless generator.is_a?(Array) && !generator.empty? &&
             generator.all? { |row| row.is_a?(Array) && row.size == generator.size }
        raise InvalidModelError, 'generator must be a non-empty square matrix'
      end
      generator.each { |row| row.each { |value| numeric_entry(value) } }
      scale = generator.flatten.map(&:abs).max
      scale = 1.0 if scale.zero?
      tol = TOLERANCE * scale
      generator.each_with_index do |row, i|
        row.each_with_index do |value, j|
          if i == j
            if value >= 0
              raise InvalidModelError, "generator diagonal entry (#{i}, #{i}) must be negative"
            end
          elsif value < -tol
            raise InvalidModelError, "generator off-diagonal entry (#{i}, #{j}) must be non-negative"
          end
        end
        row_sum = row.inject(0.0, :+)
        if row_sum > tol
          raise InvalidModelError,
                "generator row #{i} sums to #{row_sum}: rows of a transient generator must sum to <= 0"
        end
      end
      generator.map { |row| row.map(&:to_f) }
    end

    def validate_alpha(alpha, n)
      unless alpha.is_a?(Array) && alpha.size == n
        raise InvalidModelError, "alpha must be a vector of length #{n}"
      end
      alpha.each do |value|
        numeric_entry(value)
        if value < -TOLERANCE
          raise InvalidModelError, 'alpha entries must be non-negative'
        end
      end
      total = alpha.inject(0.0, :+)
      if total > 1.0 + TOLERANCE
        raise InvalidModelError, "alpha sums to #{total}: an initial probability vector cannot exceed one"
      end
      alpha.map(&:to_f)
    end

    def numeric_entry(value)
      raise InvalidModelError, 'alpha and generator entries must be numeric' unless value.is_a? Numeric
      value
    end
  end
end
