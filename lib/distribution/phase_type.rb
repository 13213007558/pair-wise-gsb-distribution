module Distribution
  class PhaseType
    attr_reader :initial_probabilities, :generator, :initial_absorption_mass, :tolerance

    def initialize(initial_probabilities, generator, tolerance = 1e-12)
      @initial_probabilities = initial_probabilities.map { |value| Float(value) }
      @generator = generator.map { |row| row.map { |value| Float(value) } }
      @tolerance = Float(tolerance)

      validate!
    end

    def order
      generator.size
    end

    def survival(time)
      validate_time!(time)
      transient_probability(state_probabilities(time))
    end

    alias_method :ccdf, :survival

    def density(time)
      validate_time!(time)
      state = state_probabilities(time)
      value = state.each_with_index.sum { |probability, state_index|
        probability * exit_rates[state_index]
      }
      [value, 0.0].max
    end

    alias_method :pdf, :density

    def cdf(time)
      validate_time!(time)
      value = 1.0 - survival(time)
      [[value, initial_absorption_mass].max, 1.0].min
    end

    def moment(order)
      order = Integer(order)
      raise ArgumentError, 'moment order must be positive' if order <= 0

      moments_upto(order).last
    end

    def moments_upto(order)
      order = Integer(order)
      raise ArgumentError, 'moment order must be positive' if order <= 0

      moments = []
      previous = Array.new(self.order, 1.0)
      subgenerator = matrix_scale(generator, -1.0)

      order.times do
        current = solve_linear_system(subgenerator, vector_scale(previous, moments.size + 1.0))
        moments << vector_dot(initial_probabilities, current)
        previous = current
      end

      moments
    end

    def mean
      moment(1)
    end

    def variance
      moment(2) - mean**2
    end

    def residual_life(time)
      survived = survival(time)
      if survived <= numeric_tolerance
        raise RangeError, "cannot condition on survival at #{time}: probability is #{survived}"
      end

      conditional_probabilities = state_probabilities(time).map { |probability|
        probability / survived
      }
      self.class.new(conditional_probabilities, generator, tolerance)
    end

    alias_method :residual_lifetime, :residual_life

    def state_probabilities(time)
      validate_time!(time)
      vector_matrix_product(initial_probabilities, matrix_exponential(generator, time))
    end

    private

    def validate!
      raise ArgumentError, 'initial probability vector is empty' if initial_probabilities.empty?
      raise ArgumentError, 'generator is empty' if generator.empty?
      unless generator.size == initial_probabilities.size
        raise ArgumentError, 'initial probability vector and generator must have the same order'
      end
      unless generator.all? { |row| row.is_a?(Array) && row.size == generator.size }
        raise ArgumentError, 'generator must be square'
      end
      unless initial_probabilities.all?(&:finite?) && generator.all? { |row| row.all?(&:finite?) }
        raise ArgumentError, 'probabilities and generator entries must be finite'
      end
      unless tolerance.finite? && tolerance > 0
        raise ArgumentError, 'tolerance must be positive'
      end

      largest_magnitude = generator.flatten.map(&:abs).max.to_f
      @numeric_tolerance = tolerance * [1.0, largest_magnitude, order.to_f].max
      validate_initial_probabilities!
      validate_generator!
      validate_absorption!
    end

    def validate_initial_probabilities!
      clean_probabilities = initial_probabilities.map do |probability|
        if probability < -numeric_tolerance
          raise ArgumentError, 'initial probabilities cannot be negative'
        end

        probability.abs <= numeric_tolerance ? 0.0 : probability
      end
      initial_probabilities.replace(clean_probabilities)

      transient_mass = initial_probabilities.sum
      if transient_mass > 1.0 + numeric_tolerance
        raise ArgumentError, 'initial transient probability cannot exceed one'
      end

      @initial_absorption_mass =
        if transient_mass < numeric_tolerance
          1.0
        else
          [[1.0 - transient_mass, 0.0].max, 1.0].min
        end
    end

    def validate_generator!
      generator.each_with_index do |row, row_index|
        unless row[row_index] <= numeric_tolerance
          raise ArgumentError, 'generator diagonal entries must be non-positive'
        end

        row.each_with_index do |value, column_index|
          if column_index != row_index && value < -numeric_tolerance
            raise ArgumentError, 'generator off-diagonal entries must be non-negative'
          end
        end

        if row.sum > numeric_tolerance
          raise ArgumentError, 'generator rows must sum to zero or a negative absorption rate'
        end
      end
    end

    def validate_absorption!
      can_absorb = Array.new(order, false)
      predecessors = Array.new(order) { [] }

      generator.each_with_index do |row, row_index|
        if -row.sum > numeric_tolerance
          can_absorb[row_index] = true
        end

        row.each_with_index do |rate, column_index|
          predecessors[column_index] << row_index if row_index != column_index && rate > numeric_tolerance
        end
      end

      pending = can_absorb.each_index.select { |index| can_absorb[index] }
      until pending.empty?
        state_index = pending.shift
        predecessors[state_index].each do |predecessor|
          next if can_absorb[predecessor]

          can_absorb[predecessor] = true
          pending << predecessor
        end
      end

      unless can_absorb.all?
        raise ArgumentError, 'every phase must be able to reach absorption'
      end
    end

    def validate_time!(time)
      time = Float(time)
      raise ArgumentError, 'time must be non-negative' unless time.finite? && time >= 0
      time
    end

    def transient_probability(state)
      probability = state.sum
      probability = 0.0 if probability.abs <= numeric_tolerance
      [[probability, 0.0].max, 1.0].min
    end

    def exit_rates
      @exit_rates ||= generator.map { |row| [-row.sum, 0.0].max }
    end

    def matrix_exponential(matrix, time = 1.0)
      scaled_time = time.to_f
      infinity_norm = matrix.map { |row| row.map(&:abs).sum }.max.to_f
      scaling = infinity_norm.zero? ? 0 : [0, Math.log2(infinity_norm * scaled_time / 0.5)].max.ceil
      scaled_matrix = matrix.map { |row| row.map { |value| value * scaled_time / 2**scaling } }

      result = pade_exponential(scaled_matrix)
      scaling.times { result = matrix_multiply(result, result) }
      result
    end

    def pade_exponential(matrix)
      coefficients = [
        1.0,
        0.5,
        3.0 / 26.0,
        5.0 / 312.0,
        5.0 / 3432.0,
        1.0 / 11440.0,
        1.0 / 308880.0,
        1.0 / 17297280.0
      ]
      identity = identity_matrix(matrix.size)
      power2 = matrix_multiply(matrix, matrix)
      power4 = matrix_multiply(power2, power2)
      power6 = matrix_multiply(power4, power2)

      numerator = matrix_multiply(
        matrix,
        matrix_add(
          matrix_scale(identity, coefficients[1]),
          matrix_add(
            matrix_scale(power2, coefficients[3]),
            matrix_add(
              matrix_scale(power4, coefficients[5]),
              matrix_scale(power6, coefficients[7])
            )
          )
        )
      )
      denominator = matrix_add(
        matrix_scale(identity, coefficients[0]),
        matrix_add(
          matrix_scale(power2, coefficients[2]),
          matrix_add(
            matrix_scale(power4, coefficients[4]),
            matrix_scale(power6, coefficients[6])
          )
        )
      )

      solve_linear_system(
        matrix_add(denominator, numerator, -1.0),
        matrix_add(denominator, numerator)
      )
    end

    def identity_matrix(size)
      Array.new(size) do |row|
        Array.new(size) { |column| row == column ? 1.0 : 0.0 }
      end
    end

    def matrix_multiply(left, right)
      Array.new(left.size) do |row|
        Array.new(right.first.size) do |column|
          left.each_index.sum { |index| left[row][index] * right[index][column] }
        end
      end
    end

    def matrix_add(left, right, right_multiplier = 1.0)
      Array.new(left.size) do |row|
        Array.new(left.first.size) do |column|
          left[row][column] + right_multiplier * right[row][column]
        end
      end
    end

    def matrix_scale(matrix, multiplier)
      matrix.map { |row| row.map { |value| value * multiplier } }
    end

    def vector_scale(vector, multiplier)
      vector.map { |value| value * multiplier }
    end

    def vector_matrix_product(vector, matrix)
      Array.new(matrix.first.size) do |column|
        vector.each_index.sum { |row| vector[row] * matrix[row][column] }
      end
    end

    def vector_dot(left, right)
      left.each_index.sum { |index| left[index] * right[index] }
    end

    def solve_linear_system(matrix, right_hand_side)
      left = matrix.map(&:dup)
      rhs = right_hand_side.map(&:dup)
      size = left.size
      solve_tolerance = Float::EPSILON * [1.0, left.flatten.map(&:abs).max].max * size

      size.times do |pivot_column|
        pivot_row = (pivot_column...size).max_by { |row| left[row][pivot_column].abs }
        if left[pivot_row][pivot_column].abs <= solve_tolerance
          raise ArgumentError, 'generator has a non-absorbing phase'
        end

        left[pivot_row], left[pivot_column] = left[pivot_column], left[pivot_row]
        rhs[pivot_row], rhs[pivot_column] = rhs[pivot_column], rhs[pivot_row]

        ((pivot_column + 1)...size).each do |row|
          factor = left[row][pivot_column] / left[pivot_column][pivot_column]
          (pivot_column...size).each { |column| left[row][column] -= factor * left[pivot_column][column] }
          if rhs.first.is_a?(Array)
            rhs[row].each_index { |column| rhs[row][column] -= factor * rhs[pivot_column][column] }
          else
            rhs[row] -= factor * rhs[pivot_column]
          end
        end
      end

      solution = Array.new(size) { Array.new(rhs.first.is_a?(Array) ? rhs.first.size : 1) { 0.0 } }
      size.times do |pivot_column|
        row = size - pivot_column - 1
        if rhs.first.is_a?(Array)
          rhs[row].each_with_index do |value, column|
            known = ((row + 1)...size).sum { |index| left[row][index] * solution[index][column] }
            solution[row][column] = (value - known) / left[row][row]
          end
        else
          known = ((row + 1)...size).sum { |index| left[row][index] * solution[index][0] }
          solution[row][0] = (rhs[row] - known) / left[row][row]
        end
      end

      rhs.first.is_a?(Array) ? solution : solution.map(&:first)
    end

    attr_reader :numeric_tolerance
  end
end
