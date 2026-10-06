module Distribution
  # Continuous-time phase-type distribution.
  #
  # A phase-type random variable is the absorption time of a continuous-time
  # Markov chain started with transient-state probabilities +alpha+. The
  # sub-generator +t+ has non-negative off-diagonal transition rates and
  # non-positive row sums. An initial atom at absorption may be supplied as
  # +initial_absorbing+; states that cannot reach absorption contribute mass at
  # infinity.
  module PhaseType
    DEFAULT_TOLERANCE = 1e-12
    TAYLOR_TARGET = 0.5
    TAYLOR_MAX_ITERATIONS = 1000

    def self.new(alpha, t, initial_absorbing = nil, tolerance: DEFAULT_TOLERANCE)
      Ruby_.new(alpha, t, initial_absorbing, tolerance: tolerance)
    end

    class Ruby_
      attr_reader :alpha, :t, :initial_absorbing_probability,
                  :absorption_probability, :infinite_probability,
                  :reachable_states, :closed_states, :tolerance, :solver_info

      def initialize(alpha, t, initial_absorbing = nil, tolerance: DEFAULT_TOLERANCE)
        @tolerance = validate_finite_number(tolerance, 'tolerance').to_f
        raise ArgumentError, 'tolerance must be positive' unless @tolerance > 0

        @t = validate_matrix(t)
        @n = @t.size
        @alpha = validate_alpha(alpha, initial_absorbing)

        classify_states
        build_reachable_matrices
        solve_absorption_probabilities
      end

      def pdf(x, details: false, tolerance: @tolerance)
        ensure_numeric_time(x)
        return value_with_details(0.0, finite_information(:pdf, x, tolerance), details) if negative_time?(x)

        if infinite_time?(x)
          info = base_information(:pdf, Float::INFINITY).merge(error_bound: 0.0)
          return value_with_details(0.0, info, details)
        end

        time = finite_time(x)
        state, exp_info = state_distribution(time, tolerance)
        value = dot(state, @exit_rates)
        value = 0.0 if value < 0 && value.abs <= tolerance
        info = exp_info.merge(
          method: :pdf,
          time: time,
          error_bound: exp_info[:error_bound] * @exit_rates.sum.to_f
        )
        value_with_details(value, info, details)
      end

      def cdf(x, details: false, tolerance: @tolerance)
        ensure_numeric_time(x)
        return value_with_details(0.0, finite_information(:cdf, x, tolerance), details) if negative_time?(x)

        if infinite_time?(x)
          info = base_information(:cdf, Float::INFINITY).merge(error_bound: 0.0)
          return value_with_details(@absorption_probability, info, details)
        end

        time = finite_time(x)
        state, info = distribution_from_states(time, tolerance, :cdf)
        transient_mass = state.sum
        value = 1.0 - transient_mass
        value = [[value, @initial_absorbing_probability].max, @absorption_probability].min
        value_with_details(value, info, details)
      end

      def sf(x, details: false, tolerance: @tolerance)
        ensure_numeric_time(x)
        return value_with_details(1.0, finite_information(:sf, x, tolerance), details) if negative_time?(x)

        if infinite_time?(x)
          info = base_information(:sf, Float::INFINITY).merge(error_bound: 0.0)
          return value_with_details(@infinite_probability, info, details)
        end

        time = finite_time(x)
        state, info = distribution_from_states(time, tolerance, :sf)
        value = state.sum
        value = [[value, @infinite_probability].max, 1.0 - @initial_absorbing_probability].min
        value_with_details(value, info, details)
      end

      def quantile(p, details: false, tolerance: 1e-10, max_iterations: 100)
        probability = validate_probability(p)
        return quantile_result(0.0, quantile_info(:zero, probability, 0.0, 0.0, 0, tolerance), details) if probability <= @initial_absorbing_probability

        if probability >= @absorption_probability
          info = quantile_info(:beyond_absorption_probability, probability, nil, nil, 0, tolerance)
          return quantile_result(Float::INFINITY, info, details)
        end

        low = 0.0
        high = [1.0, conditional_mean * 2.0].max
        iterations = 0

        function_tolerance = Float::EPSILON / 10.0

        while cdf(high, tolerance: function_tolerance) < probability && high.infinite? == false && iterations < 1024
          high *= 2.0
          iterations += 1
        end

        if cdf(high, tolerance: function_tolerance) < probability
          info = quantile_info(:no_finite_bracket, probability, low, high, iterations, tolerance)
          return quantile_result(Float::INFINITY, info, details)
        end

        bisection_iterations = 0
        loop do
          middle = low + (high - low) / 2.0
          bisection_iterations += 1
          if cdf(middle, tolerance: function_tolerance) < probability
            low = middle
          else
            high = middle
          end

          error = (high - low) / [1.0, high.abs].max
          break if error <= tolerance || bisection_iterations >= max_iterations || high == low
        end

        status = bisection_iterations >= max_iterations ? :max_iterations : :converged
        info = quantile_info(status, probability, low, high, iterations + bisection_iterations, tolerance)
        quantile_result((low + high) / 2.0, info, details)
      end

      alias_method :p_value, :quantile

      def finite_moment(order, details: false)
        k = validate_order(order)
        value, solver = finite_raw_moment(k)
        info = base_information(:moment, nil).merge(
          order: k,
          finite_moment: value,
          solver_info: solver
        )
        value_with_details(value, info, details)
      end

      def moment(order, conditional: true, details: false)
        k = validate_order(order)
        finite = finite_moment(k)
        if conditional && @absorption_probability.zero? && !k.zero?
          raise ArgumentError, 'conditional moments are unavailable when absorption probability is zero'
        end

        value = if k.zero?
                  conditional ? 1.0 : @absorption_probability
                elsif conditional
                  finite / @absorption_probability
                elsif @infinite_probability > 0
                  Float::INFINITY
                else
                  finite
                end

        finite_info = finite_moment(k, details: true)
        info = base_information(:moment, nil).merge(
          order: k,
          conditional: conditional,
          finite_moment: finite,
          absorption_probability: @absorption_probability,
          solver_info: finite_info[:solver_info]
        )
        value_with_details(value, info, details)
      end

      def mean(conditional: true, details: false)
        moment(1, conditional: conditional, details: details)
      end

      def variance(conditional: true, details: false)
        if conditional == false && @infinite_probability > 0
          info = base_information(:moment, nil).merge(order: 2, conditional: false)
          return value_with_details(Float::INFINITY, info, details)
        end
        if conditional && @absorption_probability.zero?
          raise ArgumentError, 'conditional variance is unavailable when absorption probability is zero'
        end

        second = moment(2, conditional: conditional)
        first = moment(1, conditional: conditional)
        value = second - first * first
        info = base_information(:moment, nil).merge(order: 2, conditional: conditional, variance: value)
        value_with_details(value, info, details)
      end

      private

      def validate_finite_number(value, name)
        number = Float(value)
        raise ArgumentError, "#{name} must be finite" unless number.finite?
        number
      rescue ArgumentError, TypeError
        raise ArgumentError, "#{name} must be a number"
      end

      def validate_matrix(matrix)
        rows = matrix.to_a
        raise ArgumentError, 'sub-generator must be non-empty' if rows.empty?
        rows.map.with_index do |row, i|
          values = row.to_a
          if values.size != rows.size
            raise ArgumentError, "sub-generator must be square (row #{i} has #{values.size} entries)"
          end
          values.map { |value| validate_finite_number(value, 'sub-generator entry') }
        end
      rescue NoMethodError
        raise ArgumentError, 'sub-generator must be a matrix represented by arrays'
      end

      def validate_alpha(alpha, initial_absorbing)
        values = alpha.to_a
        unless values.size == @n
          raise ArgumentError, "initial probability length (#{values.size}) must match matrix size (#{@n})"
        end

        probabilities = values.map do |value|
          number = validate_finite_number(value, 'initial probability')
          raise ArgumentError, 'initial probabilities must be non-negative' if number < -@tolerance
          number
        end

        probabilities.each_with_index do |value, i|
          probabilities[i] = 0.0 if value < 0 && value.abs <= @tolerance
        end

        transient_mass = probabilities.sum
        scale = [1.0, transient_mass].max
        epsilon = @tolerance * scale

        if initial_absorbing.nil?
          raise ArgumentError, 'initial transient probabilities cannot exceed 1' if transient_mass > 1.0 + epsilon
          transient_mass = 1.0 if transient_mass > 1.0
          @initial_absorbing_probability = 1.0 - transient_mass
        else
          atom = validate_finite_number(initial_absorbing, 'initial absorbing probability')
          raise ArgumentError, 'initial absorbing probability must be non-negative' if atom < -@tolerance
          atom = [atom, 0.0].max
          total = transient_mass + atom
          raise ArgumentError, 'initial probabilities must sum to 1' unless (total - 1.0).abs <= @tolerance * [1.0, total].max
          factor = total.zero? ? 1.0 : 1.0 / total
          probabilities.map! { |value| value * factor }
          @initial_absorbing_probability = atom * factor
        end

        probabilities
      rescue NoMethodError
        raise ArgumentError, 'initial probabilities must respond to to_a'
      end

      def classify_states
        @t.each_with_index do |row, i|
          scale = [1.0, row.map { |value| value.abs }.sum].max
          epsilon = @tolerance * scale

          row.each_with_index do |value, j|
            next if i == j
            if value < -epsilon
              raise ArgumentError, "sub-generator off-diagonal entry (#{i}, #{j}) must be non-negative"
            elsif value < 0
              row[j] = 0.0
            end
          end

          if row[i] > epsilon
            raise ArgumentError, "sub-generator diagonal entry (#{i}, #{i}) must be non-positive"
          elsif row[i] > 0
            row[i] = 0.0
          end

          row_sum = row.sum
          row_scale = [1.0, row.map { |value| value.abs }.sum].max
          if row_sum > @tolerance * row_scale
            raise ArgumentError, "sub-generator row #{i} sum must be non-positive"
          end
        end

        @exit_rates = @t.map { |row| -row.sum }

        absorbing = @exit_rates.map.with_index { |rate, i| rate > 0 ? i : nil }.compact
        reverse = Array.new(@n) { [] }
        @t.each_with_index do |row, i|
          row.each_with_index { |rate, j| reverse[j] << i if rate > 0 }
        end

        @reachable_states = []
        queued = Array.new(@n, false)
        queue = absorbing.dup
        absorbing.each { |i| queued[i] = true }

        until queue.empty?
          state = queue.shift
          @reachable_states << state
          reverse[state].each do |previous|
            next if queued[previous]
            queued[previous] = true
            queue << previous
          end
        end
        @reachable_states.sort!
        @closed_states = (0...@n).to_a - @reachable_states
      end

      def build_reachable_matrices
        @reachable_alpha = @reachable_states.map { |i| @alpha[i] }
        @reachable_exit_rates = @reachable_states.map { |i| @exit_rates[i] }
        @minus_t_reachable = @reachable_states.map do |i|
          @reachable_states.map { |j| -@t[i][j] }
        end
      end

      def solve_absorption_probabilities
        if @reachable_states.empty?
          @absorption_probability = 0.0
          @infinite_probability = 1.0
          @solver_info = {
            method: :lu_solve,
            status: :no_reachable_states,
            relative_residual: 0.0,
            pivot_ratio: 1.0,
            matrix_one_norm: 0.0,
            solution_one_norm: 0.0,
            right_hand_side_one_norm: 0.0,
            reachable_states: [],
            closed_states: @closed_states,
            reachable_initial_mass: 0.0,
            closed_initial_mass: @alpha.sum,
            eventual_absorption_by_state: []
          }
          return
        end

        eventual_from_reachable, solve_info = lu_solve(
          @minus_t_reachable, @reachable_exit_rates
        )
        eventual_from_reachable = eventual_from_reachable.map { |value| clamp_probability(value) }

        reachable_mass = @reachable_alpha.sum
        closed_initial_mass = @closed_states.sum { |i| @alpha[i] }
        transient_absorbed = dot(@reachable_alpha, eventual_from_reachable)
        @absorption_probability = clamp_probability(
          @initial_absorbing_probability + transient_absorbed
        )
        @infinite_probability = clamp_probability(1.0 - @absorption_probability)
        @solver_info = solve_info.merge(
          reachable_states: @reachable_states,
          closed_states: @closed_states,
          reachable_initial_mass: reachable_mass,
          closed_initial_mass: closed_initial_mass,
          eventual_absorption_by_state: eventual_from_reachable
        )
      end

      def lu_solve(matrix, vector)
        n = matrix.size
        a = matrix.map(&:dup)
        b = vector.dup
        permutation = Array.new(n) { |i| i }
        pivot_ratios = []

        n.times do |col|
          pivot = (col...n).max_by { |row| a[row][col].abs }
          pivot_value = a[pivot][col].abs
          column_scale = (col...n).map { |row| a[row][col].abs }.max.to_f
          pivot_ratios << (column_scale.zero? ? 0.0 : pivot_value / column_scale)

          if pivot_value <= Float::EPSILON * [1.0, one_norm_matrix(a)].max
            info = lu_information(matrix, vector, vector, Float::INFINITY, pivot_ratios.min, :singular_subgenerator)
            raise ArgumentError, 'reachable sub-generator is numerically singular'
          end

          a[col], a[pivot] = a[pivot], a[col]
          b[col], b[pivot] = b[pivot], b[col]
          permutation[col], permutation[pivot] = permutation[pivot], permutation[col]

          ((col + 1)...n).each do |row|
            factor = a[row][col] / a[col][col]
            a[row][col] = factor
            ((col + 1)...n).each do |j|
              a[row][j] -= factor * a[col][j]
            end
          end
        end

        x = forward_backward(a, b)
        2.times do
          residual = vector.each_with_index.map do |value, i|
            value - dot(matrix[i], x)
          end
          correction = forward_backward(a, residual.values_at(*permutation))
          x = x.zip(correction).map { |value, delta| value + delta }
        end

        residual = vector.each_with_index.map { |value, i| value - dot(matrix[i], x) }
        scale = [1.0, one_norm_matrix(matrix) * l1_norm(x)].max
        relative_residual = l1_norm(residual) / scale
        near_singular = relative_residual > 1e-10 || pivot_ratios.min < 1e-12
        status = near_singular ? :near_singular : :converged
        [x, lu_information(matrix, vector, x, relative_residual, pivot_ratios.min, status)]
      end

      def finite_raw_moment(order)
        return [@absorption_probability, @solver_info] if order.zero?
        return [0.0, @solver_info] if @reachable_states.empty?
        values, infos = raw_moments(order)
        [values[order], infos.last]
      end

      def forward_backward(lu, vector)
        n = lu.size
        y = Array.new(n, 0.0)
        n.times do |i|
          y[i] = vector[i] - (0...i).sum { |j| lu[i][j] * y[j] }
        end

        x = Array.new(n, 0.0)
        (n - 1).downto(0) do |i|
          x[i] = (y[i] - ((i + 1)...n).sum { |j| lu[i][j] * x[j] }) / lu[i][i]
        end
        x
      end

      def lu_information(matrix, vector, solution, relative_residual, pivot_ratio, status)
        {
          method: :lu_solve,
          status: status,
          relative_residual: relative_residual,
          pivot_ratio: pivot_ratio,
          matrix_one_norm: one_norm_matrix(matrix),
          solution_one_norm: l1_norm(solution),
          right_hand_side_one_norm: l1_norm(vector)
        }
      end

      def raw_moments(order)
        moments = { 0 => Array.new(@reachable_states.size, 1.0) }
        infos = []
        (1..order).each do |k|
          values, info = lu_solve(@minus_t_reachable, moments[k - 1])
          moments[k] = values.map { |value| value * k }
          infos << info
        end
        values = {}
        moments.each do |k, state_values|
          values[k] = dot(@reachable_alpha, state_values)
        end
        [values, infos]
      end

      def conditional_mean
        return 0.0 if @absorption_probability <= @initial_absorbing_probability
        moment(1, conditional: true)
      end

      def state_distribution(time, tolerance)
        return [initial_state_distribution, base_information(:matrix_exponential, time).merge(error_bound: 0.0, iterations: 0)] if time.zero?

        t_reachable = @reachable_states.map do |i|
          @reachable_states.map { |j| @t[i][j] }
        end
        norm = one_norm_matrix(t_reachable)
        scaled_time = norm * time
        if norm.zero? || scaled_time.zero? || @reachable_states.empty?
          return [initial_state_distribution, base_information(:matrix_exponential, time).merge(error_bound: 0.0, iterations: 0)]
        end

        scale = [0, Math.log2(scaled_time / TAYLOR_TARGET).ceil].max
        dt = time * Math.ldexp(1.0, -scale)
        scaled_matrix = t_reachable.map { |row| row.map { |value| value * dt } }
        matrix_tolerance = [tolerance / 100.0, 1e-15].min
        expm, truncation_error, iteration_count = scaled_matrix_exponential(
          scaled_matrix, matrix_tolerance
        )

        rounding_error = 0.0
        scale.times do
          expm, squaring_error = square_matrix(expm)
          rounding_error += squaring_error
        end
        reachable_state = vector_matrix(@reachable_alpha, expm)
        state = Array.new(@n, 0.0)
        @reachable_states.each_with_index { |state_index, i| state[state_index] = reachable_state[i] }
        @closed_states.each { |state_index| state[state_index] = @alpha[state_index] }
        error_bound = one_norm_matrix(expm) * truncation_error + rounding_error

        info = base_information(:matrix_exponential, time).merge(
          status: error_bound <= tolerance ? :converged : :max_iterations,
          error_bound: error_bound,
          iterations: [iteration_count, scale],
          scaling_steps: scale
        )
        [state, info]
      end

      def initial_state_distribution
        @alpha.dup
      end

      def scaled_matrix_exponential(matrix, tolerance)
        n = matrix.size
        exponential = identity_matrix(n)
        term = identity_matrix(n)
        k = 1

        loop do
          scaled_matrix = matrix.map { |row| row.map { |value| value / k } }
          term = matrix_product(term, scaled_matrix)
          exponential = matrix_add(exponential, term)
          error = one_norm_matrix(term)
          return [exponential, error, k] if error <= tolerance || k >= TAYLOR_MAX_ITERATIONS
          k += 1
        end
      end

      def identity_matrix(n)
        Array.new(n) do |i|
          Array.new(n) { |j| i == j ? 1.0 : 0.0 }
        end
      end

      def matrix_add(left, right)
        left.zip(right).map do |left_row, right_row|
          left_row.zip(right_row).map { |x, y| x + y }
        end
      end

      def matrix_product(left, right)
        left.map do |row|
          (0...right.size).map do |column|
            (0...right.size).sum do |k|
              row[k] * right[k][column]
            end
          end
        end
      end

      def square_matrix(matrix)
        product = matrix_product(matrix, matrix)
        rounding_error = matrix.size * Float::EPSILON * one_norm_matrix(matrix)**2
        [product, rounding_error]
      end

      def vector_matrix(vector, matrix)
        matrix.transpose.map { |column| dot(vector, column) }
      end

      def distribution_from_states(time, tolerance, method)
        state, info = state_distribution(time, tolerance)
        [state, info.merge(method: method, error_bound: info[:error_bound])]
      end

      def finite_information(method, x, tolerance)
        base_information(method, x.to_f).merge(status: :outside_support, error_bound: 0.0, tolerance: tolerance)
      end

      def base_information(method, x)
        {
          method: method,
          time: x,
          status: :converged,
          initial_absorbing_probability: @initial_absorbing_probability,
          absorption_probability: @absorption_probability,
          infinite_probability: @infinite_probability
        }
      end

      def quantile_info(status, probability, low, high, iterations, tolerance)
        {
          method: :quantile,
          status: status,
          probability: probability,
          lower_bracket: low,
          upper_bracket: high,
          iterations: iterations,
          tolerance: tolerance,
          absorption_probability: @absorption_probability,
          infinite_probability: @infinite_probability
        }
      end

      def quantile_result(value, info, details)
        value_with_details(value, info, details)
      end

      def value_with_details(value, info, details)
        details ? { value: value, **info } : value
      end

      def validate_order(order)
        k = Integer(order)
        raise ArgumentError, 'moment order must be non-negative' if k < 0
        k
      rescue ArgumentError, TypeError
        raise ArgumentError, 'moment order must be an integer'
      end

      def validate_probability(value)
        p = Float(value)
        raise ArgumentError, 'probability must be between 0 and 1' unless p >= 0 && p <= 1
        p
      end

      def negative_time?(x)
        Float(x) < 0
      rescue ArgumentError, TypeError
        raise ArgumentError, 'time must be numeric'
      end

      def ensure_numeric_time(x)
        raise ArgumentError, 'time must not be NaN' if x.respond_to?(:nan?) && x.nan?

        begin
          time = Float(x)
        rescue ArgumentError, TypeError
          raise ArgumentError, 'time must be numeric'
        end
        raise ArgumentError, 'time must not be NaN' if time.nan?
        time
      end

      def infinite_time?(x)
        Float(x).infinite? == 1
      rescue ArgumentError, TypeError
        raise ArgumentError, 'time must be numeric'
      end

      def finite_time(x)
        time = Float(x)
        raise ArgumentError, 'time must not be NaN' if time.nan?
        raise ArgumentError, 'time must be finite or Infinity' if time.infinite?
        time
      rescue ArgumentError, TypeError
        raise ArgumentError, 'time must be numeric'
      end

      def dot(left, right)
        sum = 0.0
        left.each_with_index { |value, i| sum += value * right[i] }
        sum
      end

      def l1_norm(vector)
        vector.sum { |value| value.abs }
      end

      def one_norm_matrix(matrix)
        return 0.0 if matrix.empty?
        n = matrix.size
        (0...n).map do |column|
          matrix.sum { |row| row[column].abs }
        end.max.to_f
      end

      def clamp_probability(value)
        return 0.0 if value < 0 && value.abs <= @tolerance
        return 1.0 if value > 1 && value - 1.0 <= @tolerance
        [[value, 0.0].max, 1.0].min
      end
    end
  end
end
