module Distribution
  module PhaseType
    module Ruby_
      class Distribution
        attr_reader :alpha, :t, :initial_absorbed_probability,
                    :eventual_absorption_probability,
                    :infinite_mass, :reachable_states, :closed_states,
                    :reduced_alpha, :reduced_t, :absorption_probabilities

        def initialize(alpha, t, initial_absorbed_probability: nil,
                       tolerance: 1e-10)
          @tolerance = tolerance
          @alpha = validate_vector(alpha, 'alpha')
          @t = validate_matrix(t)
          check_dimensions
          validate_alpha(initial_absorbed_probability)
          validate_generator
          build_reachable_system
          calculate_eventual_mass
        end

        def order
          @alpha.length
        end

        def pdf(x, with_info: false, target_error: 1e-12)
          validate_time(x)
          return wrap(0.0, simple_info(x), with_info) if x < 0 || x == 0

          result = evolve(x, target_error, :pdf)
          wrap(result[:value], result[:info], with_info)
        end

        def cdf(x, with_info: false, target_error: 1e-12)
          validate_time(x)
          return wrap(0.0, simple_info(x), with_info) if x < 0
          return wrap(@initial_absorbed_probability, simple_info(x), with_info) if x == 0

          result = evolve(x, target_error, :cdf)
          wrap(result[:value], result[:info], with_info)
        end

        def sf(x, with_info: false, target_error: 1e-12)
          validate_time(x)
          return wrap(1.0, simple_info(x), with_info) if x < 0
          return wrap(1.0 - @initial_absorbed_probability,
                      simple_info(x), with_info) if x == 0

          result = evolve(x, target_error, :sf)
          wrap(result[:value], result[:info], with_info)
        end

        def moment(n, conditional: false, with_info: false)
          n = validate_order(n)
          if n.zero?
            value = conditional ? 1.0 : @eventual_absorption_probability
            return wrap(value, { order: 0, conditional: conditional }, with_info)
          end
          transient_mass = @reduced_alpha.sum
          if transient_mass.zero?
            if conditional
              raise Error, 'conditional moments are undefined when absorption has zero probability'
            end
            info = { order: n, conditional: conditional,
                     solves: [], note: :no_finite_absorption_mass }
            return wrap(0.0, info, with_info)
          end

          vector = @reduced_alpha.dup
          moments = []
          infos = []
          n.times do |k|
            order = k + 1
            rhs = vector.map { |value| order * value }
            negative_t = @reduced_t.transpose.map { |row| row.map { |value| -value } }
            solved = solve(negative_t, rhs)
            vector = solved[:solution]
            moments << dot(vector, @absorption_probabilities)
            infos << solved[:info]
          end

          value = moments.last
          value /= @eventual_absorption_probability if conditional
          info = { order: n, conditional: conditional, solves: infos }
          wrap(value, info, with_info)
        end

        def finite_moment(n, with_info: false)
          moment(n, conditional: false, with_info: with_info)
        end

        def conditional_moment(n, with_info: false)
          moment(n, conditional: true, with_info: with_info)
        end

        def mean
          @eventual_absorption_probability < 1.0 ? Float::INFINITY : moment(1)
        end

        def quantile(p, with_info: false, tolerance: 1e-10,
                     max_iterations: 200, target_error: 1e-12)
          validate_probability(p)
          if p <= @initial_absorbed_probability
            return wrap(0.0, quantile_info(0, p, 0, tolerance, 0, true, :zero_time), with_info)
          end
          if p >= @eventual_absorption_probability
            info = quantile_info(Float::INFINITY, p, nil, tolerance, 0, true,
                                :at_or_beyond_eventual_absorption_probability)
            return wrap(Float::INFINITY, info, with_info)
          end

          lower = 0.0
          upper = 1.0
          iterations = 0
          while cdf(upper, target_error: target_error) < p
            upper *= 2.0
            if upper.infinite?
              info = quantile_info(upper, p, nil, tolerance, iterations, false,
                                  :unbracketed_quantile)
              return wrap(Float::INFINITY, info, with_info)
            end
          end

          converged = false
          max_iterations.times do
            iterations += 1
            middle = (lower + upper) / 2.0
            result = cdf(middle, with_info: true, target_error: target_error)
            value = result.first
            if value < p
              lower = middle
            else
              upper = middle
            end

            width = upper - lower
            if width <= tolerance || width <= tolerance * upper.abs
              converged = true
              break
            end
          end

          result = cdf(upper, with_info: true, target_error: target_error)
          residual = (result.first - p).abs
          info = quantile_info(upper, p, [lower, upper], tolerance, iterations,
                               converged, converged ? nil : :max_iterations,
                               residual: residual, cdf_info: result.last)
          if !converged && !with_info
            raise ConvergenceError, "Quantile did not converge after #{max_iterations} iterations"
          end
          wrap(upper, info, with_info)
        end

        alias_method :p_value, :quantile

        private

        def validate_vector(values, name)
          unless values.is_a?(Array) && !values.empty?
            raise Error, "#{name} must be a non-empty array"
          end
          values.map do |value|
            number = Float(value)
            raise Error, "#{name} must contain finite numbers" unless number.finite?
            number
          end
        rescue ArgumentError, TypeError
          raise Error, "#{name} must contain numeric values"
        end

        def validate_matrix(matrix)
          unless matrix.is_a?(Array) && !matrix.empty?
            raise Error, 'T must be a non-empty array'
          end
          matrix.each do |row|
            raise Error, 'T rows must be arrays' unless row.is_a?(Array)
          end
          size = matrix.length
          unless matrix.all? { |row| row.length == size }
            raise Error, 'T must be a square matrix'
          end
          matrix.map { |row| validate_vector(row, 'T') }
        end

        def check_dimensions
          return if @alpha.length == @t.length
          raise Error, 'alpha and T must have the same dimension'
        end

        def validate_alpha(given)
          @alpha.map! do |value|
            if value < -@tolerance || value > 1 + @tolerance
              raise Error, 'alpha entries must be between 0 and 1'
            end
            [[value, 0.0].max, 1.0].min
          end

          transient_sum = @alpha.sum
          if transient_sum > 1 + @tolerance
            raise Error, 'transient probabilities cannot sum above 1'
          end

          if given.nil?
            @initial_absorbed_probability = clamp01(1.0 - transient_sum)
          else
            initial = Float(given)
            unless initial.finite? && initial >= -@tolerance && initial <= 1 + @tolerance
              raise Error, 'initial absorbed probability must be between 0 and 1'
            end
            @initial_absorbed_probability = clamp01(initial)
            total = transient_sum + @initial_absorbed_probability
            if (total - 1.0).abs > @tolerance
              raise Error, 'alpha and initial absorbed probability must sum to 1'
            end
          end
        rescue ArgumentError, TypeError
          raise Error, 'initial absorbed probability must be numeric'
        end

        def validate_generator
          scale = [1.0, @t.map { |row| row.map(&:abs).max }.max].max
          tolerance = @tolerance * scale
          numerical_zero = Float::EPSILON * scale
          @t.each_with_index do |row, i|
            row.each_with_index do |value, j|
              if i != j && value < -tolerance
                raise Error, 'off-diagonal T entries must be non-negative'
              end
              if i == j && value > tolerance
                raise Error, 'diagonal T entries must be non-positive'
              end
              if i != j && value < 0
                row[j] = 0.0
              elsif i == j && value > 0
                row[j] = 0.0
              elsif value.abs <= numerical_zero
                row[j] = 0.0
              end
            end
            row_sum = row.sum
            if row_sum > tolerance
              raise Error, 'each T row sum must be non-positive'
            end
            row[i] -= row_sum if row_sum > 0
          end
        end

        def build_reachable_system
          @reachable_states = states_reaching_absorption.sort
          @closed_states = (0...order).to_a - @reachable_states
          @reduced_t = @reachable_states.map do |i|
            @reachable_states.map { |j| @t[i][j] }
          end
          @reduced_alpha = @reachable_states.map { |i| @alpha[i] }
          @full_exit_rates = @t.map { |row| [-row.sum, 0.0].max }
          @uniformization_rate = @full_exit_rates.max
          @reduced_exit_rates = @reachable_states.map do |i|
            [-@t[i].sum, 0.0].max
          end
        end

        def states_reaching_absorption
          scale = [1.0, @t.map { |row| row.map(&:abs).max }.max].max
          numerical_zero = Float::EPSILON * scale
          reverse = Array.new(order) { [] }
          exits = Array.new(order, false)
          @t.each_with_index do |row, i|
            row.each_with_index { |value, j| reverse[j] << i if value > 0 }
            exits[i] = true if row.sum < -numerical_zero
          end
          queue = exits.each_index.select { |i| exits[i] }
          reachable = []
          until queue.empty?
            state = queue.shift
            next if reachable.include?(state)
            reachable << state
            queue.concat(reverse[state])
          end
          reachable
        end

        def calculate_eventual_mass
          if @reachable_states.empty?
            @absorption_probabilities = []
            @eventual_absorption_probability = @initial_absorbed_probability
          else
            negative_t = @reduced_t.map { |row| row.map { |value| -value } }
            solved = solve(negative_t, @reduced_exit_rates, raise_on_failure: false)
            @absorption_probabilities = solved[:solution].map { |value| clamp01(value) }
            absorbed_mass = dot(@reduced_alpha, @absorption_probabilities)
            @eventual_absorption_probability =
              clamp01(@initial_absorbed_probability + absorbed_mass)
          end
          @infinite_mass = clamp01(1.0 - @eventual_absorption_probability)
        end

        def evolve(x, target_error, kind)
          transient_mass = @alpha.sum
          if x.infinite?
            value = case kind
                    when :pdf then 0.0
                    when :sf then @infinite_mass
                    when :cdf then @eventual_absorption_probability
                    end
            return { value: value, info: { time: x, target_error: target_error,
                                          error: 0.0, iterations: 0,
                                          method: :asymptotic } }
          end

          rate = @uniformization_rate || 0.0
          mean_jumps = rate * x
          if @reachable_states.empty? || rate.zero? || @alpha.sum.zero?
            value = case kind
                    when :cdf then @initial_absorbed_probability
                    when :sf then 1.0 - @initial_absorbed_probability
                    when :pdf then 0.0
                    end
            return { value: value, info: { time: x, target_error: target_error,
                                          error: 0.0, iterations: 0,
                                          method: :trivial } }
          end

          if mean_jumps <= 200_000
            uniformized_evolution(x, target_error, kind, rate, transient_mass)
          else
            matrix_exponential_evolution(x, target_error, kind, rate, transient_mass)
          end
        end

        def uniformized_evolution(x, target_error, kind, rate, reduced_mass)
          m = rate * x
          p = transition_probability_matrix(rate)
          k_min = poisson_k_min(m, target_error)
          vector = @alpha.dup
          log_weight = -m
          if k_min.positive?
            1.upto(k_min) do |k|
              log_weight += Math.log(m) - Math.log(k)
              vector = vector_matrix_product(vector, p)
            end
          end
          weight = Math.exp(log_weight)
          transient_survival = 0.0
          absorbed_density = 0.0
          iterations = k_min
          lower_error = k_min.zero? || m <= k_min ? 0.0 : poisson_lower_tail(m, k_min)

          loop do
            mass = vector.sum
            transient_survival += weight * mass
            absorbed_density += weight * dot(vector, @full_exit_rates)
            next_ratio = m.fdiv(iterations + 1)
            upper_error = if next_ratio < 1.0
                            weight * next_ratio / (1.0 - next_ratio)
                          else
                            1.0
                          end
            error = lower_error + upper_error
            break if error <= target_error
            if mass.zero? && upper_error <= target_error
              break
            end
            if iterations > 500_000
              raise ConvergenceError, 'Uniformized evolution did not converge'
            end

            vector = vector_matrix_product(vector, p)
            weight *= next_ratio
            iterations += 1
            lower_error = 0.0
          end

          value = case kind
                  when :sf then transient_survival
                  when :cdf then @initial_absorbed_probability + reduced_mass - transient_survival
                  when :pdf then absorbed_density
                  end
          value = clamp01(value) unless kind == :pdf
          info = { time: x, target_error: target_error,
                   error: lower_error + poisson_tail_bound(m, iterations),
                   iterations: iterations - k_min + 1, method: :uniformization }
          { value: value, info: info }
        end

        def matrix_exponential_evolution(x, target_error, kind, rate, reduced_mass)
          n = order
          scale = [0, Math.log2(rate * x)].max.ceil
          dt = x / 2.0**scale
          m = rate * dt
          jump = full_transition_probability_matrix(rate)
          step = Array.new(n) do |i|
            Array.new(n) { |j| i == j ? Math.exp(-m) : 0.0 }
          end
          term = identity(n)
          weight = 1.0
          included_probability = Math.exp(-m)
          k = 0
          loop do
            k += 1
            weight *= m / k
            poisson_weight = weight * Math.exp(-m)
            term = matrix_product(term, jump)
            included_probability += poisson_weight
            n.times do |i|
              n.times { |j| step[i][j] += poisson_weight * term[i][j] }
            end
            break if poisson_weight <= target_error || k > 100
          end

          evolution = step
          scale.times { evolution = matrix_product(evolution, evolution) }
          transient = vector_matrix_product(@alpha, evolution)
          value = case kind
                  when :sf then transient.sum
                  when :cdf then @initial_absorbed_probability + reduced_mass - transient.sum
                  when :pdf then dot(transient, @full_exit_rates)
                  end
          value = clamp01(value) unless kind == :pdf
          truncation_error = [1.0 - included_probability, 0.0].max
          info = { time: x, target_error: target_error,
                   error: truncation_error, iterations: scale + k,
                   method: :scaled_matrix_exponential }
          { value: value, info: info }
        end

        def transition_probability_matrix(rate)
          full_transition_probability_matrix(rate)
        end

        def full_transition_probability_matrix(rate)
          @t.each_with_index.map do |row, i|
            row.each_with_index.map do |value, j|
              (i == j ? 1.0 : 0.0) + value / rate
            end
          end
        end

        def poisson_k_min(m, target_error)
          return 0 if m <= 40
          k = (m - Math.sqrt(2.0 * m * Math.log(1.0 / target_error))).floor
          [k, 0].max
        end

        def poisson_lower_tail(m, k)
          deviation = m - k
          Math.exp(-(deviation**2) / (2.0 * m))
        end

        def poisson_tail_bound(m, k)
          return 1.0 - Math.exp(-m) if k.zero?
          ratio = m / (k + 1.0)
          return 1.0 if ratio >= 1.0
          poisson_mass(m, k) * ratio / (1.0 - ratio)
        end

        def poisson_mass(m, k)
          Math.exp(-m + k * Math.log(m) - Math.lgamma(k + 1).first)
        end

        def vector_matrix_product(vector, matrix)
          Array.new(matrix.first.length) do |j|
            sum = 0.0
            vector.each_with_index { |value, i| sum += value * matrix[i][j] }
            sum
          end
        end

        def matrix_product(left, right)
          Array.new(left.length) do |i|
            Array.new(right.first.length) do |j|
              sum = 0.0
              left.first.length.times { |k| sum += left[i][k] * right[k][j] }
              sum
            end
          end
        end

        def identity(size)
          Array.new(size) do |i|
            Array.new(size) { |j| i == j ? 1.0 : 0.0 }
          end
        end

        def solve(matrix, rhs, raise_on_failure: true)
          a = matrix.map(&:dup)
          b = rhs.dup
          n = a.length
          scale = [1.0, a.flatten.map(&:abs).max].max

          n.times do |col|
            pivot = (col...n).max_by { |row| a[row][col].abs }
            if a[pivot][col].abs < 1e-14 * scale
              raise ConvergenceError, 'Linear system is numerically singular'
            end
            a[col], a[pivot] = a[col], a[pivot] if pivot != col
            b[col], b[pivot] = b[col], b[pivot] if pivot != col

            pivot_value = a[col][col]
            (col + 1...n).each do |row|
              factor = -a[row][col] / pivot_value
              next if factor.zero?
              a[row][col] = 0.0
              (col + 1...n).each { |j| a[row][j] += factor * a[col][j] }
              b[row] += factor * b[col]
            end
          end

          solution = Array.new(n)
          (n - 1).downto(0) do |row|
            sum = b[row]
            (row + 1...n).each { |j| sum -= a[row][j] * solution[j] }
            solution[row] = sum / a[row][row]
          end

          residual = Array.new(n) do |i|
            dot(matrix[i], solution) - rhs[i]
          end
          residual_norm = Math.sqrt(residual.sum { |value| value**2 })
          rhs_norm = Math.sqrt(rhs.sum { |value| value**2 })
          relative_error = residual_norm / [1.0, rhs_norm].max
          if relative_error > 1e-10 && raise_on_failure
            raise ConvergenceError, "Moment solve residual too large: #{relative_error}"
          end
          { solution: solution, info: { relative_residual: relative_error } }
        end

        def dot(left, right)
          sum = 0.0
          left.each_with_index { |value, i| sum += value * right[i] }
          sum
        end

        def simple_info(x)
          { time: x, error: 0.0, iterations: 0, method: :direct }
        end

        def quantile_info(value, p, bracket, tolerance, iterations, converged,
                          reason, residual: nil, cdf_info: nil)
          { value: value, p: p, bracket: bracket, tolerance: tolerance,
            iterations: iterations, converged: converged, reason: reason,
            residual: residual, cdf_info: cdf_info }.compact
        end

        def validate_time(x)
          Float(x)
        rescue ArgumentError, TypeError
          raise Error, 'time must be numeric'
        end

        def validate_probability(p)
          value = Float(p)
          unless value.finite? && value >= 0 && value <= 1
            raise Error, 'probability must be between 0 and 1'
          end
          value
        rescue ArgumentError, TypeError
          raise Error, 'probability must be numeric'
        end

        def validate_order(n)
          value = Integer(n)
          raise Error, 'moment order must be non-negative' if value < 0
          value
        rescue ArgumentError, TypeError
          raise Error, 'moment order must be an integer'
        end

        def clamp01(value)
          [[value, 0.0].max, 1.0].min
        end

        def wrap(value, info, with_info)
          with_info ? [value, info] : value
        end
      end

      def self.new(alpha, t, options = {})
        Distribution.new(alpha, t, **options)
      end
    end
  end
end
