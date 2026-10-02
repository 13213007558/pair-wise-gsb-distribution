module Distribution
  module PoissonBinomial
    module Ruby_
      class << self
        def pmf(k, probabilities)
          probabilities = validate_probabilities(probabilities)
          masses = probability_masses(probabilities)
          mass_at(masses, validate_count(k))
        end

        alias_method :pdf, :pmf

        def cdf(k, probabilities)
          probabilities = validate_probabilities(probabilities)
          masses = probability_masses(probabilities)
          k = validate_count(k).floor
          return zero_result(masses) if k < 0
          return one_result(masses) if k >= masses.length - 1

          sum_masses(masses[0..k])
        end

        def right_tail(k, probabilities)
          probabilities = validate_probabilities(probabilities)
          masses = probability_masses(probabilities)

          k = validate_count(k).floor
          return one_result(masses) if k <= 0
          return zero_result(masses) if k >= masses.length

          sum_masses(masses[k..-1])
        end

        alias_method :ccdf, :right_tail
        alias_method :survival, :right_tail

        def quantile(probability, probabilities)
          probabilities = validate_probabilities(probabilities)
          masses = probability_masses(probabilities)

          validate_quantile_probability(probability)
          return 0 if probability.zero?
          return masses.length - 1 if probability == 1

          cumulative = zero_value(masses)
          masses.each_with_index do |mass, k|
            cumulative += mass
            return k if probability <= cumulative
          end
          masses.length - 1
        end

        alias_method :p_value, :quantile
        alias_method :exact_pmf, :pmf
        alias_method :exact_pdf, :pmf
        alias_method :exact_cdf, :cdf
        alias_method :exact_right_tail, :right_tail
        alias_method :exact_ccdf, :right_tail
        alias_method :exact_quantile, :quantile
        alias_method :exact_p_value, :quantile

        private

        def validate_probabilities(probabilities)
          fail ArgumentError, 'probabilities must respond to #to_a' unless probabilities.respond_to?(:to_a)

          probabilities.to_a.map do |probability|
            unless probability.is_a?(Numeric)
              fail ArgumentError, "invalid probability #{probability.inspect}"
            end

            if probability.respond_to?(:nan?) && probability.nan?
              fail ArgumentError, 'probability must not be NaN'
            end

            unless probability >= 0 && probability <= 1
              fail ArgumentError, "probability #{probability} is outside [0, 1]"
            end

            probability
          end
        end

        def probability_masses(probabilities)
          certain = probabilities.count { |probability| probability == 1 }
          uncertain = probabilities.reject { |probability| probability == 0 || probability == 1 }
          initial_mass = probabilities.any? { |probability| probability.is_a?(Float) } ? 1.0 : 1
          masses = [initial_mass]

          uncertain.each do |probability|
            next_masses = Array.new(masses.length + 1, zero_value(masses, probability))
            masses.each_with_index do |mass, k|
              next_masses[k] += mass * (1 - probability)
              next_masses[k + 1] += mass * probability
            end
            masses = next_masses
          end

          masses.unshift(*Array.new(certain, zero_value(masses))) unless certain.zero?

          if masses.any? { |mass| mass.is_a?(Float) }
            normalize_float_masses(masses.map(&:to_f))
          else
            masses
          end
        end

        def normalize_float_masses(masses)
          total = masses.inject(0.0, :+)
          return masses if total.zero?

          masses.map! { |mass| mass / total }
          largest = masses.index(masses.max)
          masses[largest] += 1.0 - masses.inject(0.0, :+)
          masses
        end

        def zero_value(*values)
          values.any? { |value| value.is_a?(Float) } ? 0.0 : 0
        end

        def mass_at(masses, k)
          integer = k.respond_to?(:integer?) ? k.integer? : k.is_a?(Integer)
          return zero_result(masses) unless integer
          return zero_result(masses) if k < 0 || k >= masses.length

          masses[k]
        end

        def sum_masses(masses)
          if masses.any? { |mass| mass.is_a?(Float) }
            masses.map(&:to_f).inject(0.0, :+)
          else
            masses.inject(0, :+)
          end
        end

        def zero_result(masses)
          masses.any? { |mass| mass.is_a?(Float) } ? 0.0 : 0
        end

        def one_result(masses)
          masses.any? { |mass| mass.is_a?(Float) } ? 1.0 : 1
        end

        def validate_count(count)
          unless count.is_a?(Numeric)
            fail ArgumentError, "invalid count #{count.inspect}"
          end

          if count.respond_to?(:nan?) && count.nan?
            fail ArgumentError, 'count must not be NaN'
          end

          if count.respond_to?(:infinite?) && count.infinite?
            fail ArgumentError, 'count must be finite'
          end

          count
        end

        def validate_quantile_probability(probability)
          unless probability.is_a?(Numeric)
            fail ArgumentError, "invalid probability #{probability.inspect}"
          end

          if probability.respond_to?(:nan?) && probability.nan?
            fail ArgumentError, 'probability must not be NaN'
          end

          unless probability >= 0 && probability <= 1
            fail ArgumentError, "quantile probability #{probability} is outside [0, 1]"
          end
        end
      end
    end
  end
end
