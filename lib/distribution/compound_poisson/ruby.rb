module Distribution
  module CompoundPoisson
    module Ruby_
      class << self
        ProbabilityMass = Struct.new(:pmf, :cdf, :range_mass, :omitted_mass) do
          def to_h
            {
              pmf: pmf,
              cdf: cdf,
              range_mass: range_mass,
              omitted_mass: omitted_mass
            }
          end
        end

        def pdf(k, lambda_val, cost_pmf)
          pmf_values, = build_pmf(k, lambda_val, cost_pmf)
          return 0.0 unless k.is_a?(Integer) && k >= 0 && k < pmf_values.length

          pmf_values[k]
        end

        alias_method :pmf, :pdf

        def cdf(k, lambda_val, cost_pmf)
          return 0.0 unless k.is_a?(Integer) && k >= 0

          _, cumulative = build_pmf(k, lambda_val, cost_pmf)
          cumulative[k]
        end

        def distribution(lambda_val, cost_pmf, range_or_min, max_cost = nil)
          range = validate_range(range_or_min, max_cost)
          max = [range.max, 0].max
          pmf_values, cumulative = build_pmf(max, lambda_val, cost_pmf)

          selected_pmf = {}
          selected_cdf = {}
          range.each do |cost|
            selected_pmf[cost] = cost >= 0 ? pmf_values[cost] : 0.0
            selected_cdf[cost] = cost >= 0 ? cumulative[cost] : 0.0
          end

          range_mass = selected_pmf.values.reduce(0.0, :+)
          below = probability_below(range.min, cumulative)
          above = 1.0 - cumulative_probability(range.max, cumulative)
          above = 0.0 if above < 0 && above > -Float::EPSILON

          omitted = {
            below: below,
            above: above,
            total: below + above
          }

          ProbabilityMass.new(selected_pmf, selected_cdf, range_mass, omitted)
        end

        alias_method :pmf_range, :distribution

        private

        def build_pmf(max_cost, lambda_val, cost_pmf)
          probabilities = normalize_probabilities(cost_pmf)
          validate_lambda(lambda_val)
          max_cost = 0 unless max_cost.is_a?(Integer) && max_cost >= 0

          zero_probability = probabilities.fetch(0, 0.0)
          positive_probabilities = probabilities.reject { |cost, _| cost.zero? }
          pmf_values = Array.new(max_cost + 1, 0.0)
          cumulative = Array.new(max_cost + 1, 0.0)
          pmf_values[0] = Math.exp(-lambda_val * (1.0 - zero_probability))

          running_mass = pmf_values[0]
          cumulative[0] = running_mass

          (1..max_cost).each do |cost|
            mass = positive_probabilities.reduce(0.0) do |sum, (single_cost, probability)|
              next sum if single_cost > cost

              sum + single_cost * probability * pmf_values[cost - single_cost]
            end
            pmf_values[cost] = lambda_val * mass / cost
            running_mass += pmf_values[cost]
            running_mass = 1.0 if running_mass > 1.0
            cumulative[cost] = running_mass
          end

          [pmf_values, cumulative]
        end

        def normalize_probabilities(cost_pmf)
          fail ArgumentError, 'Cost probability table is empty' if cost_pmf.nil? || cost_pmf.empty?

          probabilities = {}
          entries = cost_pmf.is_a?(Array) ? cost_pmf.each_with_index.map { |probability, cost| [cost, probability] } : cost_pmf
          entries.each do |cost, probability|
            begin
              integer_cost = Integer(cost)
              value = Float(probability)
            rescue StandardError
              fail ArgumentError, 'Costs and probabilities must be numbers'
            end
            fail ArgumentError, 'Costs must be non-negative integers' if integer_cost.negative?
            fail ArgumentError, 'Probabilities must be non-negative' if value.negative?

            probabilities[integer_cost] = probabilities.fetch(integer_cost, 0.0) + value
          end

          total = probabilities.values.reduce(0.0, :+)
          fail ArgumentError, 'Probabilities must sum to one' unless total.positive? && (total - 1.0).abs <= 1e-12

          probabilities.each_with_object({}) do |(cost, probability), normalized|
            normalized_probability = probability / total
            normalized[cost] = normalized_probability unless normalized_probability.zero?
          end
        end

        def validate_lambda(lambda_val)
          value = Float(lambda_val) rescue fail(ArgumentError, 'Lambda must be a finite non-negative number')
          fail ArgumentError, 'Lambda must be non-negative' unless value >= 0.0 && value.finite?
        end

        def validate_range(range_or_min, max_cost)
          if max_cost.nil?
            fail ArgumentError, 'Range must be given' if range_or_min.nil?

            min_cost = Integer(range_or_min.min) rescue fail(ArgumentError, 'Range must contain integer endpoints')
            max_cost = Integer(range_or_min.max) rescue fail(ArgumentError, 'Range must contain integer endpoints')
          else
            min_cost = Integer(range_or_min)
            max_cost = Integer(max_cost)
          end

          fail ArgumentError, 'Range minimum cannot exceed maximum' if min_cost > max_cost
          min_cost..max_cost
        end

        def probability_below(min_cost, cumulative)
          return 1.0 if min_cost > cumulative.length - 1

          cumulative_probability(min_cost - 1, cumulative)
        end

        def cumulative_probability(cost, cumulative)
          return 0.0 if cost.negative?
          return cumulative.last if cost >= cumulative.length

          cumulative[cost]
        end
      end
    end
  end
end
