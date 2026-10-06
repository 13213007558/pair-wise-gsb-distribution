module Distribution
  class Mixture
    Component = Struct.new(:source, :weight, :arguments, :options) do
      def kind
        return :point if source.is_a?(Distribution::PointMass)
        return :mixed if source.is_a?(Distribution::Mixture)
        return options[:type].to_sym if options[:type]
        return :discrete if known_discrete_source?

        if source.respond_to?(:mass_at)
          source.respond_to?(:density) ? :mixed : :discrete
        else
          :continuous
        end
      end

      def known_discrete_source?
        [Distribution::Binomial, Distribution::Hypergeometric, Distribution::Poisson].include?(source)
      end

      def discrete?
        kind == :discrete || kind == :mixed || kind == :point
      end

      def density(x)
        return 0.0 unless kind == :continuous || kind == :mixed

        invoke(source.respond_to?(:density) ? :density : :pdf, x).to_f
      end

      def cdf(x)
        clamp_unit(invoke(:cdf, x).to_f)
      end

      def cdf_before(x)
        cdf(x) - mass_at(x)
      end

      def mass_at(x)
        case kind
        when :point
          source.location == x ? 1.0 : 0.0
        when :discrete, :mixed
          if source.respond_to?(:mass_at)
            invoke(:mass_at, x).to_f
          else
            invoke(:pdf, x).to_f
          end
        else
          0.0
        end
      end

      def quantile(p)
        return source.quantile(p) if source.is_a?(Distribution::Mixture)

        invoke(source.respond_to?(:quantile) ? :quantile : :p_value, p)
      end

      def support
        return options[:support] if options.key?(:support)
        return source.support if source.respond_to?(:support)

        inferred_support
      end

      def finite_events
        values =
          if options[:events]
            Array(options[:events])
          elsif source.respond_to?(:events)
            Array(source.events)
          elsif kind == :point
            [source.location]
          elsif kind == :discrete
            discrete_support_events
          else
            continuous_support_events
          end

        values.select { |value| value.is_a?(Numeric) && value.finite? }.uniq.sort
      end

      private

      def invoke(method, value)
        args = arguments_for(method)
        if source.respond_to?(method)
          source.public_send(method, value, *args)
        elsif source.respond_to?(:call)
          source.call(value)
        else
          fail NoMethodError, "#{source.inspect} does not implement #{method}"
        end
      end

      def arguments_for(method)
        args = arguments.is_a?(Hash) ? arguments.fetch(method, []) : arguments
        args.nil? ? [] : Array(args)
      end

      def inferred_support
        return nil unless source.respond_to?(:quantile) || source.respond_to?(:p_value)

        lower = quantile(0.0)
        upper = quantile(1.0)
        return nil unless valid_bound?(lower) && valid_bound?(upper)

        lower..upper
      rescue StandardError
        nil
      end

      def valid_bound?(value)
        value.is_a?(Numeric) && !(value.is_a?(Float) && value.nan?)
      end

      def discrete_support_events
        case support
        when Array then support
        when Range
          support.to_a
        else
          []
        end
      rescue RangeError, TypeError
        []
      end

      def continuous_support_events
        values =
          case support
          when Array then support.minmax
          when Range then [support.min, support.max]
          else
            []
          end

        values.select { |value| value.is_a?(Numeric) && value.finite? }
      end

      def clamp_unit(value)
        return 0.0 if value < 0 && value > -1e-12
        return 1.0 if value > 1 && value < 1 + 1e-12

        value
      end
    end

    attr_reader :components, :point_masses, :total_weight

    def initialize(components = [], point_masses: {}, weights: nil, normalize: true, normalize_weights: normalize)
      entries = Array(components)
      entries = entries.each_with_index.map { |component, index| [component, weights[index]] } if weights

      @raw_components = []
      merged_points = Hash.new(0)

      entries.each { |entry| merge_entry(entry, merged_points) }
      coerce_point_masses(point_masses).each { |location, weight| merged_points[location] += weight }

      reject_invalid_weights(merged_points.values, @raw_components.map { |component| component[:weight] })
      @total_weight = merged_points.values.sum + @raw_components.sum { |component| component[:weight] }
      fail ArgumentError, 'mixture weights must sum to a positive value' unless @total_weight.positive?

      unless normalize_weights || (@total_weight.to_f - 1.0).abs <= 1e-12
        fail ArgumentError, 'unnormalized mixture weights must sum to 1'
      end

      scale = normalize_weights ? 1.0 / @total_weight : 1.0
      @point_masses = merged_points.transform_values { |weight| weight * scale }
      @components = @raw_components.map do |component|
        Component.new(
          component[:source],
          component[:weight] * scale,
          component[:arguments],
          component[:options]
        )
      end
    end

    def density(x)
      components.sum(0.0) { |component| component.weight * component.density(x) }
    end

    def mass_at(x)
      explicit = point_masses.sum(0.0) { |location, weight| location == x ? weight : 0.0 }
      component_mass = components.sum(0.0) { |component| component.weight * component.mass_at(x) }
      explicit + component_mass
    end

    def cdf(x)
      value = point_masses.sum(0.0) { |location, weight| location <= x ? weight : 0.0 }
      value += components.sum(0.0) { |component| component.weight * component.cdf(x) }
      clamp_probability(value)
    end

    def cdf_before(x)
      cdf(x) - mass_at(x)
    end

    def quantile(p)
      probability = validate_probability(p)
      return lower_support if probability.zero?
      return upper_support if probability == 1.0

      events.each do |event|
        before = cdf_before(event)
        return solve_interval(previous_event(event), event, probability) if probability <= before
        return event if probability <= cdf(event)
      end

      solve_interval(events.last || -Float::INFINITY, Float::INFINITY, probability)
    end

    alias_method :p_value, :quantile
    alias_method :inverse_cdf, :quantile
    alias_method :pdf, :density

    def events
      merge_equal_locations(components.flat_map(&:finite_events).concat(point_masses.keys))
    end

    def support
      lower = lower_support
      upper = upper_support
      lower..upper if lower && upper
    end

    private

    def merge_entry(entry, merged_points)
      case entry
      when Hash
        merge_hash(entry, merged_points)
      when Numeric
        add_point(entry, 1, merged_points)
      when Distribution::PointMass
        add_point(entry.location, 1, merged_points)
      when Array
        merge_tuple(entry, merged_points)
      else
        add_component(entry, 1, [])
      end
    end

    def merge_hash(entry, merged_points)
      if entry.key?(:location) || entry.key?(:point) || entry.key?(:at) || entry.key?(:value)
        location = entry.fetch(:location) { entry.fetch(:point) { entry.fetch(:at) { entry[:value] } } }
        add_point(location, entry.fetch(:weight, 1), merged_points)
      else
        source = entry.fetch(:distribution) { entry.fetch(:component) { entry[:source] } }
        arguments = entry[:args] || entry[:params] || entry[:arguments] || []
        options = entry.slice(:type, :support, :events)
        weight = entry.fetch(:weight, 1)
        if source.is_a?(Distribution::PointMass)
          add_point(source.location, weight, merged_points)
        else
          add_component(source, weight, arguments, options)
        end
      end
    end

    def merge_tuple(entry, merged_points)
      source, weight = entry
      arguments = entry[2] || []
      options = entry[3].is_a?(Hash) ? entry[3].slice(:type, :support, :events) : {}

      if source.is_a?(Numeric)
        add_point(source, weight || 1, merged_points)
      elsif source.is_a?(Distribution::PointMass)
        add_point(source.location, weight || 1, merged_points)
      else
        add_component(source, weight || 1, arguments, options)
      end
    end

    def add_point(location, weight, merged_points)
      validate_weight(weight)
      unless location.is_a?(Numeric) && finite_numeric?(location)
        fail ArgumentError, 'point-mass locations must be finite numeric values'
      end

      existing = merged_points.keys.find { |candidate| candidate == location }
      key =
        if existing && integer_location?(location) && !integer_location?(existing)
          location
        else
          existing || location
        end
      merged_points.delete(existing) if existing && key != existing
      merged_points[key] += weight
    end

    def merge_equal_locations(locations)
      locations.select { |location| location.is_a?(Numeric) }.sort_by(&:to_f)
               .chunk_while { |left, right| left == right }.map(&:first)
    end

    def integer_location?(location)
      location.is_a?(Integer)
    end

    def add_component(source, weight, arguments, options = {})
      validate_weight(weight)
      @raw_components << {
        source: source,
        weight: weight,
        arguments: arguments,
        options: options
      }
    end

    def coerce_point_masses(point_masses)
      case point_masses
      when Hash
        point_masses.to_a
      when Distribution::PointMass
        [[point_masses.location, 1]]
      else
        Array(point_masses).map do |point|
          point.is_a?(Distribution::PointMass) ? [point.location, 1] : Array(point)
        end
      end
    end

    def reject_invalid_weights(*weights)
      weights.flatten.each { |weight| validate_weight(weight) }
    end

    def validate_weight(weight)
      unless weight.is_a?(Numeric) && finite_numeric?(weight) && weight >= 0
        fail ArgumentError, 'mixture weights must be non-negative numeric values'
      end
    end

    def clamp_probability(value)
      return 0.0 if value < 0 && value > -1e-12
      return 1.0 if value > 1 && value < 1 + 1e-12

      value
    end

    def validate_probability(p)
      unless p.is_a?(Numeric) && finite_numeric?(p) && p >= 0 && p <= 1
        fail RangeError, "quantile probability must be in [0, 1], got #{p.inspect}"
      end

      p.to_f
    end

    def finite_numeric?(value)
      !value.respond_to?(:nan?) || !value.nan?
    end

    def previous_event(event)
      index = events.index(event)
      index.zero? ? -Float::INFINITY : events[index - 1]
    end

    def solve_interval(left, right, p)
      active = active_components(left, right)
      total_probability = active.sum(0.0) { |component| interval_probability(component, left, right) }
      return left if total_probability.zero?

      if active.one?
        component = active.first
        return left if !left.infinite? && p <= cdf(left)
        start = left.infinite? ? 0.0 : component.cdf(left)
        target = start + (p - (left.infinite? ? 0.0 : cdf(left))) / component.weight
        value = component.quantile(clamp_component_probability(target))
        return clamp_interval_value(value, left, right)
      end

      bisect_interval(left, right, p)
    end

    def active_components(left, right)
      components.select do |component|
        next false if component.kind == :point

        interval_probability(component, left, right) > 1e-14
      end
    end

    def interval_probability(component, left, right)
      finish = right.infinite? ? 1.0 : component.cdf_before(right)
      start = left.infinite? ? 0.0 : component.cdf(left)
      [finish - start, 0.0].max
    end

    def clamp_component_probability(value)
      return 0.0 if value < 0 && value > -1e-12
      return 1.0 if value > 1 && value < 1 + 1e-12

      [[value, 0.0].max, 1.0].min
    end

    def clamp_interval_value(value, left, right)
      value = [value, left].max unless left.infinite?
      value = [value, right].min unless right.infinite?
      value
    end

    def bisect_interval(left, right, p)
      low, high = bisection_bounds(left, right)
      80.times do
        middle = (low + high) / 2.0
        if cdf(middle) >= p
          high = middle
        else
          low = middle
        end
      end
      high
    end

    def bisection_bounds(left, right)
      low =
        if left.infinite?
          candidate = -1.0
          candidate *= 2 while cdf(candidate).positive?
          candidate
        else
          left
        end

      high =
        if right.infinite?
          candidate = [1.0, events.last.to_f + 1.0].max
          candidate *= 2 while cdf(candidate) < 1.0
          candidate
        else
          right
        end

      [low, high]
    end

    def lower_support
      return -Float::INFINITY if components.any? { |component| unbounded_bound?(component, :lower) }

      bounds = components.filter_map { |component| finite_support_bound(component, :lower) }
      locations = events
      candidates = bounds + locations
      finite = candidates.select(&:finite?)
      finite.min || -Float::INFINITY
    end

    def upper_support
      return Float::INFINITY if components.any? { |component| unbounded_bound?(component, :upper) }

      bounds = components.filter_map { |component| finite_support_bound(component, :upper) }
      locations = events
      candidates = bounds + locations
      finite = candidates.select(&:finite?)
      finite.max || Float::INFINITY
    end

    def unbounded_bound?(component, side)
      support = component.support
      bound = side == :lower ? support&.min : support&.max
      bound.nil? || bound.infinite?
    end

    def finite_support_bound(component, side)
      support = component.support
      bound = side == :lower ? support&.min : support&.max
      bound if bound.is_a?(Numeric) && bound.finite?
    end

    class PointMass
      attr_reader :location

      def initialize(location)
        @location = location
      end

      def density(_x)
        0.0
      end

      def mass_at(x)
        x == location ? 1.0 : 0.0
      end

      def cdf(x)
        location <= x ? 1.0 : 0.0
      end

      def quantile(p)
        fail RangeError, "quantile probability must be in [0, 1], got #{p.inspect}" unless p.is_a?(Numeric) && p >= 0 && p <= 1

        location
      end

      alias_method :p_value, :quantile
      alias_method :pdf, :mass_at

      def support
        location..location
      end

      def events
        [location]
      end
    end
  end

  PointMass = Mixture::PointMass
end
