module Distribution
  # Adaptive rejection sampling for differentiable, log-concave densities on
  # an interval. The supplied log density does not need to be normalized.
  class AdaptiveRejectionSampler
    class Error < StandardError; end
    class NonLogConcaveError < Error; end
    class NonIntegrableEnvelopeError < Error; end
    class SamplingLimitError < Error; end

    Node = Struct.new(:x, :h, :d)
    Segment = Struct.new(:node, :left, :right, :log_area)

    EPSILON = 1e-12

    attr_reader :abscissae, :attempts, :rejections, :segment_log_areas

    # @param log_density [#call] log of the unnormalized target density
    # @param derivative [#call] derivative of +log_density+
    # @param support [Range, Array] interval containing the target support
    # @param starting_points [Array<Numeric>] initial tangent abscissae
    # @param random [#rand] caller-controlled uniform random source
    # @param max_attempts [Integer] maximum proposals before one sample fails
    def initialize(log_density:, derivative:, support:, starting_points:, random:, max_attempts:)
      @log_density = log_density
      @derivative = derivative
      @lower, @upper = parse_support(support)
      @random = random
      @max_attempts = validate_positive_integer(max_attempts, :max_attempts)

      points = Array(starting_points).map(&:to_f).sort
      validate_starting_points(points)

      @nodes = points.map { |x| evaluate(x) }
      validate_log_concavity
      validate_outer_derivatives
      rebuild_envelope

      @rejections = 0
      @attempts = 0
      @stopped = false
    end

    # Draws one independent sample. Rejected proposals are added as knots and
    # the piecewise exponential envelope areas are rebuilt.
    def sample
      raise NonLogConcaveError, 'sampler has stopped because the target is not log-concave' if @stopped

      @attempts = 0
      loop do
        @attempts += 1
        if @attempts > @max_attempts
          raise SamplingLimitError, "exceeded #{@max_attempts} sampling attempts"
        end

        x = proposal
        h = @log_density.call(x)
        segment = segment_at(x)
        upper = tangent_value(segment.node, x)
        ensure_finite(h, "log density at proposal #{x}")

        unless h <= upper + value_tolerance(h, upper)
          stop!('proposal lies above its tangent upper envelope')
        end
        upper = [upper, h].max

        squeezed, squeeze = squeeze_value(x)
        uniform = @random.rand
        if squeezed
          if squeeze > h + value_tolerance(squeeze, h)
            stop!('squeezing chord lies above the log density')
          end
          squeeze = [squeeze, h].min
          return x if uniform <= Math.exp(squeeze - h)
        end

        return x if uniform <= Math.exp(h - upper)

        derivative = @derivative.call(x)
        ensure_finite(derivative, "log density derivative at proposal #{x}")
        add_rejected_point(x, h, derivative)
      end
    end

    private

    def parse_support(support)
      if support.is_a?(Range)
        lower = support.begin
        upper = support.end
      elsif support.is_a?(Array) && support.length == 2
        lower, upper = support
      else
        fail ArgumentError, 'support must be a Range or a two-element [lower, upper] array'
      end

      lower = -Float::INFINITY if lower.nil?
      upper = Float::INFINITY if upper.nil?
      lower = lower.to_f
      upper = upper.to_f

      unless lower < upper
        fail ArgumentError, 'support lower bound must be strictly less than upper bound'
      end

      [lower, upper]
    end

    def validate_positive_integer(value, name)
      unless value.is_a?(Integer) && value.positive?
        fail ArgumentError, "#{name} must be a positive integer"
      end
      value
    end

    def validate_starting_points(points)
      if @lower.infinite? && @upper.infinite? && points.length < 2
        fail ArgumentError, 'at least two starting points are required on an unbounded interval'
      elsif points.empty?
        fail ArgumentError, 'at least one starting point is required'
      end

      points.each_cons(2) do |left, right|
        if right - left <= EPSILON * [1.0, left.abs, right.abs].max
          fail ArgumentError, 'starting points must be distinct and ordered'
        end
      end

      points.each do |x|
        unless x > @lower && x < @upper
          fail ArgumentError, "starting point #{x} must be strictly inside the support"
        end
      end
    end

    def evaluate(x)
      h = @log_density.call(x)
      d = @derivative.call(x)
      ensure_finite(h, "log density at #{x}")
      ensure_finite(d, "log density derivative at #{x}")
      Node.new(x, h, d)
    end

    def ensure_finite(value, description)
      unless value.is_a?(Numeric) && value.respond_to?(:finite?) && value.finite?
        stop!("#{description} must be finite")
      end
      value.to_f
    end

    def validate_log_concavity
      @nodes.each_cons(2) do |left, right|
        if right.d > left.d + derivative_tolerance(left.d, right.d)
          stop!("derivative increased from #{left.d} to #{right.d}")
        end

        intersection = intersection(left, right)
        unless intersection >= left.x && intersection <= right.x
          stop!('tangent intersection lies outside the adjacent abscissae')
        end

        check_tangent_envelope(left, right)
        check_tangent_envelope(right, left)
        check_chord(left, right)
      end

      @nodes.each_cons(3) do |left, middle, right|
        left_secant = (middle.h - left.h) / (middle.x - left.x)
        right_secant = (right.h - middle.h) / (right.x - middle.x)
        if right_secant > left_secant + derivative_tolerance(left_secant, right_secant)
          stop!('log-density secant slopes increased')
        end
      end
    end

    def check_tangent_envelope(tangent, point)
      value = tangent_value(tangent, point.x)
      if point.h > value + value_tolerance(point.h, value)
        stop!("tangent at #{tangent.x} is below the log density at #{point.x}")
      end
    end

    def check_chord(left, right)
      midpoint = (left.x + right.x) / 2.0
      [left.x, midpoint, right.x].each do |x|
        h = x == midpoint ? @log_density.call(x) : (x == left.x ? left.h : right.h)
        chord = left.h + (right.h - left.h) * (x - left.x) / (right.x - left.x)
        if h.is_a?(Numeric) && h.finite? && h + value_tolerance(h, chord) < chord
          stop!("squeezing chord is above the log density at #{x}")
        end
      end
    end

    def validate_outer_derivatives
      if @lower.infinite? && @nodes.first.d <= 0
        @stopped = true
        raise NonIntegrableEnvelopeError, 'leftmost derivative must be positive for a finite left envelope tail'
      end
      if @upper.infinite? && @nodes.last.d >= 0
        @stopped = true
        raise NonIntegrableEnvelopeError, 'rightmost derivative must be negative for a finite right envelope tail'
      end
    end

    def rebuild_envelope
      validate_log_concavity
      validate_outer_derivatives

      boundaries = [@lower]
      @nodes.each_cons(2) { |left, right| boundaries << intersection(left, right) }
      boundaries << @upper

      @segments = @nodes.each_with_index.map do |node, index|
        left = boundaries[index]
        right = boundaries[index + 1]
        log_area = segment_log_area(node, left, right)
        Segment.new(node, left, right, log_area)
      end
      @abscissae = @nodes.map(&:x)
      @segment_log_areas = @segments.map(&:log_area)
    end

    def intersection(left, right)
      if (left.d - right.d).abs <= derivative_tolerance(left.d, right.d)
        first_intercept = left.h - left.d * left.x
        second_intercept = right.h - right.d * right.x
        unless (first_intercept - second_intercept).abs <= value_tolerance(first_intercept, second_intercept)
          stop!('parallel tangents do not form an upper envelope')
        end
        return (left.x + right.x) / 2.0
      end

      (left.h - right.h - left.d * left.x + right.d * right.x) / (right.d - left.d)
    end

    def segment_log_area(node, left, right)
      if left.infinite? && right.finite?
        unless node.d > 0
          raise NonIntegrableEnvelopeError, 'upper envelope is not integrable on the left tail'
        end
        node.h + node.d * (right - node.x) - Math.log(node.d)
      elsif left.finite? && right.infinite?
        unless node.d < 0
          raise NonIntegrableEnvelopeError, 'upper envelope is not integrable on the right tail'
        end
        node.h + node.d * (left - node.x) - Math.log(-node.d)
      else
        width = right - left
        h_at_left = node.h + node.d * (left - node.x)
        scaled_width = node.d * width
        if scaled_width.abs <= EPSILON
          h_at_left + Math.log(width)
        elsif node.d.positive?
          h_at_left + Math.log(Math.expm1(scaled_width)) - Math.log(node.d)
        else
          node.h + node.d * (right - node.x) +
            Math.log(Math.expm1(-scaled_width)) - Math.log(-node.d)
        end
      end
    rescue NonIntegrableEnvelopeError
      @stopped = true
      raise
    end

    def proposal
      segment = select_segment
      inverse_segment_cdf(segment, @random.rand)
    end

    def select_segment
      maximum_log_area = @segments.map(&:log_area).max
      weights = @segments.map { |segment| Math.exp(segment.log_area - maximum_log_area) }
      total_area = weights.inject(0.0, :+)
      target = @random.rand * total_area
      cumulative = 0.0
      weights.each_with_index do |weight, index|
        cumulative += weight
        return @segments[index] if target <= cumulative
      end
      @segments.last
    end

    def inverse_segment_cdf(segment, uniform)
      v = open_unit_interval(uniform)
      node = segment.node
      left = segment.left
      right = segment.right
      d = node.d

      if left.infinite?
        right + Math.log(v) / d
      elsif right.infinite?
        left + Math.log(1.0 - v) / d
      elsif (d * (right - left)).abs <= EPSILON
        left + (right - left) * v
      else
        a = d * (left - node.x)
        b = d * (right - node.x)
        if d.positive?
          log_term = b + Math.log(v + (1.0 - v) * Math.exp(a - b))
        else
          log_term = a + Math.log(1.0 - v + v * Math.exp(b - a))
        end
        node.x + log_term / d
      end
    end

    def open_unit_interval(value)
      value = value.to_f
      return Float::EPSILON unless value.positive?
      return 1.0 - Float::EPSILON unless value < 1.0
      value
    end

    def segment_at(x)
      segment = @segments.find { |candidate| x >= candidate.left && x <= candidate.right }
      return segment unless segment.nil?
      stop!("proposal #{x} is outside the tangent upper envelope")
    end

    def tangent_value(node, x)
      node.h + node.d * (x - node.x)
    end

    def squeeze_value(x)
      right_index = @abscissae.index { |knot| knot >= x }
      return [false, nil] if right_index.nil? || right_index.zero?

      left = @nodes[right_index - 1]
      right = @nodes[right_index]
      if x == right.x
        [true, right.h]
      else
        value = left.h + (right.h - left.h) * (x - left.x) / (right.x - left.x)
        [true, value]
      end
    end

    def add_rejected_point(x, h, derivative)
      @rejections += 1
      index = @abscissae.bsearch_index { |knot| knot >= x } || @abscissae.length
      if index < @abscissae.length &&
         (x - @abscissae[index]).abs <= x_tolerance(x, @abscissae[index])
        return
      end

      @nodes.insert(index, Node.new(x, h, derivative))
      rebuild_envelope
    end

    def value_tolerance(*values)
      EPSILON * ([1.0] + values.map(&:abs)).max
    end

    def derivative_tolerance(*values)
      EPSILON * ([1.0] + values.map(&:abs)).max
    end

    def x_tolerance(*values)
      EPSILON * ([1.0] + values.map(&:abs)).max
    end

    def stop!(message)
      @stopped = true
      raise NonLogConcaveError, message
    end
  end
end
