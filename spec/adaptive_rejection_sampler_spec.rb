require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

class ScriptedRandom
  def initialize(values, fallback = Random.new(1))
    @values = values
    @fallback = fallback
  end

  def rand
    @values.empty? ? @fallback.rand : @values.shift
  end
end

class FixedRandom
  def initialize(value)
    @value = value
  end

  def rand
    @value
  end
end

describe Distribution::AdaptiveRejectionSampler do
  def normal_sampler(seed)
    described_class.new(
      log_density: ->(x) { -0.5 * x * x },
      derivative: ->(x) { -x },
      support: -Float::INFINITY..Float::INFINITY,
      starting_points: [-1.0, 0.5, 2.0],
      random: Random.new(seed),
      max_attempts: 100
    )
  end

  def sample_statistics(samples)
    mean = samples.sum / samples.size
    variance = samples.sum { |x| (x - mean)**2 } / (samples.size - 1)
    [mean, Math.sqrt(variance)]
  end

  it 'samples a standardized normal from an unnormalized log density' do
    sampler = normal_sampler(20_261_006)
    samples = Array.new(3000) { sampler.sample }
    mean, standard_deviation = sample_statistics(samples)

    expect(mean).to be_within(0.2).of(0.0)
    expect(standard_deviation).to be_within(0.18).of(1.0)
    expect(samples).to all(be_a(Float))
  end

  it 'samples Gamma(3, 1) on a half-unbounded support' do
    sampler = described_class.new(
      log_density: ->(x) { 2 * Math.log(x) - x },
      derivative: ->(x) { 2.0 / x - 1.0 },
      support: 0.0..Float::INFINITY,
      starting_points: [1.0, 4.0],
      random: Random.new(20_261_006),
      max_attempts: 100
    )
    samples = Array.new(3000) { sampler.sample }
    mean, standard_deviation = sample_statistics(samples)

    expect(mean).to be_within(0.3).of(3.0)
    expect(standard_deviation).to be_within(0.25).of(Math.sqrt(3))
    expect(samples).to all(be > 0.0)
  end

  it 'is reproducible with a fixed seed' do
    first = Array.new(20) { normal_sampler(42).sample }
    second = Array.new(20) { normal_sampler(42).sample }

    expect(first).to eq(second)
  end

  it 'updates knots and segment areas after a rejection' do
    random = ScriptedRandom.new([0.9, 0.99, 0.99, 0.9, 0.99, 0.1])
    sampler = normal_sampler(1)
    sampler.instance_variable_set(:@random, random)

    expect { sampler.sample }.to change(sampler, :rejections).by(1)
    expect(sampler.abscissae).to include(a_value_within(1e-12).of(3.5525850929940455))
    expect(sampler.segment_log_areas.size).to eq(sampler.abscissae.size)
  end

  it 'stops on a non-log-concave curvature violation' do
    log_density = lambda do |x|
      if x >= -1.0 && x <= 1.0
        -2.0 + x.abs
      else
        -x.abs
      end
    end
    derivative = lambda do |x|
      if x >= -1.0 && x <= 1.0
        x.negative? ? -1.0 : 1.0
      else
        x.negative? ? 1.0 : -1.0
      end
    end
    sampler = described_class.new(
      log_density: log_density,
      derivative: derivative,
      support: -Float::INFINITY..Float::INFINITY,
      starting_points: [-3.0, -2.0, 2.0, 3.0],
      random: ScriptedRandom.new([0.3, 0.315379628793682, 0.99]),
      max_attempts: 100
    )

    expect { sampler.sample }.to raise_error(described_class::NonLogConcaveError)
    expect(sampler.abscissae).to eq([-3.0, -2.0, 2.0, 3.0])
    expect(sampler.segment_log_areas.size).to eq(4)
    expect { sampler.sample }.to raise_error(described_class::NonLogConcaveError)
  end

  it 'stops when a squeeze chord lies above an evaluated proposal' do
    log_density = ->(x) { -(x * x) - 5.0 * Math.exp(-((x - 0.15) / 0.01)**2) }
    derivative = lambda do |x|
      distance = x - 0.15
      (-2.0 * x) + (100_000.0 * distance * Math.exp(-(distance / 0.01)**2))
    end
    sampler = described_class.new(
      log_density: log_density,
      derivative: derivative,
      support: -Float::INFINITY..Float::INFINITY,
      starting_points: [-2.0, 2.0],
      random: ScriptedRandom.new([0.75, 0.451_188_361_933_597_2, 0.5]),
      max_attempts: 100
    )

    expect { sampler.sample }.to raise_error(described_class::NonLogConcaveError, /squeezing chord/)
  end

  it 'stops when the initial tangents do not have integrable tails' do
    expect do
      described_class.new(
        log_density: ->(x) { x },
        derivative: ->(x) { 1.0 },
        support: 0.0..Float::INFINITY,
        starting_points: [1.0],
        random: Random.new(1),
        max_attempts: 100
      )
    end.to raise_error(described_class::NonIntegrableEnvelopeError)
  end

  it 'honors the caller supplied maximum attempt count' do
    sampler = described_class.new(
      log_density: ->(x) { -0.5 * x * x },
      derivative: ->(x) { -x },
      support: -Float::INFINITY..Float::INFINITY,
      starting_points: [-1.0, 1.0],
      random: FixedRandom.new(0.999_999),
      max_attempts: 2
    )

    expect { sampler.sample }.to raise_error(described_class::SamplingLimitError, 'exceeded 2 sampling attempts')
  end

  it 'requires positive slopes on both unbounded sides' do
    expect do
      described_class.new(
        log_density: ->(x) { -0.5 * x * x },
        derivative: ->(x) { -x },
        support: -Float::INFINITY..Float::INFINITY,
        starting_points: [1.0, 2.0],
        random: Random.new(1),
        max_attempts: 10
      )
    end.to raise_error(described_class::NonIntegrableEnvelopeError, /leftmost derivative/)
  end

  it 'requires two starting points on a doubly unbounded support' do
    expect do
      described_class.new(
        log_density: ->(x) { -0.5 * x * x },
        derivative: ->(x) { -x },
        support: -Float::INFINITY..Float::INFINITY,
        starting_points: [0.0],
        random: Random.new(1),
        max_attempts: 10
      )
    end.to raise_error(ArgumentError, /two starting points/)
  end
end
