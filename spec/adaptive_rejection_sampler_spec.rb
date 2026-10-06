require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::AdaptiveRejectionSampler do
  # Deterministic stand-in for Random: replays scripted values and then
  # keeps returning the last one, so tests can force specific proposals.
  class ScriptedRng
    def initialize(values, fallback = 0.5)
      @values = values.dup
      @last = fallback
    end

    def rand
      @last = @values.shift unless @values.empty?
      @last
    end
  end

  def normal_sampler(seed)
    described_class.new(
      log_pdf: ->(x) { -0.5 * x * x }, # unnormalized: no 1/sqrt(2*pi) needed
      dlog_pdf: ->(x) { -x },
      lower: -Float::INFINITY,
      upper: Float::INFINITY,
      points: [-1.5, 0.0, 1.5],
      rng: Random.new(seed)
    )
  end

  def sample_stats(samples)
    mean = samples.sum / samples.size
    var = samples.sum { |x| (x - mean)**2 } / samples.size
    [mean, var]
  end

  describe 'standard normal on the whole real line' do
    it 'matches mean, variance and central mass with a fixed seed' do
      samples = normal_sampler(42).sample(4000)
      mean, var = sample_stats(samples)
      central = samples.count { |x| x.abs < 1.0 }.to_f / samples.size

      expect(mean).to be_within(0.1).of(0.0)
      expect(var).to be_within(0.15).of(1.0)
      expect(central).to be_within(0.05).of(0.6827)
    end

    it 'is reproducible for a caller-supplied seeded rng' do
      expect(normal_sampler(7).sample(50)).to eq(normal_sampler(7).sample(50))
    end

    it 'returns a Float for single draws and an Array for counted draws' do
      sampler = normal_sampler(11)
      expect(sampler.sample).to be_a(Float)
      draws = sampler.sample(10)
      expect(draws.size).to eq(10)
      expect(draws).to all(be_a(Float))
    end

    it 'refines the hull as proposals are evaluated' do
      sampler = normal_sampler(5)
      initial_size = sampler.size
      sampler.sample(200)
      expect(sampler.size).to be > initial_size
    end
  end

  describe 'gamma(3, 1) on the one-sided unbounded support (0, inf)' do
    let(:sampler) do
      described_class.new(
        log_pdf: ->(x) { 2.0 * Math.log(x) - x },
        dlog_pdf: ->(x) { 2.0 / x - 1.0 },
        lower: 0.0,
        upper: Float::INFINITY,
        points: [1.0, 3.0, 6.0],
        rng: Random.new(123)
      )
    end

    it 'matches mean and variance with a fixed seed' do
      samples = sampler.sample(4000)
      mean, var = sample_stats(samples)

      expect(samples).to all(be > 0.0)
      expect(mean).to be_within(0.3).of(3.0)
      expect(var).to be_within(0.6).of(3.0)
    end
  end

  describe 'beta(2, 2) on the bounded support (0, 1)' do
    it 'stays inside the support and matches mean and variance' do
      sampler = described_class.new(
        log_pdf: ->(x) { Math.log(x) + Math.log(1.0 - x) },
        dlog_pdf: ->(x) { 1.0 / x - 1.0 / (1.0 - x) },
        lower: 0.0,
        upper: 1.0,
        points: [0.25, 0.5, 0.75],
        rng: Random.new(99)
      )
      samples = sampler.sample(4000)
      mean, var = sample_stats(samples)

      expect(samples).to all(be > 0.0)
      expect(samples).to all(be < 1.0)
      expect(mean).to be_within(0.05).of(0.5)
      expect(var).to be_within(0.02).of(0.05)
    end
  end

  describe 'log-concavity enforcement' do
    it 'refuses a heavy-tailed Student t density at initialization' do
      # t with 2 df is not log-concave: tangent intersections fall outside
      # their bracketing abscissae.
      expect do
        described_class.new(
          log_pdf: ->(x) { -1.5 * Math.log(1.0 + x * x) },
          dlog_pdf: ->(x) { -3.0 * x / (1.0 + x * x) },
          lower: -Float::INFINITY,
          upper: Float::INFINITY,
          points: [-2.0, 0.0, 2.0],
          rng: Random.new(1)
        )
      end.to raise_error(described_class::NotLogConcaveError)
    end

    it 'refuses a convex log-density' do
      expect do
        described_class.new(
          log_pdf: ->(x) { x * x },
          dlog_pdf: ->(x) { 2.0 * x },
          lower: -1.0,
          upper: 1.0,
          points: [-0.5, 0.5],
          rng: Random.new(1)
        )
      end.to raise_error(described_class::NotLogConcaveError)
    end

    it 'stops instead of sampling when the hull is violated mid-run' do
      # 1 / (1 + x^4) is log-concave around the origin but not in the
      # tails; the scripted rng forces a proposal deep into the right tail
      # where the tangent hull falls below the density.
      sampler = described_class.new(
        log_pdf: ->(x) { -Math.log(1.0 + x**4) },
        dlog_pdf: ->(x) { -4.0 * x**3 / (1.0 + x**4) },
        lower: -Float::INFINITY,
        upper: Float::INFINITY,
        points: [-1.0, 0.0, 1.0],
        rng: ScriptedRng.new([0.99, 0.9999, 0.5])
      )

      expect { sampler.sample }.to raise_error(described_class::NotLogConcaveError)
      # A failed sampler must keep failing, not resume with a broken hull.
      expect { sampler.sample }.to raise_error(described_class::Error)
    end
  end

  describe 'hull integrability' do
    it 'rejects a left tail whose tangent does not rise' do
      expect do
        described_class.new(
          log_pdf: ->(x) { -x },
          dlog_pdf: ->(_x) { -1.0 },
          lower: -Float::INFINITY,
          upper: Float::INFINITY,
          points: [0.0, 1.0],
          rng: Random.new(1)
        )
      end.to raise_error(described_class::NonIntegrableHullError)
    end

    it 'rejects a right tail whose tangent does not fall' do
      expect do
        described_class.new(
          log_pdf: ->(x) { x },
          dlog_pdf: ->(_x) { 1.0 },
          lower: -Float::INFINITY,
          upper: Float::INFINITY,
          points: [0.0, 1.0],
          rng: Random.new(1)
        )
      end.to raise_error(described_class::NonIntegrableHullError)
    end
  end

  describe 'attempt limiting' do
    it 'raises after the caller-supplied number of proposals' do
      sampler = described_class.new(
        log_pdf: ->(x) { -0.5 * x * x },
        dlog_pdf: ->(x) { -x },
        lower: -Float::INFINITY,
        upper: Float::INFINITY,
        points: [-1.0, 1.0],
        rng: ScriptedRng.new([0.999999]),
        max_attempts: 3
      )

      expect { sampler.sample }.to raise_error(described_class::MaxAttemptsExceededError)
    end
  end

  describe 'input validation' do
    it 'rejects starting points outside the support' do
      expect do
        described_class.new(
          log_pdf: ->(x) { -x * x },
          dlog_pdf: ->(x) { -2.0 * x },
          lower: 0.0,
          upper: 1.0,
          points: [0.5, 1.5],
          rng: Random.new(1)
        )
      end.to raise_error(ArgumentError)
    end

    it 'rejects non-positive max_attempts' do
      expect do
        described_class.new(
          log_pdf: ->(x) { -x * x },
          dlog_pdf: ->(x) { -2.0 * x },
          lower: -1.0,
          upper: 1.0,
          points: [-0.5, 0.5],
          rng: Random.new(1),
          max_attempts: 0
        )
      end.to raise_error(ArgumentError)
    end
  end
end
