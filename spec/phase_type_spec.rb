require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::PhaseType do
  let(:tolerance) { 1e-10 }

  describe Distribution::PhaseType::Ruby_ do
    it 'is the pure Ruby implementation used by Distribution::PhaseType.new' do
      t = [[-1.0]]
      distribution = Distribution::PhaseType.new([1.0], t)

      expect(distribution).to be_a(Distribution::PhaseType::Ruby_)
    end

    it 'reproduces an Erlang distribution for tandem exponential phases' do
      rate = 2.0
      shape = 3
      t = Array.new(shape) { Array.new(shape, 0.0) }
      (shape - 1).times do |i|
        t[i][i] = -rate
        t[i][i + 1] = rate
      end
      t[shape - 1][shape - 1] = -rate
      distribution = Distribution::PhaseType.new([1.0, 0.0, 0.0], t)

      expect(distribution.absorption_probability).to be_within(tolerance).of(1.0)
      expect(distribution.infinite_probability).to be_within(tolerance).of(0.0)
      expect(distribution.mean).to be_within(tolerance).of(shape / rate)
      expect(distribution.variance).to be_within(tolerance).of(shape / rate**2)

      [0.0, 0.25, 0.75, 1.5].each do |x|
        expected_pdf = rate * (rate * x)**(shape - 1) * Math.exp(-rate * x) / Math.gamma(shape)
        expected_cdf = Distribution::Gamma::Ruby_.cdf(rate * x, shape.to_f, 1.0)

        expect(distribution.pdf(x)).to be_within(1e-9).of(expected_pdf)
        expect(distribution.cdf(x)).to be_within(1e-9).of(expected_cdf)
        expect(distribution.sf(x)).to be_within(1e-9).of(1.0 - expected_cdf)
        expect(distribution.cdf(x) + distribution.sf(x)).to be_within(1e-9).of(1.0)
      end

      probability = distribution.cdf(0.9)
      expect(distribution.quantile(probability, tolerance: 1e-11)).to be_within(1e-8).of(0.9)
    end

    it 'identifies mass trapped forever by a closed transient component' do
      t = [
        [0.0, 0.0],
        [0.0, -1.0]
      ]
      distribution = Distribution::PhaseType.new([0.5, 0.5], t)

      expect(distribution.reachable_states).to eq([1])
      expect(distribution.closed_states).to eq([0])
      expect(distribution.absorption_probability).to be_within(tolerance).of(0.5)
      expect(distribution.infinite_probability).to be_within(tolerance).of(0.5)
      expect(distribution.cdf(0.0)).to be_within(tolerance).of(0.0)
      expect(distribution.cdf(Float::INFINITY)).to be_within(tolerance).of(0.5)
      expect(distribution.sf(Float::INFINITY)).to be_within(tolerance).of(0.5)

      [0.0, 0.5, 1.0, 2.0].each do |x|
        expect(distribution.pdf(x)).to be_within(1e-9).of(0.5 * Math.exp(-x))
        expect(distribution.cdf(x)).to be_within(1e-9).of(0.5 * (1.0 - Math.exp(-x)))
        expect(distribution.sf(x)).to be_within(1e-9).of(0.5 * (1.0 + Math.exp(-x)))
      end

      expect(distribution.finite_moment(1)).to be_within(tolerance).of(0.5)
      expect(distribution.moment(1)).to be_within(tolerance).of(1.0)
      expect(distribution.moment(1, conditional: false)).to eq(Float::INFINITY)
      expect(distribution.quantile(0.75)).to eq(Float::INFINITY)
      expect(distribution.quantile(0.75, details: true)[:status]).to eq(:beyond_absorption_probability)
    end

    it 'allows an initial atom already in the absorbing state' do
      t = [[-1.0]]
      distribution = Distribution::PhaseType.new([0.25], t, 0.75)

      expect(distribution.cdf(0.0)).to be_within(tolerance).of(0.75)
      expect(distribution.sf(0.0)).to be_within(tolerance).of(0.25)
      expect(distribution.pdf(0.0)).to be_within(tolerance).of(0.25)
      expect(distribution.quantile(0.8)).to be_within(1e-9).of(Math.log(0.25 / 0.2))
      expect(distribution.moment(1)).to be_within(tolerance).of(0.25)
    end

    it 'handles a near-defective repeated-rate generator without eigenvector inversion' do
      epsilon = 1e-4
      t = [
        [-1.0, epsilon],
        [0.0, -1.0]
      ]
      distribution = Distribution::PhaseType.new([1.0, 0.0], t)

      expect(distribution.mean).to be_within(1e-8).of(1.0 + epsilon)
      expect(distribution.variance).to be_within(1e-8).of(1.0 + 2 * epsilon - epsilon**2)

      x = 1.0
      transient_mass = Math.exp(-x) * (1.0 + epsilon * x)
      expected_pdf = Math.exp(-x) * (1.0 - epsilon + epsilon * x)

      expect(distribution.sf(x)).to be_within(1e-9).of(transient_mass)
      expect(distribution.cdf(x)).to be_within(1e-9).of(1.0 - transient_mass)
      expect(distribution.pdf(x)).to be_within(1e-9).of(expected_pdf)
    end

    it 'handles a completely closed generator' do
      distribution = Distribution::PhaseType.new([1.0], [[0.0]])

      expect(distribution.reachable_states).to eq([])
      expect(distribution.closed_states).to eq([0])
      expect(distribution.absorption_probability).to eq(0.0)
      expect(distribution.infinite_probability).to eq(1.0)
      expect(distribution.cdf(10.0)).to eq(0.0)
      expect(distribution.sf(10.0)).to eq(1.0)
      expect(distribution.finite_moment(2)).to eq(0.0)
      expect { distribution.mean }.to raise_error(ArgumentError, /conditional moments/)
    end

    it 'returns convergence and error information' do
      t = [[-2.0, 2.0], [0.0, -2.0]]
      distribution = Distribution::PhaseType.new([1.0, 0.0], t)
      result = distribution.cdf(1.0, details: true)

      expect(result[:status]).to eq(:converged)
      expect(result[:error_bound]).to be < 1e-9
      expect(result[:iterations]).to be_an(Array)
    end

    it 'handles negative time and probabilities beyond eventual absorption' do
      distribution = Distribution::PhaseType.new([0.5, 0.5], [[0.0, 0.0], [0.0, -1.0]])

      expect(distribution.pdf(-1.0)).to eq(0.0)
      expect(distribution.cdf(-1.0)).to eq(0.0)
      expect(distribution.sf(-1.0)).to eq(1.0)
      expect(distribution.quantile(0.5)).to eq(Float::INFINITY)
      expect(distribution.quantile(0.5, details: true)[:status])
        .to eq(:beyond_absorption_probability)
    end

    it 'rejects invalid explicit initial probability totals and quantile inputs' do
      expect { Distribution::PhaseType.new([0.5], [[-1.0]], 0.6) }
        .to raise_error(ArgumentError, /sum to 1/)
      expect { Distribution::PhaseType.new([1.0], [[-1.0]]).quantile(1.1) }
        .to raise_error(ArgumentError, /between 0 and 1/)
      expect { Distribution::PhaseType.new([1.0], [[-1.0]]).pdf(Float::NAN) }
        .to raise_error(ArgumentError, /NaN/)
    end

    it 'validates probabilities, off-diagonal rates, row sums, shape, and reachability' do
      expect { Distribution::PhaseType.new([1.0], [[-1.0, 0.0]]) }.to raise_error(ArgumentError, /square/)
      expect { Distribution::PhaseType.new([1.1], [[-1.0]]) }.to raise_error(ArgumentError, /exceed 1/)
      negative_off_diagonal = [
        [-1.0, -0.5],
        [0.0, -1.0]
      ]
      expect { Distribution::PhaseType.new([1.0, 0.0], negative_off_diagonal) }
        .to raise_error(ArgumentError, /off-diagonal/)

      positive_row_sum = [
        [-0.2, 0.6],
        [0.0, 0.0]
      ]
      expect { Distribution::PhaseType.new([1.0, 0.0], positive_row_sum) }
        .to raise_error(ArgumentError, /row 0 sum/)

      closed = Distribution::PhaseType.new([1.0], [[0.0]])
      expect(closed.closed_states).to eq([0])
      expect(closed.reachable_states).to eq([])
    end
  end
end
