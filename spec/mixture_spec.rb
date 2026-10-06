require File.expand_path(File.dirname(__FILE__) + "/spec_helper.rb")
require 'distribution/mixture'
require 'distribution/uniform'

# Adapter turning an existing gem distribution into a mixture component
# (an object responding to #pdf, #cdf and #quantile).
class UniformComponent
  def initialize(lower, upper)
    @lower = lower
    @upper = upper
  end

  def pdf(x)
    Distribution::Uniform::Ruby_.pdf(x, @lower, @upper)
  end

  def cdf(x)
    Distribution::Uniform::Ruby_.cdf(x, @lower, @upper)
  end

  def quantile(q)
    Distribution::Uniform::Ruby_.quantile(q, @lower, @upper)
  end
end

describe Distribution::Mixture do
  let(:uniform_a) { UniformComponent.new(2.0, 3.0) }
  let(:uniform_b) { UniformComponent.new(5.0, 6.0) }

  # Point mass at zero plus two disjoint continuous distributions.
  # F(x) = 0          for x < 0
  # F(0) = 0.5        (atom jump)
  # F rises 0.5->0.75 on [2, 3]
  # F(x) = 0.75       on [3, 5] (plateau over the gap)
  # F rises 0.75->1.0 on [5, 6]
  let(:mixture) do
    Distribution::Mixture.new(
      components: [[0.25, uniform_a], [0.25, uniform_b]],
      atoms: [[0.5, 0.0]]
    )
  end

  describe '#initialize' do
    it 'merges atoms at the same location before anything else' do
      mix = Distribution::Mixture.new(
        components: [[0.25, uniform_a], [0.25, uniform_b]],
        atoms: [[0.2, 0.0], [0.3, 0.0], [0.1, 4.0]],
        normalize: false
      )
      expect(mix.atoms).to eq([[0.0, 0.5], [4.0, 0.1]])
      expect(mix.mass_at(0.0)).to eq(0.5)
    end

    it 'accepts atoms given as a location => weight hash' do
      mix = Distribution::Mixture.new(atoms: { 1.0 => 0.5, 3.0 => 0.5 })
      expect(mix.mass_at(1.0)).to eq(0.5)
      expect(mix.mass_at(3.0)).to eq(0.5)
    end

    it 'normalizes weights to sum to 1 by default' do
      mix = Distribution::Mixture.new(
        components: [[1.0, uniform_a], [1.0, uniform_b]],
        atoms: [[2.0, 0.0]]
      )
      expect(mix.components.map(&:first)).to eq([0.25, 0.25])
      expect(mix.mass_at(0.0)).to eq(0.5)
      expect(mix.total_weight).to eq(1.0)
      expect(mix.cdf(100.0)).to eq(1.0)
    end

    it 'keeps the given weights with normalize: false' do
      mix = Distribution::Mixture.new(
        components: [[0.5, uniform_a]],
        atoms: [[1.0, 0.0]],
        normalize: false
      )
      expect(mix.total_weight).to eq(1.5)
      expect(mix.cdf(100.0)).to eq(1.5)
    end

    it 'rejects negative weights outright' do
      expect do
        Distribution::Mixture.new(components: [[-0.1, uniform_a]], atoms: [[1.1, 0.0]])
      end.to raise_error(ArgumentError)
      expect do
        Distribution::Mixture.new(components: [[1.0, uniform_a]], atoms: [[-1e-9, 0.0]])
      end.to raise_error(ArgumentError)
    end

    it 'rejects weights that cannot be normalized' do
      expect do
        Distribution::Mixture.new(components: [[0.0, uniform_a]])
      end.to raise_error(ArgumentError)
    end

    it 'rejects components without cdf/quantile' do
      expect do
        Distribution::Mixture.new(components: [[1.0, Object.new]])
      end.to raise_error(ArgumentError)
    end
  end

  describe '#cdf' do
    it 'is 0 below the support and 1 above it' do
      expect(mixture.cdf(-1.0)).to eq(0.0)
      expect(mixture.cdf(6.0)).to eq(1.0)
      expect(mixture.cdf(100.0)).to eq(1.0)
    end

    it 'jumps at the atom and is right-continuous there' do
      expect(mixture.cdf(0.0)).to eq(0.5)
      expect(mixture.cdf(-1e-12)).to eq(0.0)
    end

    it 'adds continuous mass between jumps' do
      expect(mixture.cdf(2.0)).to eq(0.5)
      expect(mixture.cdf(2.5)).to eq(0.625)
      expect(mixture.cdf(3.0)).to eq(0.75)
      expect(mixture.cdf(5.5)).to eq(0.875)
    end

    it 'stays flat on the gap between disjoint supports' do
      expect(mixture.cdf(3.5)).to eq(0.75)
      expect(mixture.cdf(4.9)).to eq(0.75)
      expect(mixture.cdf(5.0)).to eq(0.75)
    end
  end

  describe '#mass_at' do
    it 'returns the atom mass at atom locations and 0 elsewhere' do
      expect(mixture.mass_at(0.0)).to eq(0.5)
      expect(mixture.mass_at(2.5)).to eq(0.0)
      expect(mixture.mass_at(5.0)).to eq(0.0)
    end
  end

  describe '#density' do
    it 'is the weighted sum of component densities' do
      expect(mixture.density(2.5)).to eq(0.25)
      expect(mixture.density(5.5)).to eq(0.25)
    end

    it 'is 0 at atoms and inside gaps' do
      expect(mixture.density(0.0)).to eq(0.0)
      expect(mixture.density(4.0)).to eq(0.0)
    end

    it 'is aliased as #pdf' do
      expect(mixture.pdf(2.5)).to eq(mixture.density(2.5))
    end
  end

  describe '#quantile' do
    it 'returns the atom location for probabilities inside the jump' do
      expect(mixture.quantile(0.25)).to eq(0.0)
      expect(mixture.quantile(0.5)).to eq(0.0)
    end

    it 'jumps over the massless gap right after the atom' do
      expect(mixture.quantile(0.5 + 1e-9)).to be_within(1e-6).of(2.0)
      expect(mixture.quantile(0.5 + 1e-9)).to be > 2.0
    end

    it 'inverts the continuous parts' do
      expect(mixture.quantile(0.625)).to be_within(1e-9).of(2.5)
      expect(mixture.quantile(0.875)).to be_within(1e-9).of(5.5)
    end

    it 'returns the left edge of the plateau over the gap' do
      expect(mixture.quantile(0.75)).to be_within(1e-9).of(3.0)
    end

    it 'jumps across the gap when p crosses the plateau level' do
      below = mixture.quantile(0.75 - 1e-6)
      above = mixture.quantile(0.75 + 1e-6)
      expect(below).to be < 3.0
      expect(below).to be > 2.9
      expect(above).to be > 5.0
      expect(above).to be < 5.1
    end

    it 'follows the generalized inverse at the boundaries' do
      expect(mixture.quantile(0.0)).to eq(-Float::INFINITY)
      expect(mixture.quantile(1.0)).to eq(6.0)
    end

    it 'satisfies the generalized inverse property on a grid' do
      0.01.step(0.99, 0.01) do |p|
        x = mixture.quantile(p)
        expect(mixture.cdf(x)).to be >= p - 1e-12
        expect(mixture.cdf(x - 1e-8)).to be < p
      end
    end

    it 'is aliased as #p_value' do
      expect(mixture.p_value(0.625)).to eq(mixture.quantile(0.625))
    end

    it 'rejects probabilities outside [0, 1]' do
      expect { mixture.quantile(-0.1) }.to raise_error(RangeError)
      expect { mixture.quantile(1.1) }.to raise_error(RangeError)
    end
  end

  context 'with a single continuous component' do
    let(:mixture) do
      Distribution::Mixture.new(components: [[1.0, UniformComponent.new(0.0, 1.0)]])
    end

    it 'uses the component quantile exactly' do
      expect(mixture.quantile(0.25)).to eq(0.25)
      expect(mixture.cdf(0.5)).to eq(0.5)
      expect(mixture.density(0.5)).to eq(1.0)
    end
  end

  context 'with atoms only' do
    let(:mixture) do
      Distribution::Mixture.new(atoms: [[0.5, 1.0], [0.5, 3.0]])
    end

    it 'behaves like a discrete distribution' do
      expect(mixture.cdf(0.9)).to eq(0.0)
      expect(mixture.cdf(1.0)).to eq(0.5)
      expect(mixture.cdf(2.9)).to eq(0.5)
      expect(mixture.cdf(3.0)).to eq(1.0)
      expect(mixture.quantile(0.25)).to eq(1.0)
      expect(mixture.quantile(0.5)).to eq(1.0)
      expect(mixture.quantile(0.75)).to eq(3.0)
      expect(mixture.quantile(1.0)).to eq(3.0)
      expect(mixture.mass_at(3.0)).to eq(0.5)
      expect(mixture.density(1.0)).to eq(0.0)
    end
  end

  context 'with atoms splitting the continuous support' do
    # Atom between the two components: each inter-atom interval holds a
    # single varying component, so quantiles are exact component inverses.
    let(:mixture) do
      Distribution::Mixture.new(
        components: [[0.25, uniform_a], [0.25, uniform_b]],
        atoms: [[0.5, 4.0]]
      )
    end

    it 'inverts each component exactly' do
      expect(mixture.quantile(0.125)).to eq(2.5)
      expect(mixture.quantile(0.875)).to eq(5.5)
      expect(mixture.quantile(0.3)).to eq(4.0)
      expect(mixture.quantile(0.75)).to eq(4.0)
    end
  end
end
