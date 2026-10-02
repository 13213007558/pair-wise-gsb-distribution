require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::NormalMultivariate do
  let(:rng) { Random.new(12_345) }

  def cdf(lower, upper, mean, sigma, **opts)
    Distribution::NormalMultivariate.cdf(
      lower, upper, mean, sigma,
      **{ epsilon: 0.005, max_iterations: 100_000, random: rng }.merge(opts)
    )
  end

  describe '.cdf' do
    it 'returns a Result with estimate, error, achieved flag and iterations' do
      result = cdf([0, 0], [1, 1], [0, 0], [[1, 0], [0, 1]])
      expect(result).to be_a(Distribution::NormalMultivariate::Result)
      expect(result.estimate).to be_between(0.0, 1.0)
      expect(result.error).to be >= 0.0
      expect([true, false]).to include(result.achieved)
      expect(result.iterations).to be >= 1
    end

    it 'approximates the product of univariate probabilities for independent components' do
      lower = [-1.0, 0.5, -Float::INFINITY]
      upper = [1.5, 2.0, 0.0]
      sigma = [[4.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 9.0]]
      mean = [0.5, 0.0, -1.0]
      expected = (0...3).map do |i|
        sd = Math.sqrt(sigma[i][i])
        lo = (lower[i] - mean[i]) / sd
        hi = (upper[i] - mean[i]) / sd
        Distribution::Normal.cdf(hi) - Distribution::Normal.cdf(lo)
      end.inject(:*)
      result = cdf(lower, upper, mean, sigma)
      expect(result.estimate).to be_within(0.02).of(expected)
    end

    it 'matches the analytic quadrant probability for a 2D zero-mean orthogonal normal' do
      result = cdf([0, 0], [Float::INFINITY, Float::INFINITY],
                   [0, 0], [[1, 0], [0, 1]])
      expect(result.estimate).to be_within(0.02).of(0.25)
    end

    it 'matches the analytic quadrant probability for a correlated 2D zero-mean normal' do
      rho = 0.6
      expected = 0.25 + Math.asin(rho) / (2 * Math::PI)
      result = cdf([0, 0], [Float::INFINITY, Float::INFINITY],
                   [0, 0], [[1, rho], [rho, 1]])
      expect(result.estimate).to be_within(0.02).of(expected)
    end

    it 'is exact for a one-dimensional problem' do
      result = cdf([-1.0], [2.0], [0.5], [[4.0]])
      expected = Distribution::Normal.cdf((2.0 - 0.5) / 2.0) -
                 Distribution::Normal.cdf((-1.0 - 0.5) / 2.0)
      expect(result.estimate).to be_within(1e-12).of(expected)
      expect(result.error).to eq(0.0)
      expect(result.achieved).to be true
    end

    it 'accepts infinite bounds on every dimension and returns 1 for the whole space' do
      result = cdf([-Float::INFINITY, -Float::INFINITY],
                   [Float::INFINITY, Float::INFINITY],
                   [1, -2], [[1, 0.3], [0.3, 2]])
      expect(result.estimate).to eq(1.0)
      expect(result.error).to eq(0.0)
    end

    it 'accepts nil bounds as unbounded' do
      result = cdf([nil, nil], [nil, nil], [0, 0], [[1, 0], [0, 1]])
      expect(result.estimate).to eq(1.0)
    end

    it 'rejects a non positive definite covariance matrix' do
      expect do
        cdf([0, 0], [1, 1], [0, 0], [[1, 2], [2, 1]])
      end.to raise_error(ArgumentError, /positive definite/)
    end

    it 'rejects mismatched dimensions' do
      expect do
        cdf([0, 0], [1, 1], [0], [[1, 0], [0, 1]])
      end.to raise_error(ArgumentError)
    end

    it 'rejects lower bounds above upper bounds' do
      expect do
        cdf([2], [1], [0], [[1]])
      end.to raise_error(ArgumentError)
    end

    it 'is reproducible for a fixed seed' do
      args = [[0, 0], [1, 1], [0, 0], [[1, 0.4], [0.4, 1]]]
      first = Distribution::NormalMultivariate.cdf(
        *args, epsilon: 0.005, max_iterations: 10_000, random: Random.new(99)
      )
      second = Distribution::NormalMultivariate.cdf(
        *args, epsilon: 0.005, max_iterations: 10_000, random: Random.new(99)
      )
      expect(first.estimate).to eq(second.estimate)
      expect(first.error).to eq(second.error)
      expect(first.iterations).to eq(second.iterations)
    end

    it 'reports achieved=false when the iteration budget is exhausted' do
      result = cdf([0, 0], [1, 1], [0, 0], [[1, 0.5], [0.5, 1]],
                   epsilon: 1e-12, max_iterations: 50)
      expect(result.iterations).to eq(50)
      expect(result.achieved).to be false
      expect(result.error).to be > 0.0
    end

    it 'reports achieved=true when the error target is met within budget' do
      result = cdf([0, 0], [1, 1], [0, 0], [[1, 0.5], [0.5, 1]],
                   epsilon: 0.05, max_iterations: 100_000)
      expect(result.achieved).to be true
      expect(result.error).to be <= 0.05
    end
  end
end
