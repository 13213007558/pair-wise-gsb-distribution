require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::NormalMultivariate do
  def cdf(lower, upper, mean:, covariance:, **options)
    described_class.cdf(
      lower,
      upper,
      mean: mean,
      covariance: covariance,
      random: options.delete(:random) || Random.new(42),
      **options
    )
  end

  def one_dimensional_probability(lower, upper, mean, variance)
    Distribution::Normal.cdf((upper - mean) / Math.sqrt(variance)) -
      Distribution::Normal.cdf((lower - mean) / Math.sqrt(variance))
  end

  it 'estimates independent dimensions as the one-dimensional product' do
    mean = [0.5, -0.5, 1.0]
    covariance = [[4.0, 0.0, 0.0], [0.0, 9.0, 0.0], [0.0, 0.0, 16.0]]
    lower = [-1.0, -2.0, -3.0]
    upper = [2.0, 1.0, 4.0]
    expected = mean.each_with_index.map do |dimension_mean, i|
      one_dimensional_probability(lower[i], upper[i], dimension_mean, covariance[i][i])
    end.reduce(:*)

    result = cdf(lower, upper, mean: mean, covariance: covariance)

    expect(result[:estimate]).to be_within(1e-12).of(expected)
    expect(result[:error_estimate]).to eq(0.0)
    expect(result[:target_reached]).to eq(true)
  end

  it 'matches the analytic independent zero-mean orthogonal quadrant' do
    result = cdf(
      [-Float::INFINITY, -Float::INFINITY],
      [0.0, 0.0],
      mean: [0.0, 0.0],
      covariance: [[1.0, 0.0], [0.0, 1.0]]
    )

    expect(result[:estimate]).to be_within(1e-12).of(0.25)
    expect(result[:error_kind]).to eq(:monte_carlo_confidence_interval_half_width)
    expect(result[:error_description]).to include('Statistical Monte Carlo error estimate')
  end

  it 'estimates a correlated zero-mean orthogonal quadrant near its analytic value' do
    rho = 0.7
    expected = 0.25 + Math.asin(rho) / (2.0 * Math::PI)
    result = cdf(
      [-Float::INFINITY, -Float::INFINITY],
      [0.0, 0.0],
      mean: [0.0, 0.0],
      covariance: [[1.0, rho], [rho, 1.0]],
      epsilon: 0.001,
      max_iterations: 20_000
    )

    expect(result[:estimate]).to be_within(0.02).of(expected)
    expect(result[:error_estimate]).to be > 0.0
    expect(result[:error_kind]).to eq(:monte_carlo_confidence_interval_half_width)
  end

  it 'accepts infinite bounds and shifts non-zero means' do
    result = cdf(
      [1.0, -Float::INFINITY],
      [Float::INFINITY, 1.0],
      mean: [1.0, 0.5],
      covariance: [[1.0, 0.0], [0.0, 4.0]],
      epsilon: 1e-9
    )
    expected = 0.5 * one_dimensional_probability(-Float::INFINITY, 1.0, 0.5, 4.0)

    expect(result[:estimate]).to be_within(1e-12).of(expected)
  end

  it 'is reproducible with the same explicit random source seed' do
    arguments = [
      [-Float::INFINITY, -Float::INFINITY],
      [0.0, 0.0],
      {
        mean: [0.0, 0.0],
        covariance: [[1.0, 0.7], [0.7, 1.0]],
        max_iterations: 500
      }
    ]
    first = described_class.cdf(*arguments[0, 2], **arguments[2], random: Random.new(123))
    second = described_class.cdf(*arguments[0, 2], **arguments[2], random: Random.new(123))

    expect(second).to eq(first)
    expect(second[:samples]).to eq(500)
    expect(second[:target_reached]).to eq(false)
  end

  it 'rejects non-positive-definite covariance' do
    expect do
      cdf([0.0], [1.0], mean: [0.0], covariance: [[0.0]])
    end.to raise_error(ArgumentError, 'covariance must be positive definite')

    expect do
      cdf(
        [0.0, 0.0],
        [1.0, 1.0],
        mean: [0.0, 0.0],
        covariance: [[1.0, 0.0], [0.0, -1.0]]
      )
    end.to raise_error(ArgumentError, 'covariance must be positive definite')
  end

  it 'rejects malformed dimensions and bounds' do
    expect do
      cdf([0.0], [1.0], mean: [0.0, 0.0], covariance: [[1.0, 0.0], [0.0, 1.0]])
    end.to raise_error(ArgumentError, 'lower and upper must match mean dimension')

    expect do
      cdf([2.0], [1.0], mean: [0.0], covariance: [[1.0]])
    end.to raise_error(ArgumentError, 'lower bound at dimension 0 must be <= upper bound')
  end

  it 'rejects non-symmetric covariance' do
    expect do
      cdf(
        [0.0, 0.0],
        [1.0, 1.0],
        mean: [0.0, 0.0],
        covariance: [[1.0, 0.1], [0.2, 1.0]]
      )
    end.to raise_error(ArgumentError, 'covariance must be symmetric')
  end

  it 'rejects an invalid explicit random source' do
    expect do
      cdf([0.0], [1.0], mean: [0.0], covariance: [[1.0]], random: Object.new)
    end.to raise_error(ArgumentError, 'random must respond to #rand')
  end

  it 'does not report a zero-variance rare-event estimate as a reached target' do
    random = Object.new
    def random.rand
      0.5
    end
    result = cdf(
      [5.0, 100.0],
      [Float::INFINITY, 101.0],
      mean: [0.0, 0.0],
      covariance: [[1.0, 0.7], [0.7, 1.0]],
      random: random,
      max_iterations: 100
    )

    expect(result[:estimate]).to eq(0.0)
    expect(result[:error_estimate]).to be > 0.0
    expect(result[:target_reached]).to eq(false)
  end

  it 'requires an explicit random source' do
    expect do
      described_class.cdf(
        [0.0],
        [1.0],
        mean: [0.0],
        covariance: [[1.0]]
      )
    end.to raise_error(ArgumentError)
  end
end
