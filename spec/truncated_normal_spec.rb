require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::TruncatedNormal do
  let(:narrow_tail) { described_class.new(12.0, 12.01) }

  it 'connects to the untruncated normal distribution' do
    distribution = described_class.new

    expect(distribution.mean).to eq(0.0)
    expect(distribution.variance).to eq(1.0)
    expect(distribution.pdf(0.0)).to be_within(1e-15).of(1.0 / Math.sqrt(2 * Math::PI))
    expect(distribution.pdf(5.0)).to be_within(1e-15).of(Distribution::Normal.pdf(5.0))
    expect(distribution.cdf(1.96)).to be_within(1e-12).of(Distribution::Normal.cdf(1.96))
    expect(distribution.quantile(0.975)).to be_within(1e-10).of(1.959963984540054)
  end

  it 'supports one-sided truncation' do
    lower = described_class.new(0.0, Float::INFINITY)
    upper = described_class.new(-Float::INFINITY, 0.0)

    expect(lower.mean).to be_within(1e-12).of(Math.sqrt(2.0 / Math::PI))
    expect(upper.mean).to be_within(1e-12).of(-Math.sqrt(2.0 / Math::PI))
    expect(lower.variance).to be_within(1e-12).of(1.0 - 2.0 / Math::PI)
    expect(upper.variance).to be_within(1e-12).of(1.0 - 2.0 / Math::PI)

    expect(lower.quantile(0.0)).to eq(0.0)
    expect(lower.quantile(0.5)).to be_within(1e-10).of(0.6744897501960817)
    expect(lower.quantile(1.0)).to eq(Float::INFINITY)
    expect(upper.quantile(0.5)).to be_within(1e-10).of(-0.6744897501960817)
  end

  it 'calculates normalized density and finite-interval moments' do
    distribution = described_class.new(-1.0, 1.0)

    expect(distribution.pdf(0.0)).to be_within(1e-12).of(0.5843685672568165)
    expect(distribution.cdf(0.0)).to be_within(1e-12).of(0.5)
    expect(distribution.mean).to be_within(1e-12).of(0.0)
    expect(distribution.variance).to be_within(1e-12).of(0.2911250947727933)
    expect(distribution.pdf(-1.0001)).to eq(0.0)
    expect(distribution.log_pdf(-1.0001)).to eq(-Float::INFINITY)
  end

  it 'handles very narrow distant tails without CDF saturation' do
    expect(narrow_tail.pdf(12.005)).to be_within(1e-8).of(99.94039183873409)
    expect(narrow_tail.mean).to be_between(12.0, 12.01)
    expect(narrow_tail.mean).to be_within(2e-8).of(12.0049)
    expect(narrow_tail.variance).to be_within(1e-8).of((0.01**2) / 12.0)
    expect(narrow_tail.quantile(0.5)).to be_within(5e-8).of(12.00485)
    expect(narrow_tail.quantile(0.001)).to be > 12.0
    expect(narrow_tail.quantile(0.999)).to be < 12.01
  end

  it 'handles narrow tails below the mean' do
    distribution = described_class.new(-12.01, -12.0)

    expect(distribution.mean).to be_within(2e-8).of(-12.0049)
    expect(distribution.variance).to be_within(1e-8).of((0.01**2) / 12.0)
    expect(distribution.quantile(0.5)).to be_within(5e-8).of(-12.00485)
  end

  it 'supports locations and scales other than the standard normal' do
    distribution = described_class.new(2.0, 2.0001, mean: 5.0, stddev: 2.0)

    expect(distribution.normal_mean).to eq(5.0)
    expect(distribution.normal_stddev).to eq(2.0)
    expect(distribution.mean).to be_within(1e-9).of(2.00005)
    expect(distribution.variance).to be_within(1e-12).of(4.0 * (0.00005**2) / 12.0)
    expect(distribution.cdf(distribution.quantile(0.375))).to be_within(1e-9).of(0.375)
  end

  [
    [-Float::INFINITY, Float::INFINITY],
    [0.0, Float::INFINITY],
    [-Float::INFINITY, 0.0],
    [-1.0, 1.0],
    [12.0, 12.01],
    [-12.01, -12.0],
    [10.0, 100.0],
    [-37.0, -36.999]
  ].each do |lower, upper|
    it "inverts CDF values for #{[lower, upper].inspect}" do
      distribution = described_class.new(lower, upper)

      [0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999].each do |probability|
        value = distribution.quantile(probability)
        expect(value).to be_between(lower, upper)
        expect(distribution.cdf(value)).to be_within(1e-8).of(probability)
      end
    end
  end

  it 'samples once through inverse transformation and records attempts' do
    probabilities = [0.0, 0.25, 0.5, 0.75, 1.0]
    values = []
    distribution = described_class.new(12.0, 12.01, rng: proc { probabilities[values.size] })

    5.times { values << distribution.sample }

    expect(values).to eq(probabilities.map { |probability| distribution.quantile(probability) })
    expect(values).to eq(values.sort)
    expect(distribution.attempts).to eq(5)
    distribution.reset_attempts
    expect(distribution.attempts).to eq(0)
  end

  it 'accepts a Random-like random source and an explicit seed' do
    random = Random.new(1234)
    first = described_class.new(-1.0, 1.0, rng: random).sample
    seeded_a = described_class.new(-1.0, 1.0, seed: 42).sample
    seeded_b = described_class.new(-1.0, 1.0, seed: 42).sample

    expect(first).to be_between(-1.0, 1.0)
    expect(seeded_a).to eq(seeded_b)
  end

  it 'rejects invalid parameters' do
    expect { described_class.new(1.0, 1.0) }.to raise_error(ArgumentError)
    expect { described_class.new(1.0, 0.0) }.to raise_error(ArgumentError)
    expect { described_class.new(0.0, 1.0, stddev: 0.0) }.to raise_error(ArgumentError)
    expect { described_class.new(0.0, 1.0, mean: Float::INFINITY) }
      .to raise_error(ArgumentError)
    expect { described_class.new(Float::INFINITY, 1.0) }.to raise_error(ArgumentError)
    expect { narrow_tail.quantile(1.1) }.to raise_error(RangeError)
  end

  it 'rejects intervals whose probability underflows in Float' do
    expect { described_class.new(100.0, 100.0001) }
      .to raise_error(RangeError, /zero probability/)
  end
end
