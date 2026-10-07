require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::TruncatedNormal do
  let(:normal_pdf) do
    lambda do |x, mu = 0.0, sigma = 1.0|
      z = (x - mu) / sigma
      Math.exp(-0.5 * z * z) / (Math.sqrt(2.0 * Math::PI) * sigma)
    end
  end

  def integral_pdf(distribution, subdivisions = 2000)
    lower = distribution.lower
    upper = distribution.upper
    return 1.0 unless lower.finite? && upper.finite?

    width = upper - lower
    step = width / subdivisions
    sum = distribution.pdf(lower) + distribution.pdf(upper)
    (1...subdivisions).each do |i|
      x = lower + i * step
      sum += (i.odd? ? 4.0 : 2.0) * distribution.pdf(x)
    end
    sum * step / 3.0
  end

  it 'matches the ordinary normal without truncation bounds' do
    distribution = described_class.new
    expect(distribution.mean).to eq(0.0)
    expect(distribution.variance).to eq(1.0)
    expect(distribution.pdf(0.75)).to be_within(1e-14).of(
      normal_pdf.call(0.75)
    )
    expect(distribution.cdf(0.75)).to be_within(1e-12).of(0.773372647623)
    expect(distribution.quantile(0.9)).to be_within(1e-10).of(1.2815515655)
  end

  it 'connects to one-sided normal truncation formulas' do
    lower = described_class.new(lower: 0.0)
    upper = described_class.new(upper: 0.0)
    lower_mean = Math.sqrt(2.0 / Math::PI)

    expect(lower.mean).to be_within(1e-12).of(lower_mean)
    expect(lower.variance).to be_within(1e-12).of(1.0 - 2.0 / Math::PI)
    expect(upper.mean).to be_within(1e-12).of(-lower_mean)
    expect(upper.variance).to be_within(1e-12).of(1.0 - 2.0 / Math::PI)
    expect(lower.pdf(1.0)).to be_within(1e-12).of(2.0 * normal_pdf.call(1.0))
    expect(upper.pdf(-1.0)).to be_within(1e-12).of(2.0 * normal_pdf.call(-1.0))
  end

  it 'normalizes and describes a very narrow far-tail interval' do
    a = 12.0
    b = 12.01
    distribution = described_class.new(lower: a, upper: b)
    denominator = 0.5 * (Math.erfc(a / Math.sqrt(2.0)) -
                          Math.erfc(b / Math.sqrt(2.0)))

    expect(distribution.mean).to be_between(a, b)
    expect(distribution.variance).to be_within(1.0e-9).of(8.3273e-6)
    expect(distribution.pdf(12.005)).to be_within(1.0e-8).of(
      normal_pdf.call(12.005) / denominator
    )
    expect(integral_pdf(distribution)).to be_within(1.0e-10).of(1.0)
    expect(distribution.quantile(0.5)).to be_between(a, b)
  end

  it 'supports narrow positive and negative far tails without saturation' do
    positive = described_class.new(lower: 100.0, upper: 100.01)
    negative = described_class.new(lower: -100.0, upper: -99.99)

    expect(positive.pdf(positive.mean)).to be_finite
    expect(negative.pdf(negative.mean)).to be_finite
    expect(positive.variance).to be_within(1.0e-8).of(7.933e-6)
    expect(negative.variance).to be_within(1.0e-8).of(7.933e-6)
    expect(negative.mean).to be_between(-100.0, -99.99)
  end

  it 'uses inverse-CDF sampling with an injected random source' do
    values = [0.0, 0.25, 0.5, 0.75, 1.0]
    index = 0
    random = Object.new
    random.define_singleton_method(:rand) do
      values[index].tap { index += 1 }
    end
    distribution = described_class.new(
      lower: 12.0, upper: 12.01, random: random
    )

    samples = Array.new(5) { distribution.sample }

    expect(distribution.call(0.25)).to be_between(12.0, 12.01)
    expect(samples.first).to eq(12.0)
    expect(samples.last).to eq(12.01)
    expect(samples).to all(be_between(12.0, 12.01))
    expect(distribution.attempts).to eq(5)
  end

  it 'accepts a callable random source and counts attempts' do
    values = [0.25, 0.75]
    distribution = described_class.new(
      lower: 0.0, random: proc { values.shift }
    )

    2.times { distribution.sample }
    expect(distribution.attempts).to eq(2)
  end

  it 'inverts bounded and unbounded CDFs accurately' do
    examples = [
      described_class.new,
      described_class.new(lower: 0.0),
      described_class.new(upper: 0.0),
      described_class.new(lower: 12.0, upper: 12.01),
      described_class.new(lower: -12.0, upper: -11.99),
      described_class.new(lower: 100.0, upper: 100.01),
      described_class.new(lower: -1.0, upper: 1.0),
      described_class.new(lower: -5.0, upper: 0.1)
    ]

    examples.each do |distribution|
      [0.001, 0.01, 0.5, 0.99, 0.999].each do |probability|
        quantile = distribution.quantile(probability)
        expect(distribution.cdf(quantile)).to be_within(1.0e-10).of(
          probability
        )
      end
    end
  end

  it 'scales and shifts with mu and sigma' do
    distribution = described_class.new(
      mu: 100.0, sigma: 2.0, lower: 124.0, upper: 124.002
    )

    expect(distribution.mean).to be_between(124.0, 124.002)
    expect(distribution.variance).to be_within(1.0e-9).of(3.3332e-7)
    expect(integral_pdf(distribution)).to be_within(1.0e-9).of(1.0)
  end

  it 'rejects invalid parameters' do
    expect { described_class.new(sigma: 0.0) }.to raise_error(ArgumentError)
    expect { described_class.new(lower: 1.0, upper: 1.0) }
      .to raise_error(ArgumentError)
    expect { described_class.new.quantile(1.1) }.to raise_error(
      ArgumentError
    )
  end
end
