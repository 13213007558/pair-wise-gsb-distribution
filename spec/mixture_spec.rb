require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')
require 'distribution/uniform'

describe Distribution::Mixture do
  let(:uniform) { Distribution::Uniform }

  def mixed_distribution
    described_class.new(
      [
        { distribution: uniform, args: [1, 2], weight: 2 },
        { distribution: uniform, args: [3, 4], weight: 2 }
      ],
      point_masses: { 0 => 4 }
    )
  end

  it 'normalizes component and explicit point-mass weights' do
    mixture = mixed_distribution

    expect(mixture.total_weight).to eq(8)
    expect(mixture.point_masses[0]).to be_within(1e-12).of(0.5)
  end

  it 'keeps continuous density separate from point masses' do
    mixture = mixed_distribution

    expect(mixture.density(0)).to eq(0.0)
    expect(mixture.density(1.5)).to be_within(1e-12).of(0.25)
    expect(mixture.density(2.5)).to eq(0.0)
    expect(mixture.density(3.5)).to be_within(1e-12).of(0.25)
    expect(mixture.pdf(1.5)).to be(mixture.density(1.5))
  end

  it 'reports atoms explicitly' do
    mixture = mixed_distribution

    expect(mixture.mass_at(0)).to be_within(1e-12).of(0.5)
    expect(mixture.mass_at(1)).to eq(0.0)
    expect(mixture.mass_at(3)).to eq(0.0)
  end

  it 'builds CDF jumps and plateaus explicitly' do
    mixture = mixed_distribution

    expect(mixture.cdf(-1)).to eq(0.0)
    expect(mixture.cdf(0)).to be_within(1e-12).of(0.5)
    expect(mixture.cdf(1)).to be_within(1e-12).of(0.5)
    expect(mixture.cdf(2)).to be_within(1e-12).of(0.75)
    expect(mixture.cdf(2.5)).to be_within(1e-12).of(0.75)
    expect(mixture.cdf(3)).to be_within(1e-12).of(0.75)
    expect(mixture.cdf(4)).to be_within(1e-12).of(1.0)
  end

  it 'returns generalized inverse quantiles at jumps and support gaps' do
    mixture = mixed_distribution

    expect(mixture.quantile(0)).to eq(0)
    expect(mixture.quantile(0.25)).to eq(0)
    expect(mixture.quantile(0.5)).to eq(0)
    expect(mixture.quantile(0.5000001)).to be_within(1e-10).of(1.0000004)
    expect(mixture.quantile(0.6)).to be_within(1e-10).of(1.4)
    expect(mixture.quantile(0.75)).to eq(2)
    expect(mixture.quantile(0.76)).to be_within(1e-10).of(3.04)
    expect(mixture.quantile(1)).to eq(4)
    expect(mixture.p_value(0.6)).to eq(mixture.quantile(0.6))
  end

  it 'merges explicit atoms at the same location before normalizing' do
    mixture = described_class.new(
      [[Distribution::PointMass.new(0), 2]],
      point_masses: { 0 => 3 }
    )

    expect(mixture.mass_at(0)).to eq(1.0)
    expect(mixture.cdf(0)).to eq(1.0)
  end

  it 'merges numerically equivalent atom locations' do
    mixture = described_class.new([0, [0.0, 2]])

    expect(mixture.point_masses.keys).to contain_exactly(0)
    expect(mixture.mass_at(0.0)).to eq(1.0)
  end

  it 'rejects negative component weights' do
    expect do
      described_class.new([{ distribution: uniform, args: [1, 2], weight: -1 }])
    end.to raise_error(ArgumentError, /non-negative/)
  end

  it 'rejects negative point-mass weights' do
    expect do
      described_class.new([], point_masses: { 0 => -0.1 })
    end.to raise_error(ArgumentError, /non-negative/)
  end

  it 'preserves weights that already sum to one when normalization is disabled' do
    mixture = described_class.new(
      [[uniform, 0.75, [1, 2]]],
      point_masses: { 0 => 0.25 },
      normalize: false
    )

    expect(mixture.mass_at(0)).to be_within(1e-12).of(0.25)
    expect(mixture.density(1.5)).to be_within(1e-12).of(0.75)

    expect do
      described_class.new(
        [[uniform, 1.5, [1, 2]]],
        point_masses: { 0 => 0.5 },
        normalize: false
      )
    end.to raise_error(ArgumentError, /sum to 1/)
  end

  it 'rejects probabilities outside [0, 1]' do
    mixture = mixed_distribution

    expect { mixture.quantile(-0.1) }.to raise_error(RangeError)
    expect { mixture.quantile(1.1) }.to raise_error(RangeError)
  end

  it 'keeps generalized inverse endpoints infinite for unbounded components' do
    mixture = described_class.new(
      [[Distribution::Normal, 0.5]],
      point_masses: { 0 => 0.5 }
    )

    expect(mixture.quantile(0)).to eq(-Float::INFINITY)
  end
end
