require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::VonMises do
  let(:vm) { described_class }
  let(:two_pi) { 2.0 * Math::PI }

  it 'treats angles separated by one turn as the same direction' do
    [0.0, 1.0, -3.7, 25.0].each do |angle|
      expect(vm.pdf(angle + two_pi, 0, 4)).to be_within(1e-14).of(vm.pdf(angle, 0, 4))
      expect(vm.pdf(angle, two_pi + 1, 4)).to be_within(1e-14).of(vm.pdf(angle, 1, 4))
    end
  end

  it 'has the uniform density at zero concentration' do
    expected = 1.0 / two_pi
    [0, Math::PI / 3, Math::PI, two_pi - 0.1, two_pi].each do |angle|
      expect(vm.pdf(angle, 2, 0)).to be_within(1e-15).of(expected)
    end
  end

  it 'exposes density as an alias for pdf' do
    expect(vm.density(0.4, 0.2, 7)).to eq(vm.pdf(0.4, 0.2, 7))
  end

  it 'uses arc length for uniform probability around the cut' do
    start = -Math::PI / 180
    length = 2.0 * Math::PI / 180

    expect(vm.arc_probability(start, length, 1, 0)).
      to be_within(1e-14).of(1.0 / 180)
  end

  it 'keeps zero-length and full-turn arcs distinct' do
    expect(vm.arc_probability(0, 0, 0, 5)).to eq(0)
    expect(vm.arc_probability(0, two_pi, 0, 5)).to eq(1)
    expect(vm.arc_probability(0.3, -two_pi, 0.2, 5)).to eq(1)
  end

  it 'computes short arcs crossing the cut by direction and length' do
    start = two_pi - 0.02
    probability = vm.arc_probability(start, 0.04, 0, 10)

    expect(probability).to be_within(1e-12).of(0.0497675830347301)
    expect(vm.arc_probability(start + 4 * two_pi, 0.04, 2 * two_pi, 10)).
      to be_within(1e-13).of(probability)
  end

  it 'uses the opposite short arc for a negative signed length' do
    probability = vm.arc_probability(two_pi - 0.02, -0.04, 0, 10)

    expect(probability).to be_within(1e-12).of(0.0493716356911301)
  end

  it 'does not overflow for a very narrow high-concentration peak' do
    kappa = 1e100
    expected = Math.sqrt(kappa / two_pi)

    tolerance = expected * 1e-12
    expect(vm.pdf(0, 0, kappa)).to be_within(tolerance).of(expected)
    expect(vm.pdf(2 * two_pi, two_pi, kappa)).to be_within(tolerance).of(expected)
    expect(vm.pdf(Math::PI, 0, kappa)).to eq(0)
  end

  it 'computes a narrow high-concentration arc probability' do
    probability = vm.arc_probability(-0.01, 0.02, 0, 10_000)

    expect(probability).to be_within(1e-8).of(0.682681426066465)
    expect(vm.arc_probability(4 * two_pi - 0.01, 0.02, 2 * two_pi, 10_000)).
      to be_within(1e-12).of(probability)
  end

  it 'builds a monotone CDF on the caller-selected cut' do
    cut = 0.5
    values = [0.5, 0.7, 0.9, 1.1, 3.0, Math::PI, 6.0, cut + two_pi]
    cdfs = values.map { |value| vm.cdf(value, 1, 20, cut) }

    expect(cdfs.first).to eq(0)
    expect(cdfs.last).to eq(1)
    expect(cdfs).to eq(cdfs.sort)
  end

  it 'puts half the probability at each side of the mean when cut is the mean' do
    expect(vm.cdf(1, 1, 20, 1)).to eq(0)
    expect(vm.cdf(1 + Math::PI, 1, 20, 1)).to be_within(1e-12).of(0.5)
    expect(vm.cdf(1 + two_pi, 1, 20, 1)).to eq(1)
  end

  it 'inverts the CDF within one selected turn' do
    cut = 0.5
    [0, 0.01, 0.25, 0.5, 0.75, 0.99, 1].each do |probability|
      x = vm.quantile(probability, 1, 50, cut)

      expect(x).to be >= cut
      expect(x).to be <= cut + two_pi
      next if probability == 0 || probability == 1

      expect(vm.cdf(x, 1, 50, cut)).to be_within(1e-10).of(probability)
    end
  end

  it 'returns uniform quantiles for zero concentration' do
    cut = -1.0
    expect(vm.quantile(0, 0, 0, cut)).to eq(cut)
    expect(vm.quantile(0.25, 0, 0, cut)).to be_within(1e-15).of(cut + two_pi / 4)
    expect(vm.quantile(1, 0, 0, cut)).to be_within(1e-15).of(cut + two_pi)
  end

  it 'preserves all corresponding probabilities after integer-turn shifts' do
    shifts = [-two_pi, two_pi, 4 * two_pi]
    shifts.each do |shift|
      expect(vm.pdf(0.3 + shift, shift, 25)).to be_within(1e-14).of(vm.pdf(0.3, 0, 25))
      expect(vm.arc_probability(-0.2 + shift, 0.5, shift, 25)).
        to be_within(1e-13).of(vm.arc_probability(-0.2, 0.5, 0, 25))
      expect(vm.cdf(0.7 + shift, shift, 25, shift)).
        to be_within(1e-13).of(vm.cdf(0.7, 0, 25, 0))
    end
  end

  it 'rejects negative concentration and non-finite arguments' do
    expect { vm.pdf(0, 0, -1) }.to raise_error(RangeError)
    expect { vm.pdf(Float::NAN, 0, 1) }.to raise_error(ArgumentError)
    expect { vm.quantile(1.2, 0, 1) }.to raise_error(RangeError)
  end

  it 'restricts CDF x values to the selected one-turn interval' do
    expect { vm.cdf(-0.01, 0, 1, 0) }.to raise_error(RangeError)
    expect { vm.cdf(two_pi + 0.01, 0, 1, 0) }.to raise_error(RangeError)
  end
end
