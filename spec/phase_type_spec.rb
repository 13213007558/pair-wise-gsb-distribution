require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::PhaseType do
  let(:erlang) do
    described_class.new(
      [1.0, 0.0],
      [[-1.0, 1.0],
       [0.0, -1.0]]
    )
  end

  let(:branch) do
    described_class.new(
      [1.0, 0.0],
      [[-2.0, 0.5],
       [0.0, -1.0]]
    )
  end

  it 'represents two serial exponential stages as Erlang(2, 1)' do
    expect(erlang.initial_absorption_mass).to eq(0.0)
    expect(erlang.survival(0)).to eq(1.0)
    expect(erlang.density(0)).to eq(0.0)
    expect(erlang.cdf(0)).to eq(0.0)

    expect(erlang.survival(1)).to be_within(1e-12).of(2.0 / Math::E)
    expect(erlang.density(1)).to be_within(1e-12).of(1.0 / Math::E)
    expect(erlang.cdf(1)).to be_within(1e-12).of(1.0 - 2.0 / Math::E)
  end

  it 'calculates Erlang raw moments by hand' do
    expect(erlang.moments_upto(3)).to eq([2.0, 6.0, 24.0])
    expect(erlang.moment(1)).to eq(2.0)
    expect(erlang.mean).to eq(2.0)
    expect(erlang.variance).to eq(2.0)
  end

  it 'conditions an Erlang model on survival' do
    residual = erlang.residual_life(1)

    residual.initial_probabilities.each_with_index do |probability, index|
      expect(probability).to be_within(1e-12).of([0.5, 0.5][index])
    end
    expect(residual.initial_absorption_mass).to eq(0.0)
    expect(residual.mean).to be_within(1e-12).of(1.5)
  end

  it 'calculates the hand values for a branched two-stage model' do
    state_one = Math.exp(-2.0)
    state_two = (Math.exp(-1.0) - Math.exp(-2.0)) / 2.0
    survival = state_one + state_two

    expect(branch.survival(1)).to be_within(1e-12).of(survival)
    expect(branch.density(1)).to be_within(1e-12).of(1.5 * state_one + state_two)
    expect(branch.cdf(1)).to be_within(1e-12).of(1.0 - survival)
  end

  it 'calculates branch moments without normalizing the branching probability' do
    branch.moments_upto(3).each_with_index do |value, index|
      expect(value).to be_within(1e-12).of([0.75, 1.25, 3.375][index])
    end
    expect(branch.mean).to be_within(1e-12).of(0.75)
    expect(branch.variance).to be_within(1e-12).of(0.6875)
  end

  it 'returns the branched residual phase distribution' do
    residual = branch.residual_life(1)
    state_one = Math.exp(-2.0)
    state_two = (Math.exp(-1.0) - Math.exp(-2.0)) / 2.0
    survival = state_one + state_two

    expected = [state_one / survival, state_two / survival]
    residual.initial_probabilities.each_with_index do |probability, index|
      expect(probability).to be_within(1e-12).of(expected[index])
    end
    expect(residual.mean).to be_within(1e-12).of(0.8655292893150026)
  end

  it 'keeps zero-time absorption mass separate' do
    defective = described_class.new(
      [0.2, 0.0],
      [[-1.0, 1.0],
       [0.0, -1.0]]
    )

    expect(defective.initial_absorption_mass).to be_within(1e-12).of(0.8)
    expect(defective.survival(0)).to be_within(1e-12).of(0.2)
    expect(defective.cdf(0)).to be_within(1e-12).of(0.8)
    expect(defective.density(0)).to eq(0.0)
    expect(defective.mean).to be_within(1e-12).of(0.4)
  end

  it 'supports immediate absorption' do
    absorbed = described_class.new([0.0], [[-1.0]])

    expect(absorbed.initial_absorption_mass).to eq(1.0)
    expect(absorbed.survival(1)).to eq(0.0)
    expect(absorbed.density(1)).to eq(0.0)
    expect(absorbed.cdf(1)).to eq(1.0)
    expect(absorbed.moments_upto(2)).to eq([0.0, 0.0])
  end

  it 'does not cap density as if it were a probability' do
    exponential = described_class.new([1.0], [[-10.0]])

    expect(exponential.density(0)).to be_within(1e-12).of(10.0)
  end

  it 'rejects a phase that cannot reach absorption' do
    expect {
      described_class.new(
        [1.0, 0.0],
        [[-1.0, 1.0],
         [0.0, 0.0]]
      )
    }.to raise_error(ArgumentError, 'every phase must be able to reach absorption')
  end

  it 'rejects an initial vector that starts entirely in a closed phase' do
    expect {
      described_class.new(
        [0.0, 1.0],
        [[-1.0, 0.0],
         [0.0, 0.0]]
      )
    }.to raise_error(ArgumentError, 'every phase must be able to reach absorption')
  end

  it 'rejects non-subgenerator signs and rows' do
    expect {
      described_class.new([1.0], [[1.0]])
    }.to raise_error(ArgumentError, 'generator diagonal entries must be non-positive')

    expect {
      described_class.new(
        [1.0, 0.0],
        [[-1.0, -1.0],
         [0.0, -1.0]]
      )
    }.to raise_error(ArgumentError, 'generator off-diagonal entries must be non-negative')

    expect {
      described_class.new(
        [1.0, 0.0],
        [[0.0, 1.0],
         [0.0, -1.0]]
      )
    }.to raise_error(ArgumentError, 'generator rows must sum to zero or a negative absorption rate')
  end

  it 'does not normalize an initial vector whose total exceeds one' do
    expect {
      described_class.new([0.6, 0.6], [[-1.0, 0.0], [0.0, -1.0]])
    }.to raise_error(ArgumentError, 'initial transient probability cannot exceed one')
  end

  it 'rejects negative times and conditioning on an impossible survival time' do
    expect { erlang.survival(-1) }.to raise_error(ArgumentError, 'time must be non-negative')
    expect { erlang.density(-1) }.to raise_error(ArgumentError, 'time must be non-negative')
    expect { erlang.cdf(-1) }.to raise_error(ArgumentError, 'time must be non-negative')
    expect { described_class.new([0.0], [[-1.0]]).residual_life(1) }
      .to raise_error(RangeError, /cannot condition on survival/)
  end
end
