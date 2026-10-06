require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::PhaseType do
  let(:rate) { 2.0 }
  let(:stages) { 3 }
  let(:erlang_t) do
    Array.new(stages) do |i|
      Array.new(stages) do |j|
        if i == j
          -rate
        elsif j == i + 1
          rate
        else
          0.0
        end
      end
    end
  end
  let(:erlang) { Distribution::PhaseType::Ruby_.new([1.0, 0.0, 0.0], erlang_t) }

  def erlang_pdf(x)
    rate**stages * x**(stages - 1) * Math.exp(-rate * x) /
      (1...stages).inject(1, :*)
  end

  def erlang_cdf(x)
    tail = (0...stages).inject(0.0) do |sum, k|
      sum + (rate * x)**k / (1..k).inject(1, :*)
    end
    1 - Math.exp(-rate * x) * tail
  end

  it 'builds a pure Ruby phase-type distribution' do
    expect(erlang).to be_a(Distribution::PhaseType::Ruby_::Distribution)
    expect(erlang.eventual_absorption_probability).to eq(1.0)
    expect(erlang.infinite_mass).to eq(0.0)
    expect(erlang.reachable_states).to eq([0, 1, 2])
  end

  it 'matches Erlang pdf, cdf and survival values' do
    [0.0, 0.25, 1.0, 2.5, 7.0].each do |x|
      expect(erlang.pdf(x)).to be_within(1e-11).of(erlang_pdf(x))
      expect(erlang.cdf(x)).to be_within(1e-11).of(erlang_cdf(x))
      expect(erlang.sf(x)).to be_within(1e-11).of(1 - erlang_cdf(x))
    end
  end

  it 'returns Erlang moments from stable linear solves' do
    expect(erlang.moment(1)).to be_within(1e-12).of(1.5)
    expect(erlang.moment(2)).to be_within(1e-12).of(3.0)
    expect(erlang.moment(3)).to be_within(1e-12).of(7.5)
    expect(erlang.mean).to be_within(1e-12).of(1.5)
  end

  it 'inverts finite Erlang quantiles and reports convergence information' do
    x = erlang.quantile(0.7)
    expect(erlang.cdf(x)).to be_within(1e-9).of(0.7)

    value, info = erlang.quantile(0.7, with_info: true)
    expect(value).to eq(x)
    expect(info[:converged]).to eq(true)
    expect(info[:iterations]).to be > 0
    expect(info[:residual]).to be <= 1e-9
  end

  it 'returns numerical evolution information' do
    value, info = erlang.cdf(1.0, with_info: true)
    expect(value).to be_within(1e-11).of(erlang_cdf(1.0))
    expect(info[:error]).to be <= 1e-11
    expect([:uniformization, :scaled_matrix_exponential]).to include(info[:method])
    expect(info[:iterations]).to be > 0
  end

  it 'stays stable for near-degenerate repeated rates' do
    near_t = [[-1.0, 1.0], [0.0, -(1.0 + 1e-10)]]
    near = Distribution::PhaseType.new([1.0, 0.0], near_t)
    delta = 1e-10
    mean = 1.0 + 1.0 / (1.0 + delta)
    second = 2.0 + 2.0 / (1.0 + delta)**2 + 2.0 / (1.0 + delta)
    expect(near.moment(1)).to be_within(1e-10).of(mean)
    expect(near.moment(2)).to be_within(1e-9).of(second)
    expect(near.cdf(0.5)).to be_within(1e-10).of(0.09020401043863158)
  end

  it 'handles very large times without eigen decomposition' do
    value, info = erlang.cdf(10_000_000.0, with_info: true)
    expect(value).to be_within(1e-10).of(1.0)
    expect(info[:error]).to be <= 1e-10
    expect(info[:method]).to eq(:scaled_matrix_exponential)
  end

  it 'recognizes mass initially in the absorbing state' do
    ph = Distribution::PhaseType.new([0.25], [[-2.0]],
                                    initial_absorbed_probability: 0.75)
    expect(ph.cdf(0)).to eq(0.75)
    expect(ph.sf(0)).to eq(0.25)
    expect(ph.pdf(0)).to eq(0.0)
    expect(ph.quantile(0.75)).to eq(0.0)
    expect(ph.cdf(1)).to be_within(1e-12).of(0.75 + 0.25 * (1 - Math.exp(-2)))
  end

  it 'identifies probability trapped forever by closed transient components' do
    ph = Distribution::PhaseType::Ruby_.new([0.5, 0.5], [[-1.0, 1.0], [0.0, 0.0]])
    expect(ph.reachable_states).to eq([])
    expect(ph.closed_states).to eq([0, 1])
    expect(ph.eventual_absorption_probability).to eq(0.0)
    expect(ph.infinite_mass).to eq(1.0)
    expect(ph.cdf(1_000)).to eq(0.0)
    expect(ph.sf(1_000)).to eq(1.0)
    expect(ph.sf(Float::INFINITY)).to eq(1.0)
    expect(ph.pdf(2)).to eq(0.0)
    expect(ph.finite_moment(1)).to eq(0.0)
    expect(ph.mean).to eq(Float::INFINITY)
  end

  it 'returns zero at zero and infinity beyond final absorption probability' do
    ph = Distribution::PhaseType::Ruby_.new([0.5, 0.5], [[-1.0, 1.0], [0.0, 0.0]])
    expect(ph.quantile(0.0)).to eq(0.0)
    value, info = ph.quantile(0.5, with_info: true)
    expect(value).to eq(Float::INFINITY)
    expect(info[:reason]).to eq(:at_or_beyond_eventual_absorption_probability)
  end

  it 'handles mixed initial absorbed, finite transient and trapped mass' do
    ph = Distribution::PhaseType.new([0.25, 0.25], [[-1.0, 0.0], [0.0, 0.0]],
                                    initial_absorbed_probability: 0.5)
    expect(ph.eventual_absorption_probability).to be_within(1e-12).of(0.75)
    expect(ph.infinite_mass).to be_within(1e-12).of(0.25)
    expect(ph.cdf(0)).to eq(0.5)
    expect(ph.conditional_moment(1)).to be_within(1e-12).of(ph.finite_moment(1) / 0.75)
    expect(ph.cdf(Float::INFINITY)).to be_within(1e-12).of(0.75)
    expect(ph.sf(Float::INFINITY)).to be_within(1e-12).of(0.25)
    expect(ph.quantile(0.75)).to eq(Float::INFINITY)
  end

  it 'keeps mass that branches from a transient state into a closed component' do
    ph = Distribution::PhaseType.new([1.0, 0.0], [[-2.0, 1.0], [0.0, 0.0]])
    expect(ph.reachable_states).to eq([0])
    expect(ph.closed_states).to eq([1])
    expect(ph.absorption_probabilities).to eq([0.5])
    expect(ph.eventual_absorption_probability).to eq(0.5)
    expect(ph.infinite_mass).to eq(0.5)
    expect(ph.cdf(1)).to be_within(1e-12).of(0.5 * (1 - Math.exp(-2)))
    expect(ph.sf(1)).to be_within(1e-12).of(0.5 + 0.5 * Math.exp(-2))
    expect(ph.pdf(1)).to be_within(1e-12).of(Math.exp(-2))
    expect(ph.moment(1)).to be_within(1e-12).of(0.25)
    expect(ph.mean).to eq(Float::INFINITY)
    expect(ph.sf(Float::INFINITY)).to eq(0.5)
  end

  it 'rejects invalid off-diagonal entries' do
    expect {
      Distribution::PhaseType.new([1.0, 0.0], [[-1.0, -0.1], [0.0, -1.0]])
    }.to raise_error(Distribution::PhaseType::Error, /off-diagonal/)
  end

  it 'rejects positive generator rows' do
    expect {
      Distribution::PhaseType.new([1.0, 0.0], [[-1.0, 2.0], [0.0, -1.0]])
    }.to raise_error(Distribution::PhaseType::Error, /row sum/)
  end

  it 'rejects dimension and probability vector errors' do
    expect {
      Distribution::PhaseType.new([1.0, 0.0], [[-1.0]])
    }.to raise_error(Distribution::PhaseType::Error, /same dimension/)
    expect {
      Distribution::PhaseType.new([1.2], [[-1.0]])
    }.to raise_error(Distribution::PhaseType::Error, /0 and 1/)
  end
end
