require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::SaddlepointApproximation do
  def gamma_approximation(shape)
    described_class.new(
      kind: :continuous,
      cgf: ->(t) { fail described_class::DomainError if t >= 1.0; -shape * Math.log(1.0 - t) },
      cgf_derivative: ->(t) { fail described_class::DomainError if t >= 1.0; shape / (1.0 - t) },
      second_derivative: ->(t) { fail described_class::DomainError if t >= 1.0; shape / (1.0 - t)**2 },
      third_derivative: ->(t) { 2.0 * shape / (1.0 - t)**3 }
    )
  end

  def binomial_approximation(n, pr)
    described_class.new(
      kind: :lattice,
      lattice_spacing: 1.0,
      cgf: ->(t) { n * Math.log(1.0 - pr + pr * Math.exp(t)) },
      cgf_derivative: ->(t) do
        denominator = 1.0 - pr + pr * Math.exp(t)
        n * pr * Math.exp(t) / denominator
      end,
      second_derivative: ->(t) do
        denominator = 1.0 - pr + pr * Math.exp(t)
        n * pr * Math.exp(t) * (1.0 - pr) / denominator**2
      end,
      third_derivative: ->(t) do
        denominator = 1.0 - pr + pr * Math.exp(t)
        z = pr * Math.exp(t) / denominator
        n * z * (1.0 - z) * (1.0 - 2.0 * z)
      end
    )
  end

  it 'marks continuous results as approximations and provides root diagnostics' do
    saddle = gamma_approximation(3.0)
    result = saddle.upper_tail(5.0)

    expect(result.status).to eq(:ok)
    expect(result).to be_approximate
    expect(result.saddlepoint).to be_within(1.0e-8).of(0.4)
    expect(result.diagnostics[:approximation]).to eq(:lugannani_rice)
    expect(result.diagnostics[:iterations]).to be > 0
    expect(result.diagnostics[:residual].abs).to be < 1.0e-8
  end

  it 'uses the mean limit instead of evaluating a 0/0 Lugannani-Rice expression' do
    saddle = gamma_approximation(3.0)
    exact = 1.0 - Distribution::Gamma.cdf(3.0, 3.0, 1.0)
    result = saddle.upper_tail(3.0)

    expect(result.saddlepoint).to eq(0.0)
    expect(result.diagnostics[:approximation]).to eq(:lugannani_rice_mean_limit)
    expect(result.probability).to be_within(1.0e-4).of(exact)
    expect(saddle.cdf(3.0).probability).to be_within(1.0e-4).of(1.0 - exact)
  end

  it 'uses the first-order mean expansion for a small displacement from the mean' do
    saddle = gamma_approximation(3.0)
    exact = 1.0 - Distribution::Gamma.cdf(3.0001, 3.0, 1.0)
    result = saddle.upper_tail(3.0001)

    expect(result.diagnostics[:approximation]).to eq(:lugannani_rice_mean_expansion)
    expect(result.probability).to be < saddle.upper_tail(3.0).probability
    expect(result.probability).to be_within(1.0e-4).of(exact)
  end

  it 'approximates Gamma tails away from the mean' do
    saddle = gamma_approximation(3.0)

    {5.0 => 2.0e-4, 10.0 => 2.0e-5, 0.5 => 2.0e-4}.each do |x, tolerance|
      exact = 1.0 - Distribution::Gamma.cdf(x, 3.0, 1.0)
      expect(saddle.upper_tail(x).probability).to be_within(tolerance).of(exact)
      expect(saddle.cdf(x).probability).to be_within(tolerance).of(1.0 - exact)
    end
  end

  it 'does not force a root when a Gamma target is outside CGF-reachable support' do
    saddle = gamma_approximation(3.0)
    result = saddle.upper_tail(0.0)

    expect(result.status).to eq(:saddlepoint_not_found)
    expect(result.probability).to be_nil
    expect(result.approximate?).to eq(false)
    expect(result.diagnostics[:message]).to match(/saddlepoint|CGF/)
  end

  it 'approximates inclusive lattice tails and maps continuous thresholds to lattice points' do
    saddle = binomial_approximation(100.0, 0.5)
    result = saddle.upper_tail(60.0)
    exact = 1.0 - Distribution::Binomial.cdf(59.0, 100.0, 0.5)

    expect(result.lattice_point).to eq(60.0)
    expect(result.diagnostics[:approximation]).to eq(:daniels_lattice)
    expect(result.probability).to be_within(1.0e-4).of(exact)

    fractional = saddle.upper_tail(59.2)
    expect(fractional.lattice_point).to eq(60.0)
    expect(fractional.requested_x).to eq(59.2)
    expect(fractional.probability).to be_within(1.0e-4).of(exact)
  end

  it 'uses a lattice mean limit with the local half-mass correction' do
    saddle = binomial_approximation(100.0, 0.5)
    exact = 1.0 - Distribution::Binomial.cdf(49.0, 100.0, 0.5)
    result = saddle.upper_tail(50.0)

    expect(result.saddlepoint).to eq(0.0)
    expect(result.diagnostics[:approximation]).to eq(:daniels_mean_limit)
    expect(result.probability).to be_within(2.0e-4).of(exact)
  end

  it 'approximates binomial CDFs in both tails while respecting lattice atoms' do
    saddle = binomial_approximation(100.0, 0.5)

    expect(saddle.cdf(60.0).probability).to be_within(1.0e-4).of(
      Distribution::Binomial.cdf(60.0, 100.0, 0.5)
    )
    expect(saddle.cdf(40.0).probability).to be_within(1.0e-4).of(
      Distribution::Binomial.cdf(40.0, 100.0, 0.5)
    )
    expect(saddle.cdf(70.0).probability).to be_within(1.0e-6).of(
      Distribution::Binomial.cdf(70.0, 100.0, 0.5)
    )
  end

  it 'accepts a custom lattice spacing and works without a supplied third derivative' do
    n = 50.0
    saddle = described_class.new(
      kind: :lattice,
      lattice_spacing: 2.0,
      cgf: ->(t) { n * Math.log(1.0 + Math.exp(2.0 * t)) - n * Math.log(2.0) },
      cgf_derivative: ->(t) { 2.0 * n / (1.0 + Math.exp(-2.0 * t)) },
      second_derivative: ->(t) do
        probability = 1.0 / (1.0 + Math.exp(-2.0 * t))
        4.0 * n * probability * (1.0 - probability)
      end
    )

    result = saddle.upper_tail(80.0)
    expect(result.lattice_point).to eq(80.0)
    expect(result.status).to eq(:ok)
    expect(result.probability).to be_between(0.0, 1.0)
  end

  it 'rejects invalid construction parameters' do
    expect do
      described_class.new(kind: :unknown, cgf: proc {}, cgf_derivative: proc {},
                          second_derivative: proc {})
    end.to raise_error(ArgumentError)

    expect do
      described_class.new(kind: :lattice, lattice_spacing: 0,
                          cgf: proc {}, cgf_derivative: proc {},
                          second_derivative: proc {})
    end.to raise_error(ArgumentError)
  end
end
