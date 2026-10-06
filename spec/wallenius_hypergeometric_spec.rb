require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::WalleniusNoncentralHypergeometric do
  let(:engine) { Distribution::WalleniusNoncentralHypergeometric::Ruby_ }

  it 'exposes the short distribution name' do
    expect(described_class::SHORTHAND).to eq('wnhypg')
  end

  it 'matches ordinary hypergeometric probabilities when both weights are equal' do
    [[3, 2, 4], [5, 7, 6], [1, 4, 3]].each do |m1, m2, n|
      lower = [0, n - m2].max
      upper = [m1, n].min
      (lower..upper).each do |k|
        ordinary = Distribution::Hypergeometric::Ruby_.pdf(k, m1, n, m1 + m2)
        expect(engine.pdf(k, m1, m2, n, 1)).to eq(ordinary)
        expect(engine.pmf(k, m1, m2, n, 1)).to eq(ordinary)
        expect(engine.exact_pdf(k, m1, m2, n, Rational(1))).to eq(ordinary)
      end
    end
  end

  it 'does not use Fisher noncentral hypergeometric weights' do
    m1 = 2
    m2 = 3
    n = 2
    odds = 2
    fisher_numerators = (0..2).map do |k|
      Math.binomial_coefficient(m1, k) *
        Math.binomial_coefficient(m2, n - k) * odds**k
    end
    fisher_pmf = fisher_numerators[1] / fisher_numerators.sum

    expect(engine.pdf(1, m1, m2, n, odds)).not_to be_within(1e-6).of(fisher_pmf)
    expect(engine.pdf(1, m1, m2, n, odds)).to be_within(1e-12).of(Rational(22, 35))
  end

  it 'enumerates every path for small populations' do
    cases = [[4, 3, 5, 2], [2, 5, 4, Rational(1, 3)], [6, 2, 5, 7]]
    cases.each do |m1, m2, n, odds|
      lower = [0, n - m2].max
      upper = [m1, n].min
      exact = {}

      (lower..upper).each do |k|
        exact[k] = engine.exact_pdf(k, m1, m2, n, odds)
        expect(engine.pdf(k, m1, m2, n, odds.to_f)).to be_within(1e-12).of(exact[k])
      end

      expect(exact.values.inject(:+)).to eq(1)
    end
  end

  it 'returns cumulative probabilities over the feasible support' do
    m1 = 5
    m2 = 4
    n = 6
    odds = 2.5
    lower = [0, n - m2].max
    upper = [m1, n].min
    cumulative = 0

    (lower..upper).each do |k|
      cumulative += engine.pdf(k, m1, m2, n, odds)
      expect(engine.cdf(k, m1, m2, n, odds)).to be_within(1e-12).of(cumulative)
    end

    expect(engine.cdf(lower - 1, m1, m2, n, odds)).to eq(0)
    expect(engine.cdf(upper, m1, m2, n, odds)).to eq(1)
  end

  it 'calculates the expectation from the probability mass' do
    m1 = 4
    m2 = 3
    n = 5
    odds = Rational(2)
    lower = [0, n - m2].max
    upper = [m1, n].min
    expected = (lower..upper).inject(0) do |total, k|
      total + k * engine.exact_pdf(k, m1, m2, n, odds)
    end

    expect(engine.mean(m1, m2, n, odds)).to eq(expected)
    expect(engine.mean(m1, m2, n, odds.to_f)).to be_within(1e-12).of(expected)
    expect(engine.mean(m1, m2, n, 2)).to be_within(1e-12).of(expected)
    expect(engine.expectation(m1, m2, n, odds)).to eq(expected)
    expect(engine.expected_value(m1, m2, n, odds)).to eq(expected)
  end

  it 'samples the full sequential process with reproducible randomness' do
    draws_a = engine.sequential_sample(6, 5, 8, 3, 99)
    draws_b = engine.sequential_sample(6, 5, 8, 3, 99)

    expect(draws_a).to eq(draws_b)
    expect(draws_a.length).to eq(8)
    expect(draws_a.all? { |draw| [0, 1].include?(draw) }).to be(true)
    expect(engine.rng(6, 5, 8, 3, 99).call).to eq(draws_a.count(1))
  end

  it 'keeps probabilities legal when a color is exhausted or the odds are extreme' do
    expect(engine.pdf(0, 2, 5, 5, 0)).to eq(1.0)
    expect(engine.pdf(2, 2, 5, 5, Float::INFINITY)).to eq(1.0)
    expect(engine.pdf(3, 3, 1, 4, Float::INFINITY)).to eq(1.0)
    expect(engine.pdf(3, 3, 1, 4, 0.0)).to eq(1.0)
    expect(engine.cdf(3, 3, 1, 4, Float::INFINITY)).to eq(1)
    [[1e300, 1], [1e-300, 0]].each do |odds, forced_first_draw|
      probabilities = (1..3).map { |k| engine.pdf(k, 3, 4, 5, odds) }
      expect(probabilities).to all(be_between(0, 1))
      expect(probabilities.sum).to be_within(1e-12).of(1)
      expect(engine.sequential_sample(3, 4, 5, odds, 1).first).to eq(forced_first_draw)
    end

    random = double('random')
    allow(random).to receive(:rand).and_return(1.0)
    expect(engine.sample(2, 1, 2, Float::INFINITY, random: random)).to eq(2)
  end

  it 'reports integral convergence status without changing the legal default PMF' do
    probability, status = engine.pdf_with_status(2, 4, 3, 5, 2.0)
    expect(status).to eq(:ok)
    expect(probability).to be_within(1e-10).of(engine.exact_pdf(2, 4, 3, 5, 2))

    _, tight_status = engine.pdf_with_status(2, 10, 7, 9, 0.125, 1e-15, 1)
    expect(tight_status).to eq(:did_not_converge)

    extreme_probability, extreme_status = engine.pdf_with_status(2, 3, 4, 5, 1e300)
    expect(extreme_status).to eq(:integral_not_evaluated)
    expect(extreme_probability).to be_between(0, 1)
  end

  it 'rejects invalid populations and odds' do
    expect { engine.pdf(0, -1, 2, 1) }.to raise_error(ArgumentError)
    expect { engine.pdf(0, 1, 1, 3) }.to raise_error(ArgumentError)
    expect { engine.pdf(0, 1, 2, 1, -1) }.to raise_error(ArgumentError)
    expect { engine.pdf(0, 1, 2, 1, Float::NAN) }.to raise_error(ArgumentError)
    expect { engine.pdf_with_status(0, 1, 2, 1, 1, -1) }.to raise_error(ArgumentError)
  end
end
