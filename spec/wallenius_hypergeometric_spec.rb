require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::WalleniusHypergeometric do
  before do
    @engine = Distribution::WalleniusHypergeometric::Ruby_
  end

  # Ground truth: enumerate every sequential draw path, where each
  # path's probability is the product of per-draw class
  # probabilities proportional to remaining_balls * weight.
  def enumerate_pmf(m, n, total, odds)
    probs = Hash.new(0.0)
    walk = lambda do |r1, r2, drawn1, left, prob|
      if left.zero?
        probs[drawn1] += prob
      else
        w1 = odds * r1
        w2 = r2.to_f
        den = w1 + w2
        walk.call(r1 - 1, r2, drawn1 + 1, left - 1, prob * w1 / den) if r1.positive?
        walk.call(r1, r2 - 1, drawn1, left - 1, prob * w2 / den) if r2.positive?
      end
    end
    walk.call(m, total - m, 0, n, 1.0)
    probs
  end

  def support(m, n, total)
    [[0, n - (total - m)].max, [n, m].min]
  end

  shared_examples_for 'wallenius engine' do
    it 'matches the enumeration of all draw paths on small populations' do
      [[3, 4, 7, 2.5], [2, 3, 5, 0.4], [4, 5, 6, 1.7], [3, 3, 8, 0.25],
       [5, 6, 9, 3.0], [1, 1, 2, 7.0]].each do |m, n, total, odds|
        truth = enumerate_pmf(m, n, total, odds)
        lo, hi = support(m, n, total)
        (lo..hi).each do |k|
          expect(@engine.pdf(k, m, n, total, odds)).to be_within(1e-8).of(truth[k])
        end
      end
    end

    it 'reduces to the central hypergeometric distribution when odds == 1' do
      central = Distribution::Hypergeometric::Ruby_
      [[5, 4, 12], [3, 3, 8], [8, 5, 20]].each do |m, n, total|
        lo, hi = support(m, n, total)
        (lo..hi).each do |k|
          expect(@engine.pdf(k, m, n, total, 1.0)).to be_within(1e-12).of(central.pdf(k, m, n, total).to_f)
        end
      end
    end

    it 'returns valid probabilities that sum to 1, even for extreme odds' do
      [1e-8, 1e-4, 0.1, 0.5, 2.0, 1e4, 1e8].each do |odds|
        m = 4
        n = 4
        total = 9
        lo, hi = support(m, n, total)
        probs = (lo..hi).map { |k| @engine.pdf(k, m, n, total, odds) }
        probs.each do |p|
          expect(p).to be >= 0.0
          expect(p).to be <= 1.0 + 1e-12
        end
        expect(probs.inject(:+)).to be_within(1e-6).of(1.0)
      end
    end

    it 'keeps probabilities valid when a class is exhausted by the sample' do
      # n > total - m forces at least n - (total - m) interesting balls
      m = 4
      n = 5
      total = 6
      lo, hi = support(m, n, total)
      expect(lo).to eq(3)
      expect(@engine.pdf(lo - 1, m, n, total, 2.0)).to eq(0.0)
      sum = (lo..hi).sum { |k| @engine.pdf(k, m, n, total, 2.0) }
      expect(sum).to be_within(1e-10).of(1.0)
      # degenerate odds put all mass on the forced extreme
      expect(@engine.pdf(lo, m, n, total, 0.0)).to eq(1.0)
      expect(@engine.pdf(hi, m, n, total, 0.0)).to eq(0.0)
      expect(@engine.pdf(hi, m, n, total, Float::INFINITY)).to eq(1.0)
    end

    it 'handles degenerate urns and samples' do
      expect(@engine.pdf(0, 0, 3, 5, 2.0)).to eq(1.0)
      expect(@engine.pdf(0, 5, 0, 5, 2.0)).to eq(1.0)
      expect(@engine.pdf(3, 5, 3, 5, 2.0)).to eq(1.0)
      expect(@engine.pdf(5, 5, 5, 5, 2.0)).to eq(1.0)
    end

    it 'cdf matches the sum of pmf and reaches 1 at the support end' do
      m = 5
      n = 4
      total = 12
      odds = 3.0
      lo, hi = support(m, n, total)
      ac = 0.0
      (lo..hi).each do |k|
        ac += @engine.pdf(k, m, n, total, odds)
        expect(@engine.cdf(k, m, n, total, odds)).to be_within(1e-8).of(ac)
      end
      expect(@engine.cdf(hi, m, n, total, odds)).to be_within(1e-9).of(1.0)
      expect(@engine.cdf(lo - 1, m, n, total, odds)).to eq(0.0)
    end

    it 'cdf with odds == 1 matches the central hypergeometric cdf' do
      central = Distribution::Hypergeometric::Ruby_
      (0..4).each do |k|
        expect(@engine.cdf(k, 5, 4, 12, 1.0)).to be_within(1e-12).of(central.cdf(k, 5, 4, 12).to_f)
      end
    end

    it 'mean equals n*m/total for odds == 1 and grows with the odds' do
      m = 5
      n = 4
      total = 12
      expect(@engine.mean(m, n, total, 1.0)).to be_within(1e-10).of(n * m.quo(total))
      means = [0.1, 0.5, 1.0, 2.0, 10.0].map { |o| @engine.mean(m, n, total, o) }
      means.each_cons(2) { |a, b| expect(b).to be > a }
      lo, hi = support(m, n, total)
      means.each { |mu| expect(mu).to be_between(lo, hi).inclusive }
    end

    it 'mean matches the enumeration of all draw paths' do
      m = 3
      n = 4
      total = 7
      odds = 2.5
      truth = enumerate_pmf(m, n, total, odds)
      expected = truth.sum { |k, p| k * p }
      expect(@engine.mean(m, n, total, odds)).to be_within(1e-8).of(expected)
    end

    it 'reports :ok status when the integral converges' do
      r = @engine.pdf_with_status(2, 5, 4, 12, 3.0)
      expect(r.status).to eq(:ok)
      expect(r).to be_success
      expect(r.value).to be_within(1e-12).of(@engine.pdf(2, 5, 4, 12, 3.0))
      expect(@engine.cdf_with_status(2, 5, 4, 12, 3.0).status).to eq(:ok)
      expect(@engine.mean_with_status(5, 4, 12, 3.0).status).to eq(:ok)
    end

    it 'reports :not_converged when the integration cannot converge' do
      r = @engine.pdf_with_status(2, 5, 4, 12, 3.0, tol: 1e-30, max_depth: 2)
      expect(r.status).to eq(:not_converged)
      expect(r).not_to be_success
      expect(@engine.cdf_with_status(2, 5, 4, 12, 3.0, tol: 1e-30, max_depth: 2).status).to eq(:not_converged)
      expect(@engine.mean_with_status(5, 4, 12, 3.0, tol: 1e-30, max_depth: 2).status).to eq(:not_converged)
    end

    it 'samples sequentially within the support and reproduces the mean' do
      m = 4
      n = 4
      total = 9
      odds = 3.0
      lo, hi = support(m, n, total)
      samples = Array.new(4000) { |i| @engine.rng(m, n, total, odds, 1234 + i) }
      expect(samples.min).to be >= lo
      expect(samples.max).to be <= hi
      sample_mean = samples.sum.quo(samples.size)
      expect(sample_mean).to be_within(0.1).of(@engine.mean(m, n, total, odds))
    end

    it 'rng is deterministic for a fixed seed and handles exhausted classes' do
      expect(@engine.rng(4, 4, 9, 2.0, 42)).to eq(@engine.rng(4, 4, 9, 2.0, 42))
      # n equals the whole population: everything is drawn
      expect(@engine.rng(4, 9, 9, 2.0, 42)).to eq(4)
      # sample larger than the minority class: draws continue after exhaustion
      expect(@engine.rng(1, 4, 5, 3.0, 42)).to be_between(0, 1).inclusive
      # degenerate odds
      expect(@engine.rng(4, 4, 9, 0.0, 42)).to eq(0)
      expect(@engine.rng(4, 4, 9, Float::INFINITY, 42)).to eq(4)
    end

    it 'validates its parameters' do
      expect { @engine.pdf(1, 6, 2, 5, 1.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 3, 6, 5, 1.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 3, 2, 5, -1.0) }.to raise_error(ArgumentError)
      expect { @engine.mean(3, 2, 5, -0.5) }.to raise_error(ArgumentError)
      expect { @engine.rng(3, 2, 5, -2.0) }.to raise_error(ArgumentError)
    end
  end

  describe Distribution::WalleniusHypergeometric::Ruby_ do
    it_should_behave_like 'wallenius engine'
  end

  it 'exposes the api through the module' do
    expect(Distribution::WalleniusHypergeometric.pdf(2, 5, 4, 12, 3.0)).to be_within(1e-12).of(@engine.pdf(2, 5, 4, 12, 3.0))
    expect(Distribution::WalleniusHypergeometric.cdf(2, 5, 4, 12, 3.0)).to be_within(1e-12).of(@engine.cdf(2, 5, 4, 12, 3.0))
    expect(Distribution::WalleniusHypergeometric.mean(5, 4, 12, 3.0)).to be_within(1e-12).of(@engine.mean(5, 4, 12, 3.0))
  end

  describe Distribution::Shorthand do
    include Distribution::Shorthand
    it 'exposes the whypg shortcuts' do
      expect(Distribution::Shorthand.instance_methods.map(&:to_sym)).to include(:whypg_pdf, :whypg_cdf, :whypg_r)
      expect(whypg_pdf(2, 5, 4, 12, 3.0)).to be_within(1e-12).of(@engine.pdf(2, 5, 4, 12, 3.0))
    end
  end
end
