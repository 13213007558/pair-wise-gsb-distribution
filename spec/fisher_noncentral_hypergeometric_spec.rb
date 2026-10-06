require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::FisherNoncentralHypergeometric do
  before do
    @engine = Distribution::FisherNoncentralHypergeometric::Ruby_
  end

  # Exact reference probabilities computed by enumerating every
  # n-subset of the population and weighting each subset by odds^(#type-1),
  # using exact Rational arithmetic.
  def exhaustive_probabilities(m1, m2, n, odds)
    w = Hash.new(0)
    (0...(m1 + m2)).to_a.combination(n) do |subset|
      k = subset.count { |i| i < m1 }
      w[k] += odds**k
    end
    total = w.values.inject(0) { |s, v| s + v }
    w.transform_values { |v| v / total }
  end

  shared_examples_for 'fisher noncentral hypergeometric engine' do
    it 'matches exhaustive weighted-subset enumeration on small populations' do
      [[3, 4, 3, Rational(5, 2)], [4, 2, 3, Rational(1, 3)], [5, 5, 4, Rational(7, 1)],
       [2, 6, 5, Rational(3, 4)]].each do |m1, m2, n, odds|
        exact = exhaustive_probabilities(m1, m2, n, odds)
        lo = [0, n - m2].max
        hi = [n, m1].min
        (lo..hi).each do |k|
          expect(@engine.pdf(k, m1, m2, n, odds.to_f)).to be_within(1e-12).of(exact[k].to_f)
        end
        # probabilities outside the support vanish
        expect(@engine.pdf(lo - 1, m1, m2, n, odds.to_f)).to eq(0.0)
        expect(@engine.pdf(hi + 1, m1, m2, n, odds.to_f)).to eq(0.0)
      end
    end

    it 'returns probabilities summing to one over the support' do
      [[10, 20, 7, 0.25], [30, 40, 25, 3.5], [100, 100, 50, 1.7]].each do |m1, m2, n, odds|
        lo = [0, n - m2].max
        hi = [n, m1].min
        total = (lo..hi).inject(0.0) { |s, k| s + @engine.pdf(k, m1, m2, n, odds) }
        expect(total).to be_within(1e-12).of(1.0)
      end
    end

    it 'reduces to the hypergeometric distribution when odds == 1' do
      [[80, 9920, 100], [10, 15, 7], [50, 50, 25]].each do |m1, m2, n|
        lo = [0, n - m2].max
        hi = [n, m1].min
        (lo..hi).each do |k|
          expected = Distribution::Hypergeometric.pdf(k, m1, n, m1 + m2).to_f
          expect(@engine.pdf(k, m1, m2, n, 1.0)).to be_within(1e-12).of(expected)
        end
      end
    end

    it 'satisfies the swap symmetry pmf(k; m1,m2,n,odds) == pmf(n-k; m2,m1,n,1/odds)' do
      [[5, 7, 4, 3.0], [8, 3, 6, 0.4], [20, 30, 15, 2.5]].each do |m1, m2, n, odds|
        lo = [0, n - m2].max
        hi = [n, m1].min
        (lo..hi).each do |k|
          swapped = @engine.pdf(n - k, m2, m1, n, 1.0 / odds)
          expect(@engine.pdf(k, m1, m2, n, odds)).to be_within(1e-12).of(swapped)
        end
        expect(@engine.mean(m1, m2, n, odds)).to be_within(1e-10).of(n - @engine.mean(m2, m1, n, 1.0 / odds))
        expect(@engine.variance(m1, m2, n, odds)).to be_within(1e-10).of(@engine.variance(m2, m1, n, 1.0 / odds))
      end
    end

    it 'handles the degenerate cases n == 0 and n == m1 + m2' do
      [0.1, 1.0, 10.0].each do |odds|
        # zero draws: point mass on k = 0
        expect(@engine.pdf(0, 5, 7, 0, odds)).to eq(1.0)
        expect(@engine.cdf(0, 5, 7, 0, odds)).to eq(1.0)
        expect(@engine.sf(0, 5, 7, 0, odds)).to eq(0.0)
        expect(@engine.quantile(0.5, 5, 7, 0, odds)).to eq(0)
        expect(@engine.mean(5, 7, 0, odds)).to eq(0.0)
        expect(@engine.variance(5, 7, 0, odds)).to eq(0.0)
        # drawing the whole population: point mass on k = m1
        expect(@engine.pdf(5, 5, 7, 12, odds)).to eq(1.0)
        expect(@engine.pdf(4, 5, 7, 12, odds)).to eq(0.0)
        expect(@engine.cdf(4, 5, 7, 12, odds)).to eq(0.0)
        expect(@engine.sf(5, 5, 7, 12, odds)).to eq(0.0)
        expect(@engine.quantile(0.99, 5, 7, 12, odds)).to eq(5)
        expect(@engine.mean(5, 7, 12, odds)).to eq(5.0)
        expect(@engine.variance(5, 7, 12, odds)).to eq(0.0)
      end
    end

    it 'preserves tails and moments under extreme odds' do
      m1 = 10
      m2 = 10
      n = 8
      # huge odds: essentially all mass on the largest possible k
      expect(@engine.pdf(8, m1, m2, n, 1e300)).to be_within(1e-12).of(1.0)
      expect(@engine.cdf(7, m1, m2, n, 1e300)).to be_within(1e-12).of(0.0)
      expect(@engine.mean(m1, m2, n, 1e300)).to be_within(1e-6).of(8.0)
      expect(@engine.variance(m1, m2, n, 1e300)).to be_within(1e-6).of(0.0)
      # tiny odds: essentially all mass on the smallest possible k
      expect(@engine.pdf(0, m1, m2, n, 1e-300)).to be_within(1e-12).of(1.0)
      expect(@engine.sf(0, m1, m2, n, 1e-300)).to be_within(1e-12).of(0.0)
      expect(@engine.mean(m1, m2, n, 1e-300)).to be_within(1e-6).of(0.0)
      expect(@engine.variance(m1, m2, n, 1e-300)).to be_within(1e-6).of(0.0)
      # large but moderate odds: non-degenerate tails must still sum to one
      [1e6, 1e-6, 1e100, 1e-100].each do |odds|
        total = (0..8).inject(0.0) { |s, k| s + @engine.pdf(k, m1, m2, n, odds) }
        expect(total).to be_within(1e-12).of(1.0)
        mean = @engine.mean(m1, m2, n, odds)
        expect(mean).to be >= 0.0
        expect(mean).to be <= 8.0
        expect(@engine.variance(m1, m2, n, odds)).to be >= 0.0
        expect(@engine.cdf(8, m1, m2, n, odds)).to eq(1.0)
      end
    end

    it 'computes cdf, sf and quantile consistently with the pdf' do
      m1 = 8
      m2 = 12
      n = 6
      odds = 2.0
      lo = [0, n - m2].max
      hi = [n, m1].min
      ac = 0.0
      (lo..hi).each do |k|
        ac += @engine.pdf(k, m1, m2, n, odds)
        expect(@engine.cdf(k, m1, m2, n, odds)).to be_within(1e-12).of(ac)
        expect(@engine.cdf(k, m1, m2, n, odds) + @engine.sf(k, m1, m2, n, odds)).to be_within(1e-12).of(1.0)
      end
      expect(@engine.cdf(lo - 1, m1, m2, n, odds)).to eq(0.0)
      expect(@engine.sf(hi, m1, m2, n, odds)).to eq(0.0)
      [0.01, 0.25, 0.5, 0.75, 0.99, 1.0].each do |pr|
        q = @engine.quantile(pr, m1, m2, n, odds)
        expect(@engine.cdf(q, m1, m2, n, odds)).to be >= pr
        expect(@engine.cdf(q - 1, m1, m2, n, odds)).to be < pr if q > lo
      end
    end

    it 'matches exact mean and variance from exhaustive enumeration' do
      m1 = 3
      m2 = 4
      n = 3
      odds = Rational(5, 2)
      exact = exhaustive_probabilities(m1, m2, n, odds)
      mu = exact.inject(0) { |s, (k, p)| s + k * p }
      var = exact.inject(0) { |s, (k, p)| s + (k - mu)**2 * p }
      expect(@engine.mean(m1, m2, n, odds.to_f)).to be_within(1e-12).of(mu.to_f)
      expect(@engine.variance(m1, m2, n, odds.to_f)).to be_within(1e-12).of(var.to_f)
    end

    it 'rejects invalid parameters' do
      expect { @engine.pdf(1, -1, 5, 2, 1.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 5, 5, 11, 1.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 5, 5, -1, 1.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 5, 5, 2, -0.5) }.to raise_error(ArgumentError)
    end
  end

  describe Distribution::FisherNoncentralHypergeometric::Ruby_ do
    it_should_behave_like 'fisher noncentral hypergeometric engine'
  end

  it 'exposes the methods on the Distribution module' do
    expect(Distribution::FisherNoncentralHypergeometric.pdf(2, 5, 5, 4, 2.0)).to be_within(1e-12).of(@engine.pdf(2, 5, 5, 4, 2.0))
    expect(Distribution::FisherNoncentralHypergeometric.cdf(2, 5, 5, 4, 2.0)).to be_within(1e-12).of(@engine.cdf(2, 5, 5, 4, 2.0))
    expect(Distribution::FisherNoncentralHypergeometric.sf(2, 5, 5, 4, 2.0)).to be_within(1e-12).of(@engine.sf(2, 5, 5, 4, 2.0))
    expect(Distribution::FisherNoncentralHypergeometric.quantile(0.5, 5, 5, 4, 2.0)).to eq(@engine.quantile(0.5, 5, 5, 4, 2.0))
    expect(Distribution::FisherNoncentralHypergeometric.mean(5, 5, 4, 2.0)).to be_within(1e-12).of(@engine.mean(5, 5, 4, 2.0))
    expect(Distribution::FisherNoncentralHypergeometric.variance(5, 5, 4, 2.0)).to be_within(1e-12).of(@engine.variance(5, 5, 4, 2.0))
  end
end
