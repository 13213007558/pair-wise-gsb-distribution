require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::KolmogorovSmirnov do
  shared_examples_for 'finite-sample Kolmogorov-Smirnov engine' do
    it 'provides the finite-n support boundaries' do
      expect(@engine.cdf(0.0, 1)).to eq(0.0)
      expect(@engine.cdf(0.5, 1)).to eq(0.0)
      expect(@engine.cdf(1.0, 1)).to eq(1.0)
      expect(@engine.sf(0.0, 1)).to eq(1.0)
      expect(@engine.sf(0.5, 1)).to eq(1.0)
      expect(@engine.sf(1.0, 1)).to eq(0.0)
    end

    it 'matches the directly integrated n=1 CDF' do
      expect(@engine.cdf(0.75, 1)).to be_within(1e-15).of(0.5)
      expect(@engine.sf(0.75, 1)).to be_within(1e-15).of(0.5)
      [0.6, 0.8, 0.99].each do |d|
        expect(@engine.cdf(d, 1)).to be_within(1e-14).of(2 * d - 1)
        expect(@engine.sf(d, 1)).to be_within(1e-14).of(2 * (1 - d))
      end
    end

    it 'matches the directly integrated n=2 piecewise CDF' do
      expect(@engine.cdf(0.25, 2)).to eq(0.0)
      expect(@engine.cdf(0.4, 2)).to be_within(1e-14).of(8 * (0.4 - 0.25)**2)
      expect(@engine.cdf(0.5, 2)).to be_within(1e-14).of(0.5)
      expect(@engine.cdf(0.75, 2)).to be_within(1e-14).of(0.875)
      expect(@engine.sf(0.75, 2)).to be_within(1e-14).of(0.125)
      expect(@engine.cdf(0.99, 2)).to be_within(1e-14).of(1 - 2 * 0.01**2)
    end

    it 'matches independent exact values for medium n' do
      values = {
        [5, 0.2] => [0.0384, 0.9616],
        [5, 0.5] => [0.888, 0.112],
        [5, 0.6] => [0.96992, 0.03008],
        [10, 0.15] => [0.04603473, 0.95396527],
        [10, 0.5] => [0.99222259, 0.00777741],
        [20, 0.06] => [6.152693947325927e-6, 0.9999938473060527],
        [20, 0.7] => [0.9999999997607605, 2.392395040880019e-10],
        [20, 0.9] => [1.0, 2.0006866455078037e-20]
      }

      values.each do |(n, d), (cdf, sf)|
        expect(@engine.cdf(d, n)).to be_within([1e-15, cdf.abs * 1e-12].max).of(cdf)
        expect(@engine.sf(d, n)).to be_within([1e-16, sf.abs * 1e-12].max).of(sf)
      end
    end

    it 'keeps the small right tail accurate near the upper support' do
      expect(@engine.sf(0.99, 5)).to be_within(1e-24).of(2e-10)
      expect(@engine.sf(0.99, 20)).to be_within(1e-34).of(2e-40)
      expect(@engine.sf(0.99, 100)).to be_within(1e-210).of(2e-200)
    end

    it 'is monotone in d' do
      [1, 2, 5, 10, 20, 50].each do |n|
        values = 0.upto(30).map { |i| 0.5 / n + i * (1.0 - 0.5 / n) / 30 }
        cdfs = values.map { |d| @engine.cdf(d, n) }
        sfs = values.map { |d| @engine.sf(d, n) }
        expect(cdfs.each_cons(2).all? { |left, right| right >= left - 1e-15 }).to be(true)
        expect(sfs.each_cons(2).all? { |left, right| right <= left + 1e-15 }).to be(true)
      end
    end

    it 'returns quantiles that bracket the requested probability' do
      [0.0, 0.0384, 0.5, 0.9, 0.999, 1.0].each do |p|
        q = @engine.quantile(p, 5)
        expect(@engine.cdf(q, 5)).to be >= p - 1e-14
        expect(@engine.sf(q, 5)).to be >= (1 - p) - 1e-14
      end
    end

    it 'rejects invalid sample sizes and probabilities' do
      expect { @engine.cdf(0.5, 0) }.to raise_error(ArgumentError)
      expect { @engine.cdf(0.5, 1.5) }.to raise_error(ArgumentError)
      expect { @engine.quantile(-0.1, 5) }.to raise_error(RangeError)
      expect { @engine.quantile(1.1, 5) }.to raise_error(RangeError)
    end
  end

  describe 'singleton' do
    before { @engine = Distribution::KolmogorovSmirnov }
    it_should_behave_like 'finite-sample Kolmogorov-Smirnov engine'
  end

  describe Distribution::KolmogorovSmirnov::Ruby_ do
    before { @engine = Distribution::KolmogorovSmirnov::Ruby_ }
    it_should_behave_like 'finite-sample Kolmogorov-Smirnov engine'
  end
end
