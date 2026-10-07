require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::Skellam do
  def log_poisson_pmf(k, mean)
    return -Float::INFINITY if k < 0
    return k.zero? ? 0.0 : -Float::INFINITY if mean.zero?
    k * Math.log(mean) - mean - Math.lgamma(k + 1).first
  end

  def log_add(first, second)
    return second if first == -Float::INFINITY
    return first if second == -Float::INFINITY
    maximum = [first, second].max
    maximum + Math.log(Math.exp(first - maximum) + Math.exp(second - maximum))
  end

  def reference_pmf(k, mu1, mu2)
    return 0.0 unless k.is_a?(Integer)
    lower = [0, -k].max
    log_probability = -Float::INFINITY
    lower.upto(40) do |n|
      log_probability = log_add(
        log_probability,
        log_poisson_pmf(n, mu2) + log_poisson_pmf(n + k, mu1)
      )
    end
    Math.exp(log_probability)
  end

  def reference_cdf(k, mu1, mu2)
    (-30..k).sum { |n| reference_pmf(n, mu1, mu2) }
  end

  def reference_survival(k, mu1, mu2)
    (k + 1..30).sum { |n| reference_pmf(n, mu1, mu2) }
  end

  shared_examples_for 'skellam engine' do
    it 'matches an explicit convolution of small independent Poisson variables' do
      [[0.5, 0.25], [1.0, 2.0], [2.5, 0.75], [4.0, 1.0]].each do |mu1, mu2|
        (-12..12).each do |k|
          probability = reference_pmf(k, mu1, mu2)
          expect(@engine.pmf(k, mu1, mu2)).to be_within(1e-12).of(probability)
          expect(Math.exp(@engine.log_pmf(k, mu1, mu2))).to be_within(1e-12).of(probability)
          expect(@engine.cdf(k, mu1, mu2)).to be_within(1e-10).of(reference_cdf(k, mu1, mu2))
          expect(@engine.survival(k, mu1, mu2)).to be_within(1e-10).of(reference_survival(k, mu1, mu2))
        end
      end
    end

    it 'mirrors when the Poisson means are swapped' do
      [[1.5, 0.5], [10.0, 0.01], [0.2, 4.0]].each do |mu1, mu2|
        (-10..10).each do |k|
          expect(@engine.pmf(k, mu1, mu2)).to be_within(1e-14).of(@engine.pmf(-k, mu2, mu1))
          expect(@engine.cdf(k, mu1, mu2)).to be_within(1e-12).of(@engine.survival(-k - 1, mu2, mu1))
          expect(@engine.survival(k, mu1, mu2)).to be_within(1e-12).of(@engine.cdf(-k - 1, mu2, mu1))
        end
      end
    end

    it 'supports strongly asymmetric means with direct tails' do
      expect(@engine.cdf(0, 10_000.0, 0.01)).to eq(0.0)
      expect(@engine.survival(0, 10_000.0, 0.01)).to eq(1.0)
      expect(@engine.cdf(0, 0.01, 10_000.0)).to eq(1.0)
      expect(@engine.survival(0, 0.01, 10_000.0)).to eq(0.0)
      expect(@engine.survival(-1, 10_000.0, 0.01)).to eq(1.0)
      expect(@engine.cdf(-1, 0.01, 10_000.0)).to eq(1.0)
      expect(@engine.survival(-1, 0.01, 10_000.0)).to eq(0.0)
    end

    it 'reduces to a one-sided Poisson when one mean is zero' do
      [0.5, 2.0, 4.0].each do |mean|
        expect(@engine.pmf(2, mean, 0.0)).to be_within(1e-14).of(Distribution::Poisson.pdf(2, mean))
        expect(@engine.pmf(-2, 0.0, mean)).to be_within(1e-14).of(Distribution::Poisson.pdf(2, mean))
        expect(@engine.pmf(-1, mean, 0.0)).to eq(0.0)
        expect(@engine.pmf(1, 0.0, mean)).to eq(0.0)
        expect(@engine.cdf(-1, 0.0, mean)).to be_within(1e-14).of(1.0 - Distribution::Poisson.pdf(0, mean))
        expect(@engine.survival(0, mean, 0.0)).to be_within(1e-14).of(1.0 - Distribution::Poisson.cdf(0, mean))
      end
    end

    it 'has point mass at zero when both means are zero' do
      expect(@engine.pmf(0, 0.0, 0.0)).to eq(1.0)
      expect(@engine.log_pmf(0, 0.0, 0.0)).to eq(0.0)
      expect(@engine.pmf(-1, 0.0, 0.0)).to eq(0.0)
      expect(@engine.log_pmf(-1, 0.0, 0.0)).to eq(-Float::INFINITY)
      expect(@engine.cdf(-1, 0.0, 0.0)).to eq(0.0)
      expect(@engine.cdf(0, 0.0, 0.0)).to eq(1.0)
      expect(@engine.survival(-1, 0.0, 0.0)).to eq(1.0)
      expect(@engine.survival(0, 0.0, 0.0)).to eq(0.0)
    end

    it 'treats non-integer values by floor for CDF and gives zero mass between integers' do
      expect(@engine.pmf(1.2, 1.0, 1.0)).to eq(0.0)
      expect(@engine.log_pmf(1.2, 1.0, 1.0)).to eq(-Float::INFINITY)
      expect(@engine.cdf(1.7, 1.0, 1.0)).to eq(@engine.cdf(1, 1.0, 1.0))
      expect(@engine.survival(1.7, 1.0, 1.0)).to eq(@engine.survival(1, 1.0, 1.0))
    end

    it 'returns the smallest integer with CDF at or above the probability' do
      [[1.5, 0.5], [0.5, 3.0], [4.0, 0.0], [0.0, 2.0]].each do |mu1, mu2|
        [0.05, 0.5, 0.95].each do |probability|
          q = @engine.quantile(probability, mu1, mu2)
          expect(@engine.cdf(q - 1, mu1, mu2)).to be < probability
          expect(@engine.cdf(q, mu1, mu2)).to be >= probability
        end
      end
    end

    it 'validates means and probabilities' do
      expect { @engine.log_pmf(0, -1, 1) }.to raise_error(ArgumentError)
      expect { @engine.cdf(0, 1, -1) }.to raise_error(ArgumentError)
      expect { @engine.quantile(-0.1, 1, 1) }.to raise_error(ArgumentError)
      expect { @engine.quantile(1.1, 1, 1) }.to raise_error(ArgumentError)
    end

    it 'supports p_value and pdf aliases' do
      expect(@engine.pdf(0, 1.0, 2.0)).to be_within(1e-14).of(@engine.pmf(0, 1.0, 2.0))
      expect(@engine.p_value(0.5, 1.0, 2.0)).to eq(@engine.quantile(0.5, 1.0, 2.0))
    end

    it 'provides shorthand methods' do
      shorthand = Class.new { include Distribution::Shorthand }.new
      expect(shorthand.skel_pdf(0, 1.0, 2.0)).to be_within(1e-14).of(@engine.pmf(0, 1.0, 2.0))
      expect(shorthand.skel_lpmf(0, 1.0, 2.0)).to be_within(1e-14).of(@engine.log_pmf(0, 1.0, 2.0))
      expect(shorthand.skel_sf(0, 1.0, 2.0)).to be_within(1e-14).of(@engine.survival(0, 1.0, 2.0))
    end
  end

  describe 'singleton' do
    before { @engine = Distribution::Skellam }
    it_should_behave_like 'skellam engine'
  end

  describe Distribution::Skellam::Ruby_ do
    before { @engine = Distribution::Skellam::Ruby_ }
    it_should_behave_like 'skellam engine'
  end
end
