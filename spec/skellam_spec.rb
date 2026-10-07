require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::Skellam do
  def log_sum_exp(values)
    values = values.reject { |value| value == -Float::INFINITY }
    return -Float::INFINITY if values.empty?
    maximum = values.max
    maximum + Math.log(values.inject(0.0) { |sum, value| sum + Math.exp(value - maximum) })
  end

  def poisson_log_pmf(k, mean)
    return -Float::INFINITY if k < 0
    return k.zero? ? 0.0 : -Float::INFINITY if mean.zero?
    -mean + k * Math.log(mean) - Math.lgamma(k + 1).first
  end

  def poisson_bound(mean)
    (mean + 40 * Math.sqrt(mean) + 80).ceil
  end

  def convolution_terms(mean1, mean2)
    upper1 = poisson_bound(mean1)
    upper2 = poisson_bound(mean2)
    terms = Hash.new { |hash, key| hash[key] = [] }
 
    0.upto(upper1) do |x|
      0.upto(upper2) do |y|
        terms[x - y] << poisson_log_pmf(x, mean1) + poisson_log_pmf(y, mean2)
      end
    end
 
    terms
  end
 
  def reference_probabilities(mean1, mean2)
    terms = convolution_terms(mean1, mean2)
    probabilities = Hash.new(0.0)
    terms.each do |difference, log_terms|
      probabilities[difference] = Math.exp(log_sum_exp(log_terms))
    end
    probabilities
  end

  def reference_log_pmf(k, mean1, mean2)
    upper2 = poisson_bound(mean2)
    upper1 = [poisson_bound(mean1), k + upper2].max
    lower_x = [0, k].max
    terms = []

    lower_x.upto(upper1) do |x|
      y = x - k
      terms << poisson_log_pmf(x, mean1) + poisson_log_pmf(y, mean2)
    end

    log_sum_exp(terms)
  end

  def reference_cdf(k, mean1, mean2)
    reference_probabilities(mean1, mean2).select { |difference, _| difference <= k }.values.sum
  end

  def reference_survival(k, mean1, mean2)
    reference_probabilities(mean1, mean2).select { |difference, _| difference > k }.values.sum
  end

  shared_examples_for 'skellam engine' do
    [[0.5, 1.5], [1.0, 2.0], [2.0, 0.7], [0.3, 5.0], [4.0, 4.0]].each do |mean1, mean2|
      (-10..10).each do |k|
        it "calculates log_pmf for #{k} with means #{mean1}, #{mean2}" do
          expect(@engine.log_pmf(k, mean1, mean2)).to(
            be_within(1e-12).of(reference_log_pmf(k, mean1, mean2))
          )
        end

        it "calculates cdf for #{k} with means #{mean1}, #{mean2}" do
          expect(@engine.cdf(k, mean1, mean2)).to(
            be_within(1e-11).of(reference_cdf(k, mean1, mean2))
          )
        end

        it "calculates survival for #{k} with means #{mean1}, #{mean2}" do
          expect(@engine.survival(k, mean1, mean2)).to(
            be_within(1e-11).of(reference_survival(k, mean1, mean2))
          )
        end
      end
    end

    it 'supports very asymmetric means and tiny right tails' do
      probability = @engine.survival(10, 0.2, 20.0)
      expect(probability).to be_within(1e-14).of(reference_survival(10, 0.2, 20.0))
      expect(probability).to be < 1e-8
    end

    it 'mirrors pmf when the means are exchanged' do
      (-12..12).each do |k|
        expect(@engine.log_pmf(k, 0.4, 6.0)).to(
          be_within(1e-13).of(@engine.log_pmf(-k, 6.0, 0.4))
        )
      end
    end

    it 'mirrors cdf and survival when the means are exchanged' do
      (-12..12).each do |k|
        expect(@engine.cdf(k, 0.4, 6.0)).to(
          be_within(1e-12).of(@engine.survival(-k - 1, 6.0, 0.4))
        )
        expect(@engine.survival(k, 0.4, 6.0)).to(
          be_within(1e-12).of(@engine.cdf(-k - 1, 6.0, 0.4))
        )
      end
    end

    it 'handles both means equal to zero as a point mass' do
      expect(@engine.log_pmf(0, 0, 0)).to eq(0.0)
      expect(@engine.log_pmf(1, 0, 0)).to eq(-Float::INFINITY)
      expect(@engine.cdf(-1, 0, 0)).to eq(0.0)
      expect(@engine.cdf(0, 0, 0)).to eq(1.0)
      expect(@engine.survival(-1, 0, 0)).to eq(1.0)
      expect(@engine.survival(0, 0, 0)).to eq(0.0)
      expect(@engine.quantile(0.5, 0, 0)).to eq(0)
    end

    it 'degenerates to a shifted Poisson when one mean is zero' do
      [0, 1, 2, 5].each do |k|
        expected = poisson_log_pmf(k, 2.5)
        expect(@engine.log_pmf(k, 2.5, 0)).to be_within(1e-14).of(expected)
      end
      [-1, 0, 1, 2].each do |k|
        expect(@engine.cdf(k, 2.5, 0)).to be_within(1e-12).of(reference_cdf(k, 2.5, 0))
        expect(@engine.survival(k, 2.5, 0)).to be_within(1e-12).of(reference_survival(k, 2.5, 0))
      end
      [0, 1, 2].each do |k|
        expect(@engine.log_pmf(-k, 0, 2.5)).to be_within(1e-14).of(poisson_log_pmf(k, 2.5))
      end
    end

    it 'returns the discrete lower-tail quantile' do
      [0.05, 0.25, 0.5, 0.75, 0.95].each do |probability|
        k = @engine.quantile(probability, 1.5, 2.5)
        expect(@engine.cdf(k, 1.5, 2.5)).to be >= probability
        expect(@engine.cdf(k - 1, 1.5, 2.5)).to be < probability
      end
    end

    it 'rejects invalid means and probabilities' do
      expect { @engine.log_pmf(0, -1, 1) }.to raise_error(ArgumentError)
      expect { @engine.cdf(0, 1, -1) }.to raise_error(ArgumentError)
      expect { @engine.quantile(1.1, 1, 1) }.to raise_error(ArgumentError)
      expect { @engine.quantile(-0.1, 1, 1) }.to raise_error(ArgumentError)
    end

    it 'returns -Infinity for non-integer pmf values' do
      expect(@engine.log_pmf(1.5, 1.0, 1.0)).to eq(-Float::INFINITY)
    end
  end

  describe 'singleton' do
    before do
      @engine = Distribution::Skellam
    end
    it_should_behave_like 'skellam engine'
  end

  describe Distribution::Skellam::Ruby_ do
    before do
      @engine = Distribution::Skellam::Ruby_
    end
    it_should_behave_like 'skellam engine'
  end
end
