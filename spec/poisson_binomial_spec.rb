require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::PoissonBinomial do
  # Exact pmf obtained by enumerating all 2^n success/failure combinations
  def brute_force_pmf(probs)
    pmf = Hash.new(0.0)
    (0...(1 << probs.size)).each do |mask|
      prob = 1.0
      successes = 0
      probs.each_with_index do |p, i|
        if mask[i] == 1
          prob *= p
          successes += 1
        else
          prob *= (1 - p)
        end
      end
      pmf[successes] += prob
    end
    pmf
  end

  shared_examples_for 'poisson binomial engine' do
    it 'should match brute-force enumeration for short lists' do
      [
        [0.1, 0.2, 0.3, 0.4],
        [0.5, 0.5, 0.5],
        [0.01, 0.99, 0.5, 0.25, 0.75],
        [0.0, 1.0, 0.3, 0.7]
      ].each do |probs|
        expected = brute_force_pmf(probs)
        n = probs.size
        (0..n).each do |k|
          obs = @engine.pdf(k, probs)
          exp = expected[k]
          expect(obs).to be_within(1e-12).of(exp), "For pdf(#{k}, #{probs}) expected #{exp}, obtained #{obs}"
        end
        (0..n).each do |k|
          exp_cdf = (0..k).inject(0.0) { |ac, i| ac + expected[i] }
          obs_cdf = @engine.cdf(k, probs)
          expect(obs_cdf).to be_within(1e-12).of(exp_cdf), "For cdf(#{k}, #{probs}) expected #{exp_cdf}, obtained #{obs_cdf}"

          exp_sf = (k..n).inject(0.0) { |ac, i| ac + expected[i] }
          obs_sf = @engine.sf(k, probs)
          expect(obs_sf).to be_within(1e-12).of(exp_sf), "For sf(#{k}, #{probs}) expected #{exp_sf}, obtained #{obs_sf}"
        end
      end
    end

    it 'should reduce to the Binomial distribution when all probabilities are equal' do
      [1, 5, 20, 100].each do |n|
        [0.25, 0.5, 0.8].each do |pr|
          probs = Array.new(n, pr)
          [0, 1, n / 2, n - 1, n].each do |k|
            exp_pdf = Distribution::Binomial.pdf(k, n, pr)
            obs_pdf = @engine.pdf(k, probs)
            expect(obs_pdf).to be_within(1e-10).of(exp_pdf), "For pdf(#{k}, #{n}x#{pr}) expected #{exp_pdf}, obtained #{obs_pdf}"
          end
          # Binomial.cdf is not defined for k = n (regularized_beta domain)
          [0, 1, n / 2, n - 1].uniq.select { |k| k < n }.each do |k|
            exp_cdf = Distribution::Binomial.cdf(k, n, pr)
            obs_cdf = @engine.cdf(k, probs)
            expect(obs_cdf).to be_within(1e-10).of(exp_cdf), "For cdf(#{k}, #{n}x#{pr}) expected #{exp_cdf}, obtained #{obs_cdf}"
          end
        end
      end
    end

    it 'should match Binomial quantiles when all probabilities are equal' do
      n = 20
      pr = 0.3
      probs = Array.new(n, pr)
      [0.05, 0.25, 0.5, 0.9, 0.99].each do |q|
        expect(@engine.p_value(q, probs)).to eq(Distribution::Binomial.p_value(q, n, pr))
      end
    end

    it 'should treat 0 and 1 probabilities as deterministic events' do
      probs = [0.0, 1.0, 1.0, 0.3]
      expect(@engine.pdf(2, probs)).to be_within(1e-15).of(0.7)
      expect(@engine.pdf(3, probs)).to be_within(1e-15).of(0.3)
      expect(@engine.pdf(0, probs)).to eq(0.0)
      expect(@engine.pdf(1, probs)).to eq(0.0)
      expect(@engine.cdf(1, probs)).to eq(0.0)
      expect(@engine.cdf(2, probs)).to be_within(1e-15).of(0.7)
      expect(@engine.cdf(3, probs)).to eq(1.0)
      expect(@engine.sf(2, probs)).to eq(1.0)
      expect(@engine.sf(3, probs)).to be_within(1e-15).of(0.3)
      expect(@engine.sf(4, probs)).to eq(0.0)
    end

    it 'should handle fully deterministic inputs' do
      expect(@engine.pdf(3, [1.0, 1.0, 1.0])).to eq(1.0)
      expect(@engine.cdf(2, [1.0, 1.0, 1.0])).to eq(0.0)
      expect(@engine.cdf(3, [1.0, 1.0, 1.0])).to eq(1.0)
      expect(@engine.sf(3, [1.0, 1.0, 1.0])).to eq(1.0)
      expect(@engine.sf(4, [1.0, 1.0, 1.0])).to eq(0.0)
      expect(@engine.pdf(0, [0.0, 0.0])).to eq(1.0)
      expect(@engine.sf(1, [0.0, 0.0])).to eq(0.0)
      expect(@engine.p_value(0.5, [1.0, 1.0, 1.0])).to eq(3)
    end

    it 'should handle an empty probability list as a point mass at 0' do
      expect(@engine.pdf(0, [])).to eq(1.0)
      expect(@engine.pdf(1, [])).to eq(0.0)
      expect(@engine.cdf(0, [])).to eq(1.0)
      expect(@engine.sf(1, [])).to eq(0.0)
    end

    it 'should return a normalized pmf whose total mass is 1' do
      probs = Array.new(200) { |i| (i + 1) / 201.0 }
      pmf = @engine.pmf_array(probs)
      expect(pmf.size).to eq(201)
      expect(pmf.inject(0.0) { |ac, m| ac + m }).to be_within(1e-12).of(1.0)
      expect(@engine.cdf(200, probs)).to eq(1.0)
    end

    it 'should keep tiny right tails accurate where 1 - cdf underflows' do
      # P(X >= 60) = 2^-60 ~ 8.7e-19; 1 - cdf(59) rounds to exactly 0.0
      probs = Array.new(60, 0.5)
      expected = 0.5**60
      expect(1.0 - @engine.cdf(59, probs)).to eq(0.0)
      expect(@engine.sf(60, probs)).to be_within(expected * 1e-6).of(expected)
      expect(@engine.sf(60, probs)).to be > 0.0
    end

    it 'should have cdf and sf complement each other' do
      probs = [0.05, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9]
      (0..probs.size).each do |k|
        expect(@engine.cdf(k - 1, probs) + @engine.sf(k, probs)).to be_within(1e-12).of(1.0)
      end
    end

    it 'should return quantiles consistent with the cdf' do
      probs = [0.1, 0.3, 0.5, 0.7, 0.9]
      (0..probs.size).each do |k|
        cdf = @engine.cdf(k, probs)
        expect(@engine.p_value(cdf, probs)).to be <= k
        expect(@engine.cdf(@engine.p_value(cdf, probs), probs)).to be >= cdf - 1e-12
      end
      expect(@engine.p_value(0.0, probs)).to eq(0)
      expect(@engine.p_value(1.0, probs)).to eq(probs.size)
    end

    it 'should return 0 or 1 outside of the support' do
      probs = [0.2, 0.4]
      expect(@engine.pdf(-1, probs)).to eq(0.0)
      expect(@engine.pdf(3, probs)).to eq(0.0)
      expect(@engine.pdf(1.5, probs)).to eq(0.0)
      expect(@engine.cdf(-1, probs)).to eq(0.0)
      expect(@engine.cdf(2, probs)).to eq(1.0)
      expect(@engine.cdf(10, probs)).to eq(1.0)
      expect(@engine.sf(0, probs)).to eq(1.0)
      expect(@engine.sf(-3, probs)).to eq(1.0)
      expect(@engine.sf(3, probs)).to eq(0.0)
    end

    it 'should raise ArgumentError on invalid probabilities' do
      expect { @engine.pdf(1, [0.5, -0.1]) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, [0.5, 1.1]) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, [0.5, Float::NAN]) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, [0.5, 'a']) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 0.5) }.to raise_error(ArgumentError)
      expect { @engine.cdf(1, [2]) }.to raise_error(ArgumentError)
      expect { @engine.sf(1, [-1]) }.to raise_error(ArgumentError)
      expect { @engine.p_value(1.5, [0.5]) }.to raise_error(ArgumentError)
      expect { @engine.p_value(-0.1, [0.5]) }.to raise_error(ArgumentError)
    end

    it 'should accept integer and rational probabilities' do
      expect(@engine.pdf(1, [1, 0])).to eq(1.0)
      expect(@engine.pdf(1, [1.quo(2), 1.quo(2)])).to be_within(1e-15).of(0.5)
    end
  end

  describe 'singleton' do
    before do
      @engine = Distribution::PoissonBinomial
    end

    it_should_behave_like 'poisson binomial engine'

    it { expect(@engine).to respond_to(:pdf) }
    it { expect(@engine).to respond_to(:cdf) }
    it { expect(@engine).to respond_to(:sf) }
    it { expect(@engine).to respond_to(:p_value) }
    it { expect(@engine).to respond_to(:quantile) }
    it { expect(@engine).to respond_to(:pmf_array) }
  end

  describe Distribution::PoissonBinomial::Ruby_ do
    before do
      @engine = Distribution::PoissonBinomial::Ruby_
    end

    it_should_behave_like 'poisson binomial engine'
  end
end
