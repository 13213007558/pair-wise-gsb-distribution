require 'spec_helper'

describe Distribution::Saddlepoint do
  # Continuous reference: Gamma(shape = 3, scale = 2).
  # K(t) = -a log(1 - b t), defined for t < 1/b.
  let(:gamma_shape) { 3.0 }
  let(:gamma_scale) { 2.0 }
  let(:gamma) do
    a = gamma_shape
    b = gamma_scale
    Distribution::Saddlepoint.new(
      cgf:        ->(t) { -a * Math.log(1 - b * t) },
      cgf_prime:  ->(t) { a * b / (1 - b * t) },
      cgf_second: ->(t) { a * b * b / (1 - b * t)**2 },
      cgf_third:  ->(t) { 2 * a * b**3 / (1 - b * t)**3 },
      domain: -Float::INFINITY...0.5
    )
  end

  # Lattice reference: Binomial(n = 20, p = 0.3), span 1.
  let(:binom_n) { 20 }
  let(:binom_p) { 0.3 }
  let(:binomial) do
    n = binom_n
    p = binom_p
    Distribution::Saddlepoint.new(
      cgf:        ->(t) { n * Math.log(1 - p + p * Math.exp(t)) },
      cgf_prime:  ->(t) { n * p * Math.exp(t) / (1 - p + p * Math.exp(t)) },
      cgf_second: ->(t) { n * p * Math.exp(t) * (1 - p) / (1 - p + p * Math.exp(t))**2 },
      lattice: 1
    )
  end

  describe 'continuous distribution (Gamma)' do
    it 'approximates the cdf in the right tail' do
      [10.0, 16.0, 25.0].each do |x|
        result = gamma.cdf(x)
        exact = Distribution::Gamma.cdf(x, gamma_shape, gamma_scale)
        expect(result.value).to be_within(1e-3).of(exact), "cdf(#{x}): #{result.value} vs #{exact}"
      end
    end

    it 'approximates the cdf in the left tail' do
      result = gamma.cdf(2.0)
      exact = Distribution::Gamma.cdf(2.0, gamma_shape, gamma_scale)
      expect(result.value).to be_within(1e-3).of(exact)
    end

    it 'approximates the survival function tail-stably' do
      result = gamma.sf(16.0)
      exact = 1 - Distribution::Gamma.cdf(16.0, gamma_shape, gamma_scale)
      expect(result.value).to be_within(1e-4).of(exact)
    end

    it 'approximates the pdf' do
      [6.0, 10.0].each do |x|
        result = gamma.pdf(x)
        exact = Distribution::Gamma.pdf(x, gamma_shape, gamma_scale)
        expect(result.value).to be_within(0.05 * exact).of(exact)
      end
    end

    it 'stays finite and accurate at the mean (removable 0/0 singularity)' do
      mean = gamma_shape * gamma_scale
      result = gamma.cdf(mean)
      expect(result.converged?).to be true
      expect(result.value).to be_finite
      expect(result.value).to be_within(1e-3).of(Distribution::Gamma.cdf(mean, gamma_shape, gamma_scale))
    end

    it 'stays accurate arbitrarily close to the mean' do
      mean = gamma_shape * gamma_scale
      [1e-3, 1e-7].each do |eps|
        result = gamma.cdf(mean + eps)
        exact = Distribution::Gamma.cdf(mean + eps, gamma_shape, gamma_scale)
        expect(result.value).to be_within(1e-3).of(exact), "cdf(mean+#{eps})"
      end
    end

    it 'works near the mean without a cgf_third callback' do
      a = gamma_shape
      b = gamma_scale
      spa = Distribution::Saddlepoint.new(
        cgf:        ->(t) { -a * Math.log(1 - b * t) },
        cgf_prime:  ->(t) { a * b / (1 - b * t) },
        cgf_second: ->(t) { a * b * b / (1 - b * t)**2 },
        domain: -Float::INFINITY...0.5
      )
      result = spa.cdf(a * b)
      expect(result.value).to be_within(1e-3).of(Distribution::Gamma.cdf(a * b, a, b))
    end

    it 'does not force a root outside the CGF domain' do
      result = gamma.cdf(-1.0)
      expect(result.converged?).to be false
      expect(result.value).to be_nil
      expect(result.messages).not_to be_empty
    end
  end

  describe 'lattice distribution (Binomial)' do
    it 'is declared as lattice' do
      expect(binomial.lattice?).to be true
      expect(gamma.lattice?).to be false
    end

    it 'approximates the cdf with continuity correction' do
      [3, 6, 8, 12, 15].each do |k|
        result = binomial.cdf(k)
        exact = Distribution::Binomial.exact_cdf(k, binom_n, binom_p)
        expect(result.value).to be_within(5e-3).of(exact), "cdf(#{k}): #{result.value} vs #{exact}"
      end
    end

    it 'approximates the cdf near the mean' do
      k = (binom_n * binom_p).floor
      result = binomial.cdf(k)
      exact = Distribution::Binomial.exact_cdf(k, binom_n, binom_p)
      expect(result.value).to be_within(5e-3).of(exact)
    end

    it 'approximates the pmf' do
      result = binomial.pmf(8)
      exact = Distribution::Binomial.pdf(8, binom_n, binom_p)
      expect(result.value).to be_within(0.05 * exact).of(exact)
    end

    it 'does not force a root beyond the support' do
      result = binomial.cdf(binom_n + 5)
      expect(result.converged?).to be false
      expect(result.value).to be_nil
      expect(result.messages).not_to be_empty
    end

    it 'raises pmf/pdf mismatches' do
      expect { gamma.pmf(1) }.to raise_error(RuntimeError)
      expect { binomial.pdf(1) }.to raise_error(RuntimeError)
    end
  end

  describe 'result object' do
    it 'is always flagged as an approximation and carries diagnostics' do
      result = gamma.cdf(10.0)
      expect(result.approximate?).to be true
      expect(result.converged?).to be true
      expect(result.saddlepoint).to be_within(1e-8).of(0.2)
      expect(result.iterations).to be > 0
      diagnostics = result.diagnostics
      expect(diagnostics[:approximate]).to be true
      expect(diagnostics[:converged]).to be true
      expect(diagnostics[:saddlepoint]).to eq(result.saddlepoint)
      expect(diagnostics[:messages]).to be_an(Array)
      expect(result.to_f).to eq(result.value)
    end
  end

  describe 'constructor validation' do
    let(:callbacks) do
      { cgf: ->(t) { t * t / 2 }, cgf_prime: ->(t) { t }, cgf_second: ->(_t) { 1.0 } }
    end

    it 'rejects a non-positive lattice span' do
      expect { Distribution::Saddlepoint.new(**callbacks, lattice: 0) }.to raise_error(ArgumentError)
    end

    it 'rejects a domain that does not contain 0' do
      expect { Distribution::Saddlepoint.new(**callbacks, domain: 1.0..2.0) }.to raise_error(ArgumentError)
    end

    it 'accepts a [min, max] pair as domain' do
      spa = Distribution::Saddlepoint.new(**callbacks, domain: [-2.0, 3.0])
      expect(spa.cdf(0.5).converged?).to be true
    end
  end
end
