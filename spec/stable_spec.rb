require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::Stable do
  let(:engine) { Distribution::Stable::Ruby_ }

  describe 'parameter conversions' do
    it 'converts S0 to S1 away from alpha one' do
      expect(engine.s1_parameters(1.5, 0.4, 2.0, 1.0))
        .to eq([1.5, 0.4, 2.0, 1.8000000000000003])
    end

    it 'converts S1 back to S0 away from alpha one' do
      s1 = engine.s1_parameters(1.5, 0.4, 2.0, 1.0)
      expect(engine.from_s1_parameters(*s1)).to eq([1.5, 0.4, 2.0, 1.0])
    end

    it 'uses the explicit alpha one logarithmic adjustment' do
      scale = Math.exp(1.0)
      s1 = engine.s1_parameters(1.0, 0.4, scale, 0.0)
      expect(s1.last).to be_within(1e-15).of(-0.6922047835458122)
    end

    it 'keeps S0 and S1 aligned under alpha one scaling' do
      scale = Math.exp(1.0)
      s0 = [1.0, 0.4, scale, 0.3]
      s1 = engine.s1_parameters(*s0)
      expect(engine.from_s1_parameters(*s1)).to eq(s0)
    end

    it 'validates the four parameters' do
      expect { engine.pdf(0, 0, 0) }.to raise_error(ArgumentError)
      expect { engine.pdf(0, 2.1, 0) }.to raise_error(ArgumentError)
      expect { engine.pdf(0, 1, 1.1) }.to raise_error(ArgumentError)
      expect { engine.pdf(0, 1, 0, 0) }.to raise_error(ArgumentError)
    end
  end

  describe 'density and distribution functions' do
    it 'matches the normal special case' do
      x = 0.7
      expect(engine.pdf(x, 2, 0)).to be_within(1e-12).of(
        Distribution::Normal.pdf(x / Math.sqrt(2)) / Math.sqrt(2)
      )
      expect(engine.cdf(x, 2, 0)).to be_within(1e-12).of(
        Distribution::Normal.cdf(x / Math.sqrt(2))
      )
    end

    it 'matches the Cauchy special case' do
      x = 0.3
      expect(engine.pdf(x, 1, 0)).to be_within(1e-12).of(1.0 / (Math::PI * (1 + x**2)))
      expect(engine.cdf(x, 1, 0)).to be_within(1e-12).of(0.5 + Math.atan(x) / Math::PI)
    end

    it 'applies location and scale' do
      alpha = 1.4
      beta = -0.3
      scale = 0.8
      location = 0.2
      z = 0.5
      x = location + scale * z
      expect(engine.pdf(x, alpha, beta, scale, location))
        .to be_within(1e-10).of(engine.pdf(z, alpha, beta) / scale)
      expect(engine.cdf(x, alpha, beta, scale, location))
        .to be_within(1e-9).of(engine.cdf(z, alpha, beta))
    end

    it 'has continuous location under S0 scaling near and at alpha one' do
      scale = 1.5
      [0.999, 1.0, 1.001].each do |alpha|
        z = 0.0
        x = 2.0 * scale * z
        expect(engine.cdf(x, alpha, 0.4, scale, 0.0, epsilon: 1e-10))
          .to be_within(0.01).of(engine.cdf(z, 1.0, 0.4, 1.0, 0.0, epsilon: 1e-10))
      end
    end

    [0.999, 1.0, 1.001].each do |alpha|
      it "is continuous at and around alpha=1 for beta=0.4 (#{alpha})" do
        near = engine.cdf(0.5, alpha, 0.4, 1.0, 0.0, epsilon: 1e-10)
        exact = engine.cdf(0.5, 1.0, 0.4, 1.0, 0.0, epsilon: 1e-10)
        expect(near).to be_within(2e-3).of(exact)
      end
    end

    [[0.99, -0.7], [1.01, -0.7], [0.99, 0.8], [1.01, 0.8]].each do |alpha, beta|
      it "keeps location continuous for alpha=#{alpha} and beta=#{beta}" do
        median = engine.cdf(0.0, alpha, beta, 1.0, 0.0, epsilon: 1e-10)
        reference = engine.cdf(0.0, 1.0, beta, 1.0, 0.0, epsilon: 1e-10)
        expect((median - reference).abs).to be < 0.02
      end
    end

    it 'returns diagnostics without raising on success' do
      result = engine.cdf(0.2, 1.5, 0.4, 1, 0, diagnostics: true, epsilon: 1e-8)
      expect(result[:converged]).to be(true)
      expect(result[:error]).to be < 1e-8
      expect(result[:evaluations]).to be > 0
    end

    it 'diagnoses non-convergence when the interval budget is exhausted' do
      options = { diagnostics: true, max_intervals: 0, epsilon: 1e-14, cutoff: 12 }
      result = engine.pdf(1.2, 1.3, 0.7, 1, 0, options)
      expect(result[:converged]).to be(false)
      expect(result[:subdivisions]).to eq(0)
      expect(result[:error]).to be >= 0
    end

    it 'raises with diagnostics attached when precision fails' do
      expect do
        engine.pdf(1.2, 1.3, 0.7, 1, 0, max_intervals: 0, epsilon: 1e-14, cutoff: 12)
      end.to raise_error(Distribution::Stable::IntegrationError) do |error|
        expect(error.diagnostics[:converged]).to be(false)
      end
    end
  end

  describe 'sampling' do
    def quantile(samples, probability)
      ordered = samples.sort
      index = (probability * (ordered.length - 1)).round
      ordered[index]
    end

    it 'accepts an injected random source' do
      source = Object.new
      values = [0.25, -Math.log(0.5)]
      index = 0
      source.define_singleton_method(:rand) { values[index].tap { index += 1 } }
      generator = engine.rng(1, 0, 1, 0, source)
      expect(generator.call).to be_within(1e-12).of(Math.tan(Math::PI * (0.25 - 0.5)))
    end

    it 'accepts a callable random source' do
      values = [0.75, 0.25].cycle
      generator = engine.rng(1, 0, 1, 0, values)
      expect(generator.call).to be_within(1e-12).of(Math.tan(Math::PI * 0.25))
    end

    it 'produces fixed-seed normal quantiles' do
      generator = engine.rng(2, 0, Math.sqrt(2), 0, Random.new(1234))
      samples = Array.new(50_000) { generator.call }
      expect(quantile(samples, 0.05)).to be_within(0.08).of(-2.326)
      expect(quantile(samples, 0.5)).to be_within(0.03).of(0.0)
      expect(quantile(samples, 0.95)).to be_within(0.08).of(2.326)
    end

    it 'produces fixed-seed Cauchy quantiles' do
      generator = engine.rng(1, 0, 1, 0, Random.new(1234))
      samples = Array.new(50_000) { generator.call }
      expect(quantile(samples, 0.05)).to be_within(0.2).of(-6.314)
      expect(quantile(samples, 0.5)).to be_within(0.03).of(0.0)
      expect(quantile(samples, 0.95)).to be_within(0.2).of(6.314)
    end

    it 'produces fixed-seed skewed quantiles near and at alpha one' do
      generator = engine.rng(1, 0.4, 1, 0, Random.new(1234))
      samples = Array.new(50_000) { generator.call }
      expect(quantile(samples, 0.05)).to be_within(0.18).of(-3.563)
      expect(quantile(samples, 0.5)).to be_within(0.04).of(0.169)
      expect(quantile(samples, 0.95)).to be_within(0.25).of(9.296)
    end

    it 'produces fixed-seed non-Gaussian skewed quantiles' do
      generator = engine.rng(1.5, -0.5, 1, 0, Random.new(1234))
      samples = Array.new(50_000) { generator.call }
      expect(quantile(samples, 0.05)).to be_within(0.18).of(-3.934)
      expect(quantile(samples, 0.5)).to be_within(0.04).of(-0.134)
      expect(quantile(samples, 0.95)).to be_within(0.18).of(2.254)
    end
  end
end
