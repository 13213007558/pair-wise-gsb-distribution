require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::Stable do
  def quantile(samples, probability)
    sorted = samples.sort
    sorted[(probability * (sorted.size - 1)).round]
  end

  describe 'S0 and S1 conversion' do
    it 'changes only location away from alpha=1' do
      converted = described_class.parameters_from_s1(
        alpha: 0.5, beta: 1.0, delta: 0.0, gamma: 2.0
      )

      expect(converted[:alpha]).to eq(0.5)
      expect(converted[:beta]).to eq(1.0)
      expect(converted[:gamma]).to eq(2.0)
      expect(converted[:delta]).to be_within(1e-12).of(2.0)
    end

    it 'uses the alpha=1 scale-linked location correction' do
      converted = described_class.parameters_from_s1(
        alpha: 1.0, beta: 0.5, delta: 0.0, gamma: 3.0
      )
      expected = 2.0 * 0.5 * 3.0 * Math.log(3.0) / Math::PI

      expect(converted[:delta]).to be_within(1e-12).of(expected)
    end

    it 'is invertible' do
      original = { alpha: 1.25, beta: -0.4, delta: 1.7, gamma: 0.8 }
      from_s1 = described_class.parameters_from_s1(**original)
      restored = described_class.parameters_to_s1(**from_s1)

      restored.each do |name, value|
        expect(value).to be_within(1e-12).of(original[name])
      end
    end
  end

  describe 'PDF and CDF' do
    it 'matches the normal special case for alpha=2' do
      x = 1.3
      normal_pdf = Distribution::Normal.pdf(x / Math.sqrt(2.0)) / Math.sqrt(2.0)
      normal_cdf = Distribution::Normal.cdf(x / Math.sqrt(2.0))

      expect(described_class.pdf(x, alpha: 2, beta: 1)).to\
        be_within(1e-10).of(normal_pdf)
      expect(described_class.cdf(x, alpha: 2, beta: -1)).to\
        be_within(1e-10).of(normal_cdf)
    end

    it 'matches the standard Cauchy special case at alpha=1, beta=0' do
      x = -0.7

      expect(described_class.pdf(x, alpha: 1)).to\
        be_within(1e-10).of(1.0 / (Math::PI * (1.0 + x * x)))
      expect(described_class.cdf(x, alpha: 1)).to\
        be_within(1e-10).of(0.5 + Math.atan(x) / Math::PI)
    end

    it 'applies S0 location and scale' do
      x = 1.8
      z = (x - 0.4) / 1.25
      base_pdf = described_class.pdf(z, alpha: 1.4, beta: 0.3,
                                     tolerance: 1e-11)
      base_cdf = described_class.cdf(z, alpha: 1.4, beta: 0.3,
                                     tolerance: 1e-11)

      expect(described_class.pdf(x, alpha: 1.4, beta: 0.3, delta: 0.4,
                                 gamma: 1.25, tolerance: 1e-11)).to\
        be_within(1e-9).of(base_pdf / 1.25)
      expect(described_class.cdf(x, alpha: 1.4, beta: 0.3, delta: 0.4,
                                 gamma: 1.25, tolerance: 1e-11)).to\
        be_within(1e-9).of(base_cdf)
    end

    it 'matches direct S0 characteristic-function values near alpha=1' do
      cases = [
        [0.99, 0.5, 0.0, 0.292787, 0.437141],
        [1.00, 0.5, 0.0, 0.292520, 0.43751],
        [1.01, 0.5, 0.0, 0.292259, 0.437887]
      ]

      cases.each do |alpha, beta, x, expected_pdf, expected_cdf|
        expect(described_class.pdf(x, alpha: alpha, beta: beta,
                                   tolerance: 1e-9)).to\
          be_within(2e-5).of(expected_pdf)
        expect(described_class.cdf(x, alpha: alpha, beta: beta,
                                   tolerance: 1e-9)).to\
          be_within(2e-5).of(expected_cdf)
      end
    end

    it 'does not jump when alpha crosses one with nonzero skew' do
      cdf_at_med = lambda do |alpha|
        described_class.cdf(0.0, alpha: alpha, beta: 0.5, tolerance: 1e-9)
      end
      pdf_at_med = lambda do |alpha|
        described_class.pdf(0.0, alpha: alpha, beta: 0.5, tolerance: 1e-9)
      end

      expect((cdf_at_med.call(0.999) - cdf_at_med.call(1.001)).abs).to be < 0.001
      expect((pdf_at_med.call(0.999) - pdf_at_med.call(1.001)).abs).to be < 0.001
    end

    it 'returns diagnostics for successful numerical integration' do
      result = described_class.cdf(
        0.0, alpha: 1.5, beta: 0.4, diagnostics: true, tolerance: 1e-8
      )

      expect(result).to include(converged: true,
                                method: a_kind_of(Symbol),
                                error: a_kind_of(Float),
                                requested_tolerance: 1e-8)
      expect(result[:error]).to be <= result[:requested_tolerance]
    end

    it 'reports non-convergence for an unattainable tolerance' do
      result = described_class.pdf(
        0.0, alpha: 1.5, beta: 0.4, diagnostics: true, tolerance: 1e-200
      )

      expect(result[:converged]).to eq(false)
      expect(result[:error]).to eq(Float::INFINITY)
      expect(result[:message]).to include('did not reach requested tolerance')
      expect {
        described_class.pdf(0.0, alpha: 1.5, beta: 0.4, tolerance: 1e-200)
      }.to raise_error(RuntimeError, /did not reach requested tolerance/)
    end

    it 'validates parameter ranges' do
      expect {
        described_class.pdf(0, alpha: 0)
      }.to raise_error(ArgumentError, /alpha/)
      expect {
        described_class.pdf(0, alpha: 1, beta: 1.1)
      }.to raise_error(ArgumentError, /beta/)
      expect {
        described_class.pdf(0, alpha: 1, gamma: 0)
      }.to raise_error(ArgumentError, /gamma/)
    end

    it 'handles the one-sided support boundary for alpha below one' do
      boundary = Math.tan(Math::PI * 0.5 / 2.0)

      expect(described_class.pdf(-boundary, alpha: 0.5, beta: 1)).to eq(0)
      expect(described_class.cdf(-boundary, alpha: 0.5, beta: 1)).to eq(0)
      expect(described_class.pdf(boundary, alpha: 0.5, beta: -1)).to eq(0)
      expect(described_class.cdf(boundary, alpha: 0.5, beta: -1)).to eq(1)
    end
  end

  describe 'sampling' do
    it 'reproduces fixed-seed normal samples' do
      random = Random.new(20261006)
      samples = Array.new(10_000) do
        described_class.rng(alpha: 2, source: random).call
      end

      expect(quantile(samples, 0.025)).to be_within(0.25).of(-1.96 * Math.sqrt(2))
      expect(quantile(samples, 0.5)).to be_within(0.12).of(0.0)
      expect(quantile(samples, 0.975)).to be_within(0.25).of(1.96 * Math.sqrt(2))
    end

    it 'reproduces fixed-seed Cauchy samples' do
      random = Random.new(20261006)
      samples = Array.new(10_000) do
        described_class.rng(alpha: 1, beta: 0, source: random).call
      end

      expect(quantile(samples, 0.25)).to be_within(0.18).of(-1.0)
      expect(quantile(samples, 0.5)).to be_within(0.12).of(0.0)
      expect(quantile(samples, 0.75)).to be_within(0.18).of(1.0)
    end

    it 'remains continuous through alpha=1 for fixed-seed skew samples' do
      medians = [0.999, 1.0, 1.001].map do |alpha|
        random = Random.new(20261006)
        samples = Array.new(10_000) do
          described_class.rng(alpha: alpha, beta: 0.5, source: random).call
        end
        quantile(samples, 0.5)
      end

      expect((medians[0] - medians[1]).abs).to be < 0.25
      expect((medians[1] - medians[2]).abs).to be < 0.25
    end

    it 'accepts a callable uniform source and applies S0 location and scale' do
      random_one = Random.new(42)
      source = -> { random_one.rand }
      base = described_class.rng(alpha: 1, beta: 0.5, source: source)
      random_two = Random.new(42)
      shifted_source = -> { random_two.rand }
      shifted = described_class.rng(alpha: 1, beta: 0.5, delta: 2.0,
                                   gamma: 3.0, source: shifted_source)

      10.times do
        expect(shifted.call).to be_within(1e-12).of(base.call * 3.0 + 2.0)
      end
    end
  end
end
