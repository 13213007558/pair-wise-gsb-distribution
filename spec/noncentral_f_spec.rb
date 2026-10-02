require 'json'
require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::NoncentralF::Ruby_ do
  let(:engine) { described_class }
  let(:references) do
    path = File.expand_path('fixtures/noncentral_f_references.json', __dir__)
    JSON.parse(File.read(path), symbolize_names: true)[:cases]
  end

  it 'matches independent high precision PDF, CDF and survival references' do
    references.each do |reference|
      args = reference.values_at(:d1, :d2, :ncp, :x)
      d1, d2, ncp, x = args
      expect(engine.pdf(x, d1, d2, ncp)).to be_within([reference[:pdf].abs * 1e-10, 1e-14].max).of(reference[:pdf])
      expect(engine.cdf(x, d1, d2, ncp)).to be_within([reference[:cdf].abs * 1e-10, 1e-12].max).of(reference[:cdf])
      expect(engine.survival(x, d1, d2, ncp)).to be_within([reference[:survival].abs * 1e-10, 1e-12].max).of(reference[:survival])
    end
  end

  it 'matches independent inverse CDF references near 0, at the center, and near 1' do
    references.select { |reference| reference.key?(:probability) }.each do |reference|
      expected = reference[:x]
      relative_tolerance = reference[:probability] == 0.5 ? 1e-8 : 5e-4
      tolerance = [expected.abs * relative_tolerance, 1e-12].max
      actual = engine.quantile(reference[:probability], reference[:d1], reference[:d2], reference[:ncp])
      expect(actual).to be_within(tolerance).of(expected)
    end
  end

  it 'delegates lambda=0 to the existing central F implementation' do
    x = 0.75
    d1 = 5
    d2 = 10
    expect(engine.pdf(x, d1, d2, 0)).to eq(Distribution::F::Ruby_.pdf(x, d1, d2))
    expect(engine.cdf(x, d1, d2, 0)).to eq(Distribution::F::Ruby_.cdf(x, d1, d2))
    expect(engine.survival(x, d1, d2, 0)).to eq(Distribution::F::Ruby_.q(x, d1, d2))
    expect(engine.quantile(0.3, d1, d2, 0)).to eq(Distribution::F::Ruby_.quantile(0.3, d1, d2))
  end

  it 'validates degrees of freedom, noncentrality and numeric controls' do
    expect { engine.pdf(1, 0, 10, 2) }.to raise_error(ArgumentError)
    expect { engine.pdf(1, 5, 0, 2) }.to raise_error(ArgumentError)
    expect { engine.pdf(1, 5, 10, -1) }.to raise_error(ArgumentError)
    expect { engine.pdf(1, 5, 10, 2, tolerance: 0) }.to raise_error(ArgumentError)
    expect { engine.pdf(1, 5, 10, 2, max_iterations: 0) }.to raise_error(ArgumentError)
    expect { engine.quantile(1.1, 5, 10, 2) }.to raise_error(ArgumentError)
  end

  it 'reports a non-completed result when the PDF calculation budget is exhausted' do
    begin
      engine.pdf(0.5, 5, 10, 2, max_iterations: 2)
    rescue Distribution::NoncentralFConvergenceError => error
      expect(error.result).to be_a(Float)
      expect(error.error).to be > 1e-12
      expect(error.iterations).to eq(2)
      expect(error.budget).to eq(2)
    else
      raise 'Expected NoncentralFConvergenceError'
    end
  end

  it 'reports a non-completed result when the CDF calculation budget is exhausted' do
    begin
      engine.cdf(0.5, 5, 10, 2, max_iterations: 2)
    rescue Distribution::NoncentralFConvergenceError => error
      expect(error.result).to be_a(Float)
      expect(error.error).to be > 1e-12
      expect(error.iterations).to eq(2)
      expect(error.budget).to eq(2)
    else
      raise 'Expected NoncentralFConvergenceError'
    end
  end

  it 'exposes the same noncentral interface through Distribution::F' do
    expect(Distribution::F.pdf(2, 5, 10, 2)).to eq(engine.pdf(2, 5, 10, 2))
    expect(Distribution::F.cdf(2, 5, 10, 2)).to eq(engine.cdf(2, 5, 10, 2))
    expect(Distribution::F.survival(2, 5, 10, 2)).to eq(engine.survival(2, 5, 10, 2))
    expect(Distribution::F.quantile(0.5, 5, 10, 2)).to be_within(1e-10).of(engine.quantile(0.5, 5, 10, 2))
  end
end
