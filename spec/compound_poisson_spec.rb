require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::CompoundPoisson do
  shared_examples_for 'compound poisson engine' do
    it 'reduces to Poisson when every request costs one' do
      lambda_val = 2.5
      result = @engine.distribution(lambda_val, { 1 => 1.0 }, 0..10)

      0.upto(10) do |cost|
        expected = Distribution::Poisson::Ruby_.pdf(cost, lambda_val)
        expect(result.pmf[cost]).to be_within(1e-13).of(expected)
        expect(result.cdf[cost]).to be_within(1e-13).of(Distribution::Poisson::Ruby_.cdf(cost, lambda_val))
        expect(@engine.pdf(cost, lambda_val, { 1 => 1.0 })).to be_within(1e-13).of(expected)
        expect(@engine.pmf(cost, lambda_val, { 1 => 1.0 })).to be_within(1e-13).of(expected)
      end
    end

    it 'puts all probability at zero when lambda is zero' do
      result = @engine.distribution(0, { 0 => 0.3, 1 => 0.4, 2 => 0.3 }, -1..3)

      expect(result.pmf[0]).to eq(1.0)
      (-1..3).reject(&:zero?).each { |cost| expect(result.pmf[cost]).to eq(0.0) }
      expect(result.cdf[-1]).to eq(0.0)
      0.upto(3) { |cost| expect(result.cdf[cost]).to eq(1.0) }
      expect(result.range_mass).to be_within(1e-15).of(1.0)
      expect(result.omitted_mass[:total]).to be_within(1e-15).of(0.0)
    end

    it 'supports zero as a possible single-request cost' do
      lambda_val = 4.0
      cost_pmf = { 0 => 0.7, 1 => 0.3 }
      result = @engine.distribution(lambda_val, cost_pmf, 0..8)

      0.upto(8) do |cost|
        expected = Distribution::Poisson::Ruby_.pdf(cost, 1.2)
        expect(result.pmf[cost]).to be_within(1e-13).of(expected)
      end
    end

    it 'keeps odd totals impossible when every request has even cost' do
      result = @engine.distribution(2.0, { 2 => 1.0 }, 0..6)

      [1, 3, 5].each { |cost| expect(result.pmf[cost]).to eq(0.0) }
      expect(result.pmf[0]).to be_within(1e-13).of(Math.exp(-2.0))
      expect(result.pmf[2]).to be_within(1e-13).of(2.0 * Math.exp(-2.0))
      expect(result.pmf[4]).to be_within(1e-13).of(2.0 * Math.exp(-2.0))
      expect(result.pmf[6]).to be_within(1e-13).of((4.0 / 3.0) * Math.exp(-2.0))
    end

    it 'reports mass inside, below, and above a bounded range' do
      lambda_val = 3.0
      cost_pmf = { 1 => 0.6, 2 => 0.4 }
      result = @engine.distribution(lambda_val, cost_pmf, 2..4)
      full = @engine.distribution(lambda_val, cost_pmf, 0..40)

      expected_below = full.cdf[1]
      expected_range = (2..4).sum { |cost| full.pmf[cost] }
      expected_above = 1.0 - full.cdf[4]

      expect(result.range_mass).to be_within(1e-12).of(expected_range)
      expect(result.omitted_mass[:below]).to be_within(1e-12).of(expected_below)
      expect(result.omitted_mass[:above]).to be_within(1e-12).of(expected_above)
      expect(result.omitted_mass[:total]).to be_within(1e-12).of(expected_below + expected_above)
      expect(result.range_mass + result.omitted_mass[:total]).to be_within(1e-11).of(1.0)
    end

    it 'accepts separate range endpoints and an array cost table' do
      result = @engine.distribution(1.0, [0.25, 0.75], 1, 3)
      other = @engine.distribution(1.0, { 0 => 0.25, 1 => 0.75 }, 1..3)

      expect(result.to_h).to eq(other.to_h)
    end

    it 'reports all selected mass as above a wholly negative range' do
      result = @engine.distribution(2.0, { 1 => 1.0 }, -3..-1)

      expect(result.range_mass).to eq(0.0)
      expect(result.omitted_mass[:below]).to eq(0.0)
      expect(result.omitted_mass[:above]).to eq(1.0)
      expect(result.omitted_mass[:total]).to eq(1.0)
    end

    it 'rejects invalid probability tables, lambda, and ranges' do
      expect { @engine.distribution(1, { -1 => 1 }, 0..2) }.to raise_error(ArgumentError)
      expect { @engine.distribution(1, { 1 => -0.5 }, 0..2) }.to raise_error(ArgumentError)
      expect { @engine.distribution(1, { 1 => 0.5 }, 0..2) }.to raise_error(ArgumentError)
      expect { @engine.distribution(-1, { 1 => 1 }, 0..2) }.to raise_error(ArgumentError)
      expect { @engine.distribution(1, { 1 => 1 }, 3..1) }.to raise_error(ArgumentError)
    end
  end

  describe 'singleton' do
    before { @engine = Distribution::CompoundPoisson }
    it_should_behave_like 'compound poisson engine'
  end

  describe Distribution::CompoundPoisson::Ruby_ do
    before { @engine = Distribution::CompoundPoisson::Ruby_ }
    it_should_behave_like 'compound poisson engine'
  end
end
