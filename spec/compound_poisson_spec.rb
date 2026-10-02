require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::CompoundPoisson do
  let(:engine) { Distribution::CompoundPoisson::Ruby_ }

  describe '.pmf_table' do
    it 'degenerates to a point mass at zero when lambda is zero' do
      res = engine.pmf_table(0, { 0 => 0.3, 2 => 0.7 }, 0..5)
      expect(res.pmf[0]).to eq 1.0
      (1..5).each { |k| expect(res.pmf[k]).to eq 0.0 }
      (0..5).each { |k| expect(res.cdf[k]).to eq 1.0 }
      expect(res.omitted_mass).to eq 0.0
    end

    it 'allows zero-cost requests' do
      # Poisson thinning: with costs in {0, 1}, the total is
      # Poisson(lambda * P(X = 1)).
      res = engine.pmf_table(2.0, { 0 => 0.5, 1 => 0.5 }, 0..10)
      (0..10).each do |k|
        expect(res.pmf[k]).to be_within(1e-12).of(Distribution::Poisson::Ruby_.pdf(k, 1.0))
      end
    end

    it 'reduces to the Poisson distribution when the cost is always one' do
      res = engine.pmf_table(3.5, { 1 => 1.0 }, 0..15)
      (0..15).each do |k|
        expect(res.pmf[k]).to be_within(1e-12).of(Distribution::Poisson::Ruby_.pdf(k, 3.5))
        expect(res.cdf[k]).to be_within(1e-12).of(Distribution::Poisson::Ruby_.cdf(k, 3.5))
      end
    end

    it 'assigns exactly zero mass to odd totals when every cost is even' do
      res = engine.pmf_table(2.5, { 2 => 0.25, 4 => 0.75 }, 0..31)
      (0..31).each do |k|
        expect(res.pmf[k]).to eq(0.0) if k.odd?
      end
      expect(res.covered_mass).to be > 0.99
    end

    it 'matches a brute-force convolution over the request count' do
      table = { 0 => 0.2, 1 => 0.3, 2 => 0.5 }
      lambda_value = 1.5
      max = 12
      expected = Array.new(max + 1, 0.0)
      conv = [1.0]
      (0..40).each do |n|
        p_n = Distribution::Poisson::Ruby_.pdf(n, lambda_value)
        (0..max).each { |s| expected[s] += p_n * (conv[s] || 0.0) }
        nxt = Array.new(max + 1, 0.0)
        (0..max).each do |s|
          next if conv[s].nil? || conv[s].zero?
          table.each do |cost, prob|
            nxt[s + cost] += conv[s] * prob if s + cost <= max
          end
        end
        conv = nxt
      end
      res = engine.pmf_table(lambda_value, table, 0..max)
      (0..max).each { |k| expect(res.pmf[k]).to be_within(1e-12).of(expected[k]) }
    end

    it 'reports omitted mass consistent with the returned pmf and cdf' do
      table = { 1 => 0.5, 3 => 0.5 }
      res = engine.pmf_table(4.0, table, 2..10)
      expect(res.covered_mass + res.omitted_mass).to be_within(1e-12).of(1.0)
      expect(res.omitted_mass).to be_within(1e-15).of(res.lower_mass + res.upper_mass)
      expect(res.lower_mass).to be_within(1e-12).of(engine.cdf(1, 4.0, table))
      expect(res.upper_mass).to be > 0
      expect(res.cdf[10]).to be_within(1e-12).of(1.0 - res.upper_mass)
    end

    it 'accepts a bare integer as the range 0..max' do
      res = engine.pmf_table(1.0, { 1 => 1.0 }, 4)
      expect(res.pmf.keys).to eq((0..4).to_a)
      expect(res.cdf.keys).to eq((0..4).to_a)
    end

    it 'validates its arguments' do
      expect { engine.pmf_table(-1, { 1 => 1.0 }, 0..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, {}, 0..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, { 1 => 0.5 }, 0..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, { -1 => 1.0 }, 0..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, { 1.5 => 1.0 }, 0..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, { 1 => -0.5, 2 => 1.5 }, 0..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, { 1 => 1.0 }, 5..3) }.to raise_error(ArgumentError)
      expect { engine.pmf_table(1, { 1 => 1.0 }, -2..3) }.to raise_error(ArgumentError)
    end
  end

  describe '.pdf and .cdf' do
    it 'are consistent with pmf_table' do
      table = { 0 => 0.1, 1 => 0.4, 5 => 0.5 }
      res = engine.pmf_table(2.3, table, 0..8)
      (0..8).each do |k|
        expect(engine.pdf(k, 2.3, table)).to be_within(1e-14).of(res.pmf[k])
        expect(engine.cdf(k, 2.3, table)).to be_within(1e-14).of(res.cdf[k])
      end
    end

    it 'returns the degenerate distribution for lambda zero' do
      expect(engine.pdf(0, 0, { 1 => 1.0 })).to eq 1.0
      expect(engine.pdf(3, 0, { 1 => 1.0 })).to eq 0.0
      expect(engine.cdf(3, 0, { 1 => 1.0 })).to eq 1.0
    end
  end

  it 'is exposed through the Distribution::CompoundPoisson facade' do
    expect(Distribution::CompoundPoisson.pdf(2, 1.0, { 1 => 1.0 }))
      .to be_within(1e-12).of(Distribution::Poisson::Ruby_.pdf(2, 1.0))
    res = Distribution::CompoundPoisson.pmf_table(1.0, { 1 => 1.0 }, 0..3)
    expect(res).to be_a(Distribution::CompoundPoisson::Result)
    expect(res[0]).to be_within(1e-12).of(Math.exp(-1))
  end
end
