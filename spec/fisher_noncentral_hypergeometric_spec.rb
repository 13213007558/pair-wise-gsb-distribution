require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::FisherNoncentralHypergeometric do
  let(:engine) { described_class }

  def support(m1, m2, n)
    lower = [0, n - m2].max
    upper = [n, m1].min
    lower..upper
  end

  def exact_distribution(m1, m2, n, odds)
    weights = Hash.new(0)
    total_weight = Rational(0)

    (0...(1 << (m1 + m2))).each do |mask|
      next unless mask.to_s(2).count('1') == n

      type1_count = (0...m1).count { |i| mask[i] == 1 }
      weight = Rational(odds) ** type1_count
      weights[type1_count] += weight
      total_weight += weight
    end

    weights.each_with_object({}) do |(k, weight), probabilities|
      probabilities[k] = weight / total_weight
    end
  end

  describe Distribution::FisherNoncentralHypergeometric::Ruby_ do
    it 'provides the complete discrete distribution interface' do
      expect(engine).to respond_to(:pmf)
      expect(engine).to respond_to(:pdf)
      expect(engine).to respond_to(:cdf)
      expect(engine).to respond_to(:sf)
      expect(engine).to respond_to(:quantile)
      expect(engine).to respond_to(:p_value)
      expect(engine).to respond_to(:mean)
      expect(engine).to respond_to(:variance)
    end

    it 'matches an exhaustive enumeration of weighted subsets' do
      odds_values = [Rational(1), Rational(2), Rational(1, 2), Rational(7, 3), Rational(3, 7)]

      (0..3).each do |m1|
        (0..3).each do |m2|
          (0..(m1 + m2)).each do |n|
            range = support(m1, m2, n)
            odds_values.each do |odds|
              exact = exact_distribution(m1, m2, n, odds)
              exact_mean = range.reduce(Rational(0)) do |sum, k|
                sum + k * exact.fetch(k)
              end
              exact_variance = range.reduce(Rational(0)) do |sum, k|
                sum + (k - exact_mean) ** 2 * exact.fetch(k)
              end

              range.each do |k|
                exact_cdf = range.reduce(Rational(0)) do |sum, j|
                  j <= k ? sum + exact.fetch(j) : sum
                end
                exact_sf = range.reduce(Rational(0)) do |sum, j|
                  j > k ? sum + exact.fetch(j) : sum
                end

                expect(engine.pmf(k, m1, m2, n, odds)).to be_within(1e-14).of(exact.fetch(k).to_f)
                expect(engine.cdf(k, m1, m2, n, odds)).to be_within(1e-14).of(exact_cdf.to_f)
                expect(engine.sf(k, m1, m2, n, odds)).to be_within(1e-14).of(exact_sf.to_f)
              end

              expect(engine.mean(m1, m2, n, odds)).to be_within(1e-13).of(exact_mean.to_f)
              expect(engine.variance(m1, m2, n, odds)).to be_within(1e-13).of(exact_variance.to_f)

              [0.0, 0.25, 0.5, 0.75, 1.0].each do |p|
                q = engine.quantile(p, m1, m2, n, odds)
                expect(range).to cover(q)
                expect(engine.cdf(q, m1, m2, n, odds)).to be >= p - 1e-13
                if p > 0 && q > range.first
                  expect(engine.cdf(q - 1, m1, m2, n, odds) - p).to be <= 1e-13
                end
              end
            end
          end
        end
      end
    end

    it 'handles zero draws and drawing the complete population' do
      [[5, 7, 0], [5, 7, 12], [0, 0, 0]].each do |m1, m2, n|
        range = support(m1, m2, n)

        expect(engine.pmf(range.first, m1, m2, n, 3)).to eq(1.0)
        expect(engine.cdf(range.first, m1, m2, n, 3)).to eq(1.0)
        expect(engine.cdf(range.first - 1, m1, m2, n, 3)).to eq(0.0)
        expect(engine.sf(range.first, m1, m2, n, 3)).to eq(0.0)
        expect(engine.sf(range.first - 1, m1, m2, n, 3)).to eq(1.0)
        expect(engine.quantile(0, m1, m2, n, 3)).to eq(range.first)
        expect(engine.quantile(1, m1, m2, n, 3)).to eq(range.last)
        expect(engine.mean(m1, m2, n, 3)).to eq(range.first.to_f)
        expect(engine.variance(m1, m2, n, 3)).to eq(0.0)
      end
    end

    it 'reduces to the ordinary hypergeometric distribution when odds is one' do
      cases = [[4, 6, 5], [0, 5, 2], [5, 0, 3], [3, 3, 0], [3, 3, 6]]

      cases.each do |m1, m2, n|
        total = m1 + m2
        range = support(m1, m2, n)

        range.each do |k|
          expect(engine.pmf(k, m1, m2, n, 1)).to be_within(1e-14).of(
            Distribution::Hypergeometric.pdf(k, m1, n, total).to_f
          )
          expect(engine.cdf(k, m1, m2, n, 1)).to be_within(1e-13).of(
            Distribution::Hypergeometric.cdf(k, m1, n, total).to_f
          )
        end

        expect(engine.mean(m1, m2, n, 1)).to be_within(1e-13).of(n * m1.to_f / total) if total > 0
        if total > 1
          variance = n * m1 * m2 * (total - n).to_f / (total * total * (total - 1))
          expect(engine.variance(m1, m2, n, 1)).to be_within(1e-12).of(variance)
        end
      end
    end

    it 'preserves the symmetry relation after exchanging the two types' do
      cases = [[3, 5, 4], [5, 2, 3], [4, 4, 0], [4, 4, 8], [1, 3, 2]]
      odds_values = [Rational(3), Rational(1, 3), 0, Float::INFINITY]

      cases.each do |m1, m2, n|
        odds_values.each do |odds|
          reciprocal = if odds == 0
                         Float::INFINITY
                       elsif odds.infinite?
                         0
                       else
                         Rational(1) / odds
                       end
          range = support(m1, m2, n)

          range.each do |k|
            swapped_k = n - k

            expect(engine.pmf(k, m1, m2, n, odds)).to be_within(1e-14).of(
              engine.pmf(swapped_k, m2, m1, n, reciprocal)
            )
            expect(engine.cdf(k, m1, m2, n, odds)).to be_within(1e-14).of(
              engine.sf(swapped_k - 1, m2, m1, n, reciprocal)
            )
            expect(engine.sf(k, m1, m2, n, odds)).to be_within(1e-14).of(
              engine.cdf(swapped_k - 1, m2, m1, n, reciprocal)
            )
          end

          expect(engine.mean(m1, m2, n, odds) + engine.mean(m2, m1, n, reciprocal)).to be_within(1e-12).of(n)
          expect(engine.variance(m1, m2, n, odds)).to be_within(1e-12).of(
            engine.variance(m2, m1, n, reciprocal)
          )
        end
      end
    end

    it 'keeps far-tail probabilities and moments finite for extreme odds' do
      expect(engine.pmf(4, 5, 5, 6, 1e100)).to be_within(1e-111).of(1e-99)
      expect(engine.cdf(4, 5, 5, 6, 1e100)).to be_within(1e-111).of(1e-99)
      expect(engine.variance(5, 5, 6, 1e100)).to be_within(1e-111).of(1e-99)
      expect(engine.mean(5, 5, 6, 1e100)).to be_within(1e-12).of(5)

      expect(engine.pmf(2, 5, 5, 6, 1e-100)).to be_within(1e-111).of(1e-99)
      expect(engine.sf(1, 5, 5, 6, 1e-100)).to be_within(1e-111).of(1e-99)
      expect(engine.variance(5, 5, 6, 1e-100)).to be_within(1e-111).of(1e-99)
      expect(engine.mean(5, 5, 6, 1e-100)).to be_within(1e-12).of(1)
    end

    it 'handles zero and infinite odds at the fixed-margin boundaries' do
      range = support(3, 5, 4)

      expect(engine.pmf(range.first, 3, 5, 4, 0)).to eq(1.0)
      expect(engine.pmf(range.last, 3, 5, 4, 0)).to eq(0.0)
      expect(engine.mean(3, 5, 4, 0)).to eq(range.first.to_f)
      expect(engine.variance(3, 5, 4, 0)).to eq(0.0)

      expect(engine.pmf(range.last, 3, 5, 4, Float::INFINITY)).to eq(1.0)
      expect(engine.pmf(range.first, 3, 5, 4, Float::INFINITY)).to eq(0.0)
      expect(engine.mean(3, 5, 4, Float::INFINITY)).to eq(range.last.to_f)
      expect(engine.variance(3, 5, 4, Float::INFINITY)).to eq(0.0)
    end

    it 'stays finite without materializing large binomial coefficients as floats' do
      m1 = 20_000
      m2 = 30_000
      n = 10_000
      k = 4_000

      aggregate_failures do
        expect(engine.pmf(k, m1, m2, n, 2.5)).to be_between(0, 1)
        expect(engine.cdf(k, m1, m2, n, 2.5) + engine.sf(k, m1, m2, n, 2.5)).to be_within(1e-12).of(1)
        expect(engine.quantile(0.5, m1, m2, n, 2.5)).to be_between(0, n)
        expect(engine.mean(m1, m2, n, 2.5)).to be_finite
        expect(engine.variance(m1, m2, n, 2.5)).to be_finite
      end
    end

    it 'rejects invalid parameters' do
      expect { engine.pmf(0, -1, 2, 1, 1) }.to raise_error(ArgumentError)
      expect { engine.pmf(0, 2, -1, 1, 1) }.to raise_error(ArgumentError)
      expect { engine.pmf(0, 2, 1, 4, 1) }.to raise_error(ArgumentError)
      expect { engine.pmf(0, 2, 1, 1, -1) }.to raise_error(ArgumentError)
      expect { engine.pmf(0.5, 2, 1, 1, 1) }.to raise_error(ArgumentError)
      expect { engine.quantile(1.5, 2, 1, 1, 1) }.to raise_error(ArgumentError)
    end
  end
end
