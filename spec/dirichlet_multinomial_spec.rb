require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::DirichletMultinomial do
  let(:distribution) { described_class.new(3, [1, 2, 3], random: 123) }

  describe '#initialize' do
    it 'copies and freezes the concentration vector' do
      alpha = [1, 2]
      dm = described_class.new(2, alpha)
      alpha[0] = 99

      expect(dm.alpha).to eq([1.0, 2.0])
      expect(dm.alpha).to be_frozen
    end

    it 'requires an integer non-negative total' do
      expect { described_class.new(-1, [1]) }.to raise_error(ArgumentError)
      expect { described_class.new(1.0, [1]) }.to raise_error(ArgumentError)
    end

    it 'requires positive finite concentrations' do
      expect { described_class.new(1, []) }.to raise_error(ArgumentError)
      expect { described_class.new(1, [0, 1]) }.to raise_error(ArgumentError)
      expect { described_class.new(1, [-1, 1]) }.to raise_error(ArgumentError)
      expect { described_class.new(1, [1, Float::NAN]) }.to raise_error(ArgumentError)
    end
  end

  describe '#log_pmf and #pmf' do
    it 'evaluates a known beta-binomial-equivalent case' do
      dm = described_class.new(2, [1, 2])

      expect(dm.pmf([2, 0])).to be_within(1e-12).of(1.0 / 6)
      expect(dm.pmf([1, 1])).to be_within(1e-12).of(1.0 / 3)
      expect(dm.pmf([0, 2])).to be_within(1e-12).of(1.0 / 2)
      expect(dm.log_pmf([1, 1])).to be_within(1e-12).of(Math.log(1.0 / 3))
    end

    it 'is normalized over the complete finite support' do
      [[1, 1], [0.5, 2.5], [1, 2, 3], [0.25, 1.25, 0.5, 2.0]].each do |alpha|
        [0, 1, 2, 5, 8].each do |n|
          dm = described_class.new(n, alpha)
          support = dm.each_count.to_a
          probability = support.inject(0.0) { |sum, counts| sum + dm.pmf(counts) }
          log_probability = Math.log(support.inject(0.0) { |sum, counts|
            sum + Math.exp(dm.log_pmf(counts))
          })

          expect(probability).to be_within(1e-12).of(1.0)
          expect(log_probability).to be_within(1e-12).of(0.0)
        end
      end
    end

    it 'returns negative infinity outside the support' do
      expect(distribution.log_pmf([1, 1, 0])).to eq(-Float::INFINITY)
      expect(distribution.log_pmf([-1, 2, 2])).to eq(-Float::INFINITY)
      expect(distribution.log_pmf([1.0, 1, 1])).to eq(-Float::INFINITY)
      expect(distribution.log_pmf([1, 1, 1, 0])).to eq(-Float::INFINITY)
    end
  end

  describe '#sample' do
    it 'returns integer counts that sum exactly to the total' do
      [0, 1, 7, 1000].each do |n|
        sample = described_class.new(n, [1, 2, 3]).sample(99)

        expect(sample.length).to eq(3)
        expect(sample).to all(be_a(Integer))
        expect(sample.sum).to eq(n)
      end
    end

    it 'uses an injected seed consistently' do
      first = described_class.new(3, [1, 2, 3], random: 7).sample
      second = described_class.new(3, [1, 2, 3]).sample(7)

      expect(first).to eq(second)
      expect(first.sum).to eq(3)
    end
  end

  describe '#each_count' do
    it 'enumerates every weak composition exactly once' do
      compositions = described_class.new(2, [1, 1, 1]).each_count.to_a

      expect(compositions).to eq([
        [0, 0, 2], [0, 1, 1], [0, 2, 0],
        [1, 0, 1], [1, 1, 0], [2, 0, 0]
      ])
    end

    it 'supports a block' do
      collected = []
      described_class.new(1, [1, 1]).each_count { |counts| collected << counts }
      expect(collected).to eq([[0, 1], [1, 0]])
    end
  end

  describe '#marginal' do
    it 'returns the beta-binomial marginal of one category' do
      marginal = distribution.marginal(0)

      expect(marginal).to be_a(described_class::BetaBinomial)
      expect(marginal.n).to eq(3)
      expect(marginal.alpha).to eq(1.0)
      expect(marginal.beta).to eq(5.0)
      expect(marginal.pmf(1)).to be_within(1e-12).of(15.0 / 56)
      expect((0..3).inject(0.0) { |sum, count| sum + marginal.pmf(count) }).to be_within(1e-12).of(1.0)
    end

    it 'matches exhaustive aggregation of one merged group' do
      dm = described_class.new(4, [1, 2, 3, 4])
      marginal = dm.marginal([0, 2])

      0.upto(4) do |merged_count|
        expected = dm.each_count.select { |counts|
          counts.values_at(0, 2).sum == merged_count
        }.inject(0.0) { |sum, counts| sum + dm.pmf(counts) }

        expect(marginal.pmf(merged_count)).to be_within(1e-12).of(expected)
      end
    end

    it 'matches exhaustive aggregation of several labeled groups' do
      dm = described_class.new(5, [1, 2, 3, 4])
      marginal = dm.marginal(first: [0, 2], second: [1])

      expect(marginal).to be_a(described_class)
      expect(marginal.alpha).to eq([4.0, 2.0, 4.0])
      expect(marginal.groups).to eq([:first, :second, :other])

      marginal.each_count do |grouped|
        expected = dm.each_count.select { |counts|
          [counts.values_at(0, 2).sum, counts[1], counts[3]] == grouped
        }.inject(0.0) { |sum, counts| sum + dm.pmf(counts) }

        expect(marginal.pmf(grouped)).to be_within(1e-12).of(expected)
      end
    end

    it 'matches exhaustive aggregation when all categories are grouped' do
      dm = described_class.new(5, [1, 2, 3, 4])
      marginal = dm.marginal([0, 2], [1, 3])

      marginal.each_count do |grouped|
        expected = dm.each_count.select { |counts|
          [counts.values_at(0, 2).sum, counts.values_at(1, 3).sum] == grouped
        }.inject(0.0) { |sum, counts| sum + dm.pmf(counts) }

        expect(marginal.pmf(grouped)).to be_within(1e-12).of(expected)
      end
    end

    it 'rejects repeated or out-of-range categories' do
      expect { distribution.marginal([0, 1], [1, 2]) }.to raise_error(ArgumentError)
      expect { distribution.marginal(3) }.to raise_error(ArgumentError)
      expect { distribution.marginal([]) }.to raise_error(ArgumentError)
      expect { distribution.marginal([0], []) }.to raise_error(ArgumentError)
    end
  end

  describe '#conditional' do
    it 'uses the remaining total and concentrations for one fixed category' do
      conditional = distribution.conditional(0 => 1)

      expect(conditional).to be_a(described_class)
      expect(conditional.n).to eq(2)
      expect(conditional.alpha).to eq([2.0, 3.0])
      expect(conditional.groups).to eq([1, 2])
      expect(distribution.n).to eq(3)
      expect(distribution.alpha).to eq([1.0, 2.0, 3.0])
    end

    it 'accepts an array with nil for random categories' do
      conditional = distribution.conditional([1, nil, 0])

      expect(conditional.n).to eq(2)
      expect(conditional.alpha).to eq([2.0])
      expect(conditional.sample).to eq([2])
    end

    it 'matches exhaustive conditional probabilities' do
      dm = described_class.new(5, [1, 2, 3, 4])
      observations = [
        { 0 => 1, 1 => 2 },
        [1, 0, nil, nil],
        [nil, 2, nil, 1]
      ]

      observations.each do |observed|
        conditional = dm.conditional(observed)
        fixed =
          if observed.is_a?(Hash)
            observed
          else
            observed.each_with_index.each_with_object({}) do |(count, index), values|
              values[index] = count unless count.nil?
            end
          end
        denominator_support = dm.each_count.select do |counts|
          fixed.all? { |index, count| counts[index] == count }
        end
        denominator = denominator_support.inject(0.0) { |sum, counts| sum + dm.pmf(counts) }

        conditional.each_count do |remaining|
          full_counts = Array.new(4)
          remaining_index = 0
          4.times do |index|
            if fixed.key?(index)
              full_counts[index] = fixed[index]
            else
              full_counts[index] = remaining[remaining_index]
              remaining_index += 1
            end
          end

          expect(conditional.pmf(remaining)).to be_within(1e-12).of(
            dm.pmf(full_counts) / denominator
          )
        end
      end
    end

    it 'supports conditioning all remaining draws' do
      conditional = distribution.conditional(0 => 2, 1 => 1)

      expect(conditional.n).to eq(0)
      expect(conditional.alpha).to eq([3.0])
      expect(conditional.sample).to eq([0])
    end

    it 'rejects impossible or incomplete conditioning' do
      expect { distribution.conditional(0 => 4) }.to raise_error(ArgumentError)
      expect { distribution.conditional(0 => -1) }.to raise_error(ArgumentError)
      expect { distribution.conditional([1, 1]) }.to raise_error(ArgumentError)
      expect { distribution.conditional(0 => 1, 1 => 1, 2 => 1) }.to raise_error(ArgumentError)
      expect { distribution.conditional(3 => 1) }.to raise_error(ArgumentError)
    end
  end

  describe described_class::BetaBinomial do
    it 'handles degenerate concentrations' do
      always_success = described_class.new(3, 2, 0)
      always_failure = described_class.new(3, 0, 2)

      expect(always_success.sample).to eq(3)
      expect(always_success.pmf(3)).to eq(1.0)
      expect(always_success.pmf(2)).to eq(0.0)
      expect(always_failure.sample).to eq(0)
      expect(always_failure.pmf(0)).to eq(1.0)
    end

    it 'is normalized and has a monotonic cdf' do
      marginal = described_class.new(8, 1.5, 3.5)
      probabilities = (0..8).map { |count| marginal.pmf(count) }

      expect(probabilities.sum).to be_within(1e-12).of(1.0)
      expect(marginal.cdf(-1)).to eq(0.0)
      expect(marginal.cdf(8)).to eq(1.0)
      cumulative = (0..7).map { |count| marginal.cdf(count) }
      expect(cumulative).to eq(cumulative.sort)
    end

    it 'returns negative infinity outside the support' do
      marginal = described_class.new(2, 1, 1)

      expect(marginal.log_pmf(-1)).to eq(-Float::INFINITY)
      expect(marginal.log_pmf(3)).to eq(-Float::INFINITY)
    end
  end

  describe 'class-level interface' do
    it 'evaluates and samples without explicitly creating a distribution' do
      log_probability = described_class.log_pmf([1, 1, 1], 3, [1, 2, 3])
      sample = described_class.sample(3, [1, 2, 3], 321)

      expect(log_probability).to eq(distribution.log_pmf([1, 1, 1]))
      expect(sample.sum).to eq(3)
    end
  end
end
