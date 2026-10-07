require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::DirichletMultinomial do
  # All compositions of n into k non-negative parts (the full support).
  def compositions(n, k)
    return [[n]] if k == 1

    (0..n).flat_map { |x| compositions(n - x, k - 1).map { |rest| [x] + rest } }
  end

  def rising(a, k)
    (0...k).inject(1.to_r) { |ac, j| ac * (a + j) }
  end

  def factorial(n)
    (1..n).inject(1, :*)
  end

  # Exact reference pmf: n!/prod(x_i!) * prod((a_i)^(x_i)) / (a0)^(n)
  # with Pochhammer symbols, computed with rationals.
  def reference_pmf(counts, alphas)
    n = counts.inject(0, :+)
    a0 = alphas.inject(0, :+)
    mult = counts.inject(factorial(n).to_r) { |ac, x| ac / factorial(x) }
    ratio = alphas.each_with_index.inject(1.to_r) do |ac, (a, i)|
      ac * rising(a, counts[i])
    end / rising(a0, n)
    (mult * ratio).to_f
  end

  let(:alphas) { [3.quo(2), 1.quo(2), 2] }
  let(:n) { 6 }
  let(:dm) { Distribution::DirichletMultinomial.new(n, alphas) }

  describe '.new' do
    it 'rejects invalid totals and concentrations' do
      expect { Distribution::DirichletMultinomial.new(-1, [1, 1]) }.to raise_error(ArgumentError)
      expect { Distribution::DirichletMultinomial.new(2.5, [1, 1]) }.to raise_error(ArgumentError)
      expect { Distribution::DirichletMultinomial.new(5, [1]) }.to raise_error(ArgumentError)
      expect { Distribution::DirichletMultinomial.new(5, [1, 0]) }.to raise_error(ArgumentError)
      expect { Distribution::DirichletMultinomial.new(5, [1, -2]) }.to raise_error(ArgumentError)
    end
  end

  describe '#log_pmf' do
    it 'matches the exact Pochhammer reference on the whole support' do
      compositions(n, 3).each do |counts|
        exp = Math.log(reference_pmf(counts, alphas))
        obs = dm.log_pmf(counts)
        expect(obs).to be_within(1e-10).of(exp), "for counts #{counts}"
      end
    end

    it 'returns -Infinity outside the support' do
      expect(dm.log_pmf([n + 1, 0, 0])).to eq(-Float::INFINITY)
      expect(dm.log_pmf([n - 1, 0, 0])).to eq(-Float::INFINITY)
      expect(dm.log_pmf([-1, n + 1, 0])).to eq(-Float::INFINITY)
      expect(dm.log_pmf([1.5, n - 1.5, 0])).to eq(-Float::INFINITY)
    end

    it 'raises on a wrongly sized count vector' do
      expect { dm.log_pmf([n, 0]) }.to raise_error(ArgumentError)
      expect { dm.log_pmf([n, 0, 0, 0]) }.to raise_error(ArgumentError)
    end

    it 'enumerates to a total probability of 1' do
      total = compositions(n, 3).inject(0.0) { |ac, c| ac + dm.pmf(c) }
      expect(total).to be_within(1e-12).of(1.0)
    end
  end

  describe '#marginal' do
    it 'matches the joint pmf summed over the pooled categories' do
      (0..2).each do |i|
        margin = dm.marginal(i)
        expect(margin.n).to eq(n)
        expect(margin.k).to eq(2)
        (0..n).each do |x|
          exp = compositions(n, 3).select { |c| c[i] == x }
                                  .inject(0.0) { |ac, c| ac + dm.pmf(c) }
          expect(margin.pmf([x, n - x])).to be_within(1e-12).of(exp)
        end
      end
    end

    it 'pools several categories consistently with merge_categories' do
      dm4 = Distribution::DirichletMultinomial.new(5, [0.5, 1.0, 1.5, 2.0])
      margin = dm4.marginal(0, 2)
      merged = dm4.merge_categories([0], [2], [1, 3])
      compositions(5, 3).each do |c|
        expect(margin.pmf(c)).to be_within(1e-12).of(merged.pmf(c))
      end
    end

    it 'rejects invalid index sets' do
      expect { dm.marginal(0, 1, 2) }.to raise_error(ArgumentError)
      expect { dm.marginal(0, 0) }.to raise_error(ArgumentError)
      expect { dm.marginal(3) }.to raise_error(ArgumentError)
    end
  end

  describe '#conditional' do
    let(:known) { { 0 => 2 } }
    let(:cond) { dm.conditional(known) }

    it 'recomputes total and concentrations instead of reusing parent caches' do
      # Populate the parent's caches before conditioning.
      expect(dm.total_concentration).to be_within(1e-12).of(4.0)
      dm.log_pmf([2, 2, 2])

      expect(cond.n).to eq(4)
      expect(cond.alphas).to eq([0.5, 2.0])
      expect(cond.total_concentration).to be_within(1e-12).of(2.5)
    end

    it 'accepts an array with nils for unknown categories' do
      expect(dm.conditional([2, nil, nil]).alphas).to eq(cond.alphas)
      expect(dm.conditional([2, nil, nil]).n).to eq(cond.n)
    end

    it 'satisfies the chain rule against the joint and the margin' do
      p_known = dm.marginal(0).pmf([2, n - 2])
      compositions(4, 2).each do |rest|
        joint = dm.pmf([2] + rest)
        expect(cond.pmf(rest) * p_known).to be_within(1e-12).of(joint)
      end
    end

    it 'enumerates to a total probability of 1' do
      total = compositions(4, 2).inject(0.0) { |ac, c| ac + cond.pmf(c) }
      expect(total).to be_within(1e-12).of(1.0)
    end

    it 'rejects impossible or total conditioning' do
      expect { dm.conditional(0 => n + 1) }.to raise_error(ArgumentError)
      expect { dm.conditional(0 => 1, 1 => 1, 2 => 1) }.to raise_error(ArgumentError)
      expect { dm.conditional(0 => -1) }.to raise_error(ArgumentError)
    end
  end

  describe '#merge_categories' do
    let(:dm4) { Distribution::DirichletMultinomial.new(5, [0.5, 1.0, 1.5, 2.0]) }
    let(:merged) { dm4.merge_categories([0, 2], [1], [3]) }

    it 'sums the concentrations of each group' do
      expect(merged.n).to eq(5)
      expect(merged.alphas).to eq([2.0, 1.0, 2.0])
    end

    it 'matches summing the joint pmf over preimages (sample-then-merge)' do
      compositions(5, 3).each do |grouped|
        exp = compositions(5, 4).select do |c|
          [c[0] + c[2], c[1], c[3]] == grouped
        end.inject(0.0) { |ac, c| ac + dm4.pmf(c) }
        expect(merged.pmf(grouped)).to be_within(1e-12).of(exp)
      end
    end

    it 'rejects groups that do not partition the categories' do
      expect { dm4.merge_categories([0, 1], [2]) }.to raise_error(ArgumentError)
      expect { dm4.merge_categories([0, 1], [1, 2], [3]) }.to raise_error(ArgumentError)
      expect { dm4.merge_categories([0, 1], [2, 4]) }.to raise_error(ArgumentError)
    end
  end

  describe '#sample' do
    it 'returns k non-negative counts summing exactly to n' do
      rng = Random.new(42)
      100.times do
        s = dm.sample(rng)
        expect(s.size).to eq(3)
        expect(s).to all(be_a(Integer).and(be >= 0))
        expect(s.inject(0, :+)).to eq(n)
      end
    end

    it 'handles n = 0' do
      zero = Distribution::DirichletMultinomial.new(0, alphas)
      expect(zero.sample).to eq([0, 0, 0])
    end

    it 'is reproducible with a seeded rng' do
      expect(dm.sample(Random.new(7))).to eq(dm.sample(Random.new(7)))
    end
  end
end

