require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::PoissonBinomial do
  def enumerated_masses(probabilities)
    choices = probabilities.inject([[]]) do |accumulated, _probability|
      accumulated.product([0, 1]).map(&:flatten)
    end

    choices.each_with_object(Hash.new(0)) do |choices_by_sensor, masses|
      successes = choices_by_sensor.inject(:+) || 0
      probability = choices_by_sensor.each_with_index.inject(1) do |value, (choice, index)|
        value * (choice == 1 ? probabilities[index] : 1 - probabilities[index])
      end
      masses[successes] += probability
    end
  end

  shared_examples_for 'poisson-binomial engine' do
    it 'matches exhaustive true/false combinations' do
      probability_sets = [
        [],
        [1.quo(3)],
        [0, 1, 1.quo(4), 3.quo(4)],
        [1.quo(7), 2.quo(7), 3.quo(7), 5.quo(7), 6.quo(7)]
      ]

      probability_sets.each do |probabilities|
        expected = enumerated_masses(probabilities)

        (0..probabilities.length).each do |k|
          expect(@engine.pmf(k, probabilities)).to eq(expected[k]), "pmf(#{k}, #{probabilities.inspect})"
        end

        cumulative = 0
        (0..probabilities.length).each do |k|
          cumulative += expected[k]
          expect(@engine.cdf(k, probabilities)).to eq(cumulative), "cdf(#{k}, #{probabilities.inspect})"
          expect(@engine.right_tail(k, probabilities)).to eq(1 - cumulative + expected[k]),
            "right_tail(#{k}, #{probabilities.inspect})"
        end
      end
    end

    it 'supports floating-point probabilities' do
      probabilities = [0.1, 0.25, 0.4, 0.7, 0.9]
      expected = enumerated_masses(probabilities)

      (0..probabilities.length).each do |k|
        expect(@engine.pmf(k, probabilities)).to be_within(1e-14).of(expected[k])
      end

      expect((0..probabilities.length).map { |k| @engine.pmf(k, probabilities) }.inject(:+)).to eq(1.0)

      longer_probabilities = Array.new(60) { |index| 0.05 + (index % 10) * 0.07 }
      longer_mass = (0..longer_probabilities.length).map { |k| @engine.pmf(k, longer_probabilities) }
      expect(longer_mass.inject(:+)).to eq(1.0)
    end

    it 'treats zero and one probabilities as certain events' do
      probabilities = [0, 1, 0, 1]

      expect(@engine.pmf(2, probabilities)).to eq(1)
      expect(@engine.pmf(1, probabilities)).to eq(0)
      expect(@engine.cdf(1, probabilities)).to eq(0)
      expect(@engine.cdf(2, probabilities)).to eq(1)
      expect(@engine.right_tail(2, probabilities)).to eq(1)
      expect(@engine.right_tail(3, probabilities)).to eq(0)
    end

    it 'uses float boundary values when float certain events are supplied' do
      probabilities = [0.0, 1.0]

      expect(@engine.pmf(1, probabilities)).to eq(1.0)
      expect(@engine.cdf(0, probabilities)).to eq(0.0)
      expect(@engine.right_tail(2, probabilities)).to eq(0.0)
    end

    it 'handles an empty probability list' do
      expect(@engine.pmf(0, [])).to eq(1)
      expect(@engine.pmf(1, [])).to eq(0)
      expect(@engine.cdf(-1, [])).to eq(0)
      expect(@engine.cdf(0, [])).to eq(1)
      expect(@engine.right_tail(0, [])).to eq(1)
      expect(@engine.right_tail(1, [])).to eq(0)
      expect(@engine.quantile(0, [])).to eq(0)
      expect(@engine.quantile(1, [])).to eq(0)
    end

    it 'returns to the binomial distribution when probabilities are equal' do
      [2, 8, 20].each do |n|
        [0.1, 0.35, 0.75].each do |probability|
          probabilities = Array.new(n, probability)

          (0..n).each do |k|
            expect(@engine.pmf(k, probabilities))
              .to be_within(1e-12).of(Distribution::Binomial::Ruby_.pdf(k, n, probability))
            expect(@engine.cdf(k, probabilities))
              .to be_within(1e-12).of(Distribution::Binomial::Ruby_.exact_cdf(k, n, probability))
          end
        end
      end
    end

    it 'calculates a tiny right tail without subtracting from one' do
      probabilities = Array.new(1000, 0.001)

      subtractive_tail = 1.0 - @engine.cdf(19, probabilities)
      expect(subtractive_tail).to be > 1e-16
      expect(@engine.right_tail(20, probabilities)).to be > 0.0
      expect(@engine.right_tail(20, probabilities)).to be < 1e-18
    end

    it 'returns integer quantiles and inverse aliases' do
      probabilities = [0.1, 0.2, 0.8]

      expect(@engine.quantile(0, probabilities)).to eq(0)
      expect(@engine.quantile(1, probabilities)).to eq(3)
      expect(@engine.p_value(@engine.cdf(1, probabilities), probabilities)).to eq(1)
      expect(@engine.exact_p_value(@engine.cdf(2, probabilities), probabilities)).to eq(2)
    end

    it 'rejects invalid input probabilities' do
      expect { @engine.pmf(0, [0.5, 1.2]) }.to raise_error(ArgumentError)
      expect { @engine.pmf(0, [0.5, -0.1]) }.to raise_error(ArgumentError)
      expect { @engine.pmf(0, [0.5, nil]) }.to raise_error(ArgumentError)
      expect { @engine.pmf(0, [Float::NAN]) }.to raise_error(ArgumentError)
      expect { @engine.quantile(-0.1, [0.5]) }.to raise_error(ArgumentError)
      expect { @engine.quantile(1.1, [0.5]) }.to raise_error(ArgumentError)
      expect { @engine.pmf(nil, [0.5]) }.to raise_error(ArgumentError)
      expect { @engine.cdf(Float::NAN, [0.5]) }.to raise_error(ArgumentError)
      expect { @engine.right_tail(Float::INFINITY, [0.5]) }.to raise_error(ArgumentError)
    end
  end

  describe 'singleton' do
    before do
      @engine = Distribution::PoissonBinomial
    end

    it_should_behave_like 'poisson-binomial engine'
  end

  describe Distribution::PoissonBinomial::Ruby_ do
    before do
      @engine = Distribution::PoissonBinomial::Ruby_
    end

    it_should_behave_like 'poisson-binomial engine'
  end

  describe Distribution::Shorthand do
    include Distribution::Shorthand

    it 'provides poisson-binomial shortcuts' do
      probabilities = [0.25, 0.75]

      expect(pbino_pmf(1, probabilities)).to eq(Distribution::PoissonBinomial.pmf(1, probabilities))
      expect(pbino_cdf(1, probabilities)).to eq(Distribution::PoissonBinomial.cdf(1, probabilities))
      expect(pbino_ccdf(1, probabilities)).to eq(Distribution::PoissonBinomial.right_tail(1, probabilities))
      expect(pbino_p(0.5, probabilities)).to eq(Distribution::PoissonBinomial.p_value(0.5, probabilities))
    end
  end
end
