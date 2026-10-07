require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')
require 'yaml'

describe Distribution::Tweedie do
  reference_cases = YAML.load_file(
    File.expand_path('../fixtures/tweedie_reference_points.yml', __FILE__),
    permitted_classes: []
  )['cases']

  shared_examples_for 'Tweedie engine' do
    it 'provides the independent probability mass at zero' do
      reference_cases.each do |test_case|
        expected = test_case['zero_probability'].to_f
        actual = @engine.zero_probability(
          test_case['mean'], test_case['dispersion'], test_case['power']
        )

        expect(actual).to be_within(1e-15 * [expected, 1e-30].max).of(expected)
      end
    end

    it 'matches saved high-precision density and log-density reference points' do
      reference_cases.each do |test_case|
        test_case['points'].each do |point|
          arguments = [
            point['x'], test_case['mean'], test_case['dispersion'],
            test_case['power']
          ]

          expect(@engine.pdf(*arguments)).to be_within(
            5e-13 * point['density'].to_f.abs
          ).of(point['density'].to_f)
          expect(@engine.log_pdf(*arguments)).to be_within(2e-12).of(
            point['log_density'].to_f
          )
        end
      end
    end

    it 'stops the count series from a controlled relative truncation bound' do
      large_count_case = reference_cases.last
      point = large_count_case['points'].first
      arguments = [
        point['x'], large_count_case['mean'], large_count_case['dispersion'],
        large_count_case['power']
      ]
      expected = point['density'].to_f

      expect(point['terms_needed']).to be > 100
      expect(@engine.pdf(*arguments, 1e-8)).to be_within(5e-8 * expected).of(expected)
      expect(@engine.pdf(*arguments, 1e-12)).to be_within(5e-13 * expected).of(expected)
      expect(@engine.log_pdf(*arguments, 1e-14)).to be_within(1e-12).of(
        point['log_density'].to_f
      )
    end

    it 'returns zero density for non-positive x and keeps zero mass separate' do
      expect(@engine.pdf(0, 2, 0.5, 1.5)).to eq(0.0)
      expect(@engine.pdf(-1, 2, 0.5, 1.5)).to eq(0.0)
      expect(@engine.log_pdf(0, 2, 0.5, 1.5)).to eq(-Float::INFINITY)
      expect(@engine.log_pdf(-1, 2, 0.5, 1.5)).to eq(-Float::INFINITY)
      expect(@engine.zero_probability(2, 0.5, 1.5)).to be > 0
    end

    it 'validates the mean, dispersion, variance power, and tolerance' do
      expect { @engine.pdf(1, 0, 1, 1.5) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 0, 1.5) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 1, 1.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 1, 2.0) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 1, 1.5, 0) }.to raise_error(ArgumentError)
    end

    it 'samples the compound Poisson-Gamma mean, variance, and zero mass' do
      mean = 2.3
      dispersion = 0.7
      power = 1.5
      sampler = @engine.rng(mean, dispersion, power, 10_472_932)
      samples = Array.new(30_000) { sampler.call }
      sample_mean = samples.sum / samples.size.to_f
      sample_variance = samples.sum { |value| (value - sample_mean)**2 } /
                        (samples.size - 1).to_f
      zero_fraction = samples.count(&:zero?).fdiv(samples.size)

      expect(samples).to all(be >= 0)
      expect(sample_mean).to be_within(0.06).of(mean)
      expect(sample_variance).to be_within(0.18).of(
        dispersion * mean**power
      )
      expect(zero_fraction).to be_within(0.006).of(
        @engine.zero_probability(mean, dispersion, power)
      )
    end

    it 'produces reproducible samples from a seed or injected Random object' do
      first = @engine.rng(2, 0.5, 1.4, 987)
      second = @engine.rng(2, 0.5, 1.4, 987)
      values_a = Array.new(8) { first.call }
      values_b = Array.new(8) { second.call }
      injected = @engine.rng(2, 0.5, 1.4, random: Random.new(987))
      values_c = Array.new(8) { injected.call }

      expect(values_a).to eq(values_b)
      expect(values_a).to eq(values_c)
    end
  end

  describe 'singleton' do
    before { @engine = Distribution::Tweedie }
    it_should_behave_like 'Tweedie engine'
  end

  describe Distribution::Tweedie::Ruby_ do
    before { @engine = Distribution::Tweedie::Ruby_ }
    it_should_behave_like 'Tweedie engine'
  end
end
