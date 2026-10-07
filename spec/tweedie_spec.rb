require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::Tweedie do
  shared_examples_for 'a Tweedie engine' do
    let(:reference_points) do
      [
        {
          mean: 1.0, dispersion: 1.0, power: 1.5, x: 1.0,
          log_pdf: -1.02861522034198259819070795856,
          pdf: 0.3575016790048706540274189,
          log_zero_mass: -2.0,
          zero_mass: 0.1353352832366126918939995
        },
        {
          mean: 2.0, dispersion: 0.5, power: 1.3, x: 0.8,
          log_pdf: -1.3115687763144795337991729446,
          pdf: 0.269397100919060499180437,
          log_zero_mass: -4.64144226489277456649933211326,
          zero_mass: 0.009643778709101008402930025
        },
        {
          mean: 1.0, dispersion: 0.2, power: 1.8, x: 1.2,
          log_pdf: -0.385388049246295924239057462684,
          pdf: 0.6801866390733235443781777,
          log_zero_mass: -25.0000000000000041633363423443,
          zero_mass: 1.388794386496396277448035e-11
        },
        {
          mean: 0.7, dispersion: 2.0, power: 1.2, x: 0.01,
          log_pdf: -12.9088084060470500932834405912,
          pdf: 2.476143109742532715414989e-6,
          log_zero_mass: -0.469849154156278410770393015906,
          zero_mass: 0.6250965543883425997989434
        },
        {
          mean: 5.0, dispersion: 0.8, power: 1.7, x: 3.0,
          log_pdf: -2.01056509843872059877699801303,
          pdf: 0.1339129792676342489296069,
          log_zero_mass: -6.75273581955317592131004373428,
          zero_mass: 0.001167680683283349693443186
        }
      ]
    end

    it 'matches saved high-precision positive-density reference points' do
      reference_points.each do |point|
        args = [point[:x], point[:mean], point[:dispersion], point[:power]]
        expect(@engine.log_pdf(*args)).to be_within(1e-11).of(point[:log_pdf])
        expect(@engine.pdf(*args)).to be_within(1e-15).of(point[:pdf])
      end
    end

    it 'matches saved high-precision zero-mass reference points' do
      reference_points.each do |point|
        args = [point[:mean], point[:dispersion], point[:power]]
        expect(@engine.log_zero_mass(*args)).to be_within(1e-11).of(point[:log_zero_mass])
        expect(@engine.zero_mass(*args)).to be_within(1e-14).of(point[:zero_mass])
        expect(@engine.pmf(0, *args)).to eq(@engine.zero_mass(*args))
        expect(@engine.pmf(1, *args)).to eq(0.0)
      end
    end

    it 'uses a separate zero-mass interface and returns zero density off the positive support' do
      expect(@engine.pdf(-1, 1, 1, 1.5)).to eq(0.0)
      expect(@engine.pdf(0, 1, 1, 1.5)).to eq(0.0)
      expect(@engine.log_pdf(-1, 1, 1, 1.5)).to eq(-Float::INFINITY)
      expect(@engine.log_pdf(0, 1, 1, 1.5)).to eq(-Float::INFINITY)
      expect(@engine.probability_zero(1, 1, 1.5)).to eq(@engine.zero_mass(1, 1, 1.5))
    end

    it 'rejects invalid parameters and truncation tolerances' do
      expect { @engine.pdf(1, 0, 1, 1.5) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 0, 1.5) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 1, 1) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 1, 2) }.to raise_error(ArgumentError)
      expect { @engine.pdf(1, 1, 1, 1.5, tolerance: 0) }.to raise_error(ArgumentError)
    end

    it 'does not use a fixed series cutoff' do
      args = [0.01, 1.0, 1.0, 1.5]
      rough = @engine.log_pdf(*args, tolerance: 1e-6)
      accurate = @engine.log_pdf(*args, tolerance: 1e-14)

      expect(accurate).not_to eq(rough)
      expect((accurate - rough).abs).to be < 1e-6
    end

    it 'samples reproducibly with compound Poisson-Gamma moments' do
      mean = 2.0
      dispersion = 0.5
      power = 1.3
      rng_one = @engine.rng(mean, dispersion, power, 775)
      rng_two = @engine.rng(mean, dispersion, power, 775)

      expect(rng_one.call).to eq(rng_two.call)

      rng = @engine.rng(mean, dispersion, power, 775)
      n = 20_000
      samples = Array.new(n) { rng.call }

      expect(samples.min).to be >= 0.0
      sample_mean = samples.sum / n.to_f
      sample_variance = samples.sum { |value| (value - sample_mean)**2 } / (n - 1).to_f

      expect(sample_mean).to be_within(0.05).of(mean)
      expect(sample_variance).to be_within(0.08).of(dispersion * mean**power)
      expect(samples.count(0).to_f / n).to be_within(0.004).of(
        @engine.zero_mass(mean, dispersion, power)
      )
    end

    it 'samples parameters close to both power bounds' do
      [[1.0, 1.0, 1.05], [1.0, 0.2, 1.95]].each do |mean, dispersion, power|
        rng = @engine.rng(mean, dispersion, power, 31_415)
        samples = Array.new(5_000) { rng.call }
        sample_mean = samples.sum / samples.size.to_f

        expect(samples.min).to be >= 0.0
        expect(sample_mean).to be_within(0.08).of(mean)
      end
    end
  end

  describe 'singleton' do
    before { @engine = Distribution::Tweedie }
    it_should_behave_like 'a Tweedie engine'
  end

  describe Distribution::Tweedie::Ruby_ do
    before { @engine = Distribution::Tweedie::Ruby_ }
    it_should_behave_like 'a Tweedie engine'
  end
end
