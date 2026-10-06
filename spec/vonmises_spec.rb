require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::Vonmises do
  TWO_PI = 2.0 * Math::PI

  shared_examples_for 'vonmises engine' do
    # Reference values below were produced by independent trapezoid
    # quadrature of the pdf (2-4 million steps), not by the Fourier
    # series implementation under test.

    it '.pdf returns the circular uniform density for kappa = 0' do
      [0.0, 1.3, -2.7, 100.0].each do |x|
        expect(@engine.pdf(x, 0.4, 0.0)).to be_within(1e-12).of(1.0 / TWO_PI)
      end
      # Near-uniform limit: tiny concentration stays stable
      expect(@engine.pdf(1.3, 0.4, 1e-12)).to be_within(1e-9).of(1.0 / TWO_PI)
    end

    it '.pdf returns correct densities' do
      # scipy.stats.vonmises.pdf(0, 1) = 0.34171048862346314
      expect(@engine.pdf(0.0, 0.0, 1.0)).to be_within(1e-12).of(0.34171048862346314)
      expect(@engine.pdf(1.5, 0.3, 2.5)).to be_within(1e-12).of(0.119693454385016)
    end

    it '.pdf treats angles across the 0/2pi boundary as neighbours' do
      # 359 deg and 1 deg are close: densities must match around 0
      expect(@engine.pdf(359.0 * Math::PI / 180, 0.0, 3.0)).to eq(@engine.pdf(Math::PI / 180, 0.0, 3.0))
      expect(@engine.pdf(-0.01, 0.0, 3.0)).to eq(@engine.pdf(0.01, 0.0, 3.0))
    end

    it '.pdf handles a mean direction beyond one period' do
      expect(@engine.pdf(1.234, 0.5 + 6.0 * Math::PI, 2.5)).to be_within(1e-13).of(@engine.pdf(1.234, 0.5, 2.5))
    end

    it '.pdf does not overflow the Bessel constant at high concentration' do
      # Peak height ~ sqrt(kappa / (2 pi)), finite for huge kappa
      expect(@engine.pdf(0.0, 0.0, 1e6)).to be_within(1e-9).of(398.942230533625775)
      expect(@engine.pdf(0.0, 0.0, 1e12)).to be_within(1e-3).of(Math.sqrt(1e12 / TWO_PI))
      # Far from the peak the density simply underflows to 0
      expect(@engine.pdf(0.5, 0.0, 1e12)).to eq(0.0)
    end

    it '.cdf is exact in the uniform limit' do
      cut = 1.1
      [1.2, 2.5, 4.4, 7.0].each do |x|
        expect(@engine.cdf(x, 0.3, 0.0, cut)).to be_within(1e-12).of((x - cut) / TWO_PI)
      end
    end

    it '.cdf returns correct values' do
      expect(@engine.cdf(1.0, 0.0, 1.0)).to be_within(1e-9).of(0.794355307434601)
      expect(@engine.cdf(2.0, 1.0, 2.5, 0.5)).to be_within(1e-9).of(0.686183381653129)
      expect(@engine.cdf(1.9, 0.7, 1.5)).to be_within(1e-9).of(0.886871481064061)
    end

    it '.cdf goes from 0 to 1 over the caller-chosen period' do
      cut = 2.34
      expect(@engine.cdf(cut, 0.7, 2.0, cut)).to eq(0.0)
      expect(@engine.cdf(cut + TWO_PI - 1e-9, 0.7, 2.0, cut)).to be_within(1e-6).of(1.0)
      # The default cut sits opposite the mean, so the median is mu
      expect(@engine.cdf(0.7, 0.7, 2.0)).to be_within(1e-12).of(0.5)
    end

    it '.cdf is monotone increasing inside the caller-chosen period' do
      cut = 5.9
      # Half-open interval: cut + 2pi is the same angle as cut
      values = (0...20).map { |i| @engine.cdf(cut + i * TWO_PI / 20.0, 0.3, 4.0, cut) }
      expect(values.each_cons(2).all? { |a, b| b >= a }).to be(true)
    end

    it '.cdf accounts for probability on both sides of the cut' do
      # Cut at mu + 0.5 sigma of a kappa = 1e6 peak: the arc from the
      # cut forward to mu + pi carries 1 - Phi(0.5) of the mass.
      expect(@engine.cdf(Math::PI, 0.0, 1e6, 5e-4)).to be_within(1e-9).of(0.308537562563778)
      # Mass before and after any interior point sums to one full circle
      cut = 0.9
      x = 2.2
      head = @engine.cdf(x, 0.4, 3.0, cut)
      tail = @engine.arc_probability(x, cut + TWO_PI - x, 0.4, 3.0)
      expect(head + tail).to be_within(1e-12).of(1.0)
    end

    it '.quantile is exact in the uniform limit' do
      cut = 1.1
      [0.0, 0.25, 0.5, 0.75, 1.0].each do |p|
        expect(@engine.quantile(p, 0.3, 0.0, cut)).to be_within(1e-12).of(cut + TWO_PI * p)
      end
    end

    it '.quantile inverts the cdf inside the caller-chosen period' do
      [[0.7, 0.5, 1.1], [1.9, 3.0, 5.9], [0.2, 50.0, 0.0]].each do |mu, kappa, cut|
        [0.03, 0.25, 0.5, 0.77, 0.99].each do |p|
          q = @engine.quantile(p, mu, kappa, cut)
          expect(q).to be >= cut
          expect(q).to be <= cut + TWO_PI
          expect(@engine.cdf(q, mu, kappa, cut)).to be_within(1e-9).of(p)
        end
      end
    end

    it '.quantile is monotone over the caller-chosen period' do
      cut = 1.0
      qs = (1..99).map { |i| @engine.quantile(i / 100.0, 0.0, 4.0, cut) }
      expect(qs.each_cons(2).all? { |a, b| b > a }).to be(true)
      expect(@engine.quantile(0.0, 0.0, 4.0, cut)).to eq(cut)
      expect(@engine.quantile(1.0, 0.0, 4.0, cut)).to eq(cut + TWO_PI)
    end

    it '.quantile raises for probabilities outside [0, 1]' do
      expect { @engine.quantile(-0.1, 0.0, 1.0) }.to raise_error(RangeError)
      expect { @engine.quantile(1.1, 0.0, 1.0) }.to raise_error(RangeError)
    end

    it '.arc_probability is exact in the uniform limit' do
      expect(@engine.arc_probability(2.9, 1.7, 0.3, 0.0)).to be_within(1e-12).of(1.7 / TWO_PI)
      expect(@engine.arc_probability(2.9, -1.7, 0.3, 0.0)).to be_within(1e-12).of(1.7 / TWO_PI)
    end

    it '.arc_probability never confuses the full circle with an empty arc' do
      expect(@engine.arc_probability(0.7, 0.0, 0.2, 5.0)).to eq(0.0)
      expect(@engine.arc_probability(0.7, TWO_PI, 0.2, 5.0)).to eq(1.0)
      expect(@engine.arc_probability(0.7, -TWO_PI, 0.2, 5.0)).to eq(1.0)
      expect(@engine.arc_probability(0.7, 3.0 * TWO_PI, 0.2, 5.0)).to eq(1.0)
      # Just under a full circle is not the full circle
      expect(@engine.arc_probability(0.7, TWO_PI - 0.01, 0.2, 5.0)).to be < 1.0
    end

    it '.arc_probability follows short arcs across the 2pi wrap' do
      # Example: a 0.3 rad arc starting at 6.1 rad crosses 2pi
      expect(@engine.arc_probability(6.1, 0.3, 0.2, 3.0)).to be_within(1e-9).of(0.179521664709187)
      # The same arc described with angles shifted by one period
      expect(@engine.arc_probability(6.1 - TWO_PI, 0.3, 0.2, 3.0)).to be_within(1e-12).of(@engine.arc_probability(6.1, 0.3, 0.2, 3.0))
      # Example: 359 deg to 1 deg is a short 2 deg arc, not a long one
      deg = Math::PI / 180
      expect(@engine.arc_probability(359.0 * deg, 2.0 * deg, 0.0, 10.0)).to be_within(1e-9).of(0.043437310472551)
      expect(@engine.arc_probability(-1.0 * deg, 2.0 * deg, 0.0, 10.0)).to be_within(1e-12).of(@engine.arc_probability(359.0 * deg, 2.0 * deg, 0.0, 10.0))
    end

    it '.arc_probability interprets the sign of the length as direction' do
      # Clockwise arc from 0.9 of length 0.4 equals the
      # counterclockwise arc from 0.5 to 0.9
      expect(@engine.arc_probability(0.9, -0.4, 1.1, 0.6)).to be_within(1e-9).of(0.100973705314692)
      expect(@engine.arc_probability(0.9, -0.4, 1.1, 0.6)).to be_within(1e-14).of(@engine.arc_probability(0.5, 0.4, 1.1, 0.6))
    end

    it '.arc_probability of complementary arcs sums to one' do
      start = 2.6
      length = 1.9
      head = @engine.arc_probability(start, length, 0.4, 3.0)
      tail = @engine.arc_probability(start + length, TWO_PI - length, 0.4, 3.0)
      expect(head + tail).to be_within(1e-12).of(1.0)
    end

    it '.arc_probability agrees with the cdf anchored at the arc start' do
      expect(@engine.arc_probability(6.1, 0.3, 0.2, 3.0)).to be_within(1e-12).of(@engine.cdf(6.4, 0.2, 3.0, 6.1))
    end

    it 'concentrates a narrow peak like a normal distribution at high kappa' do
      # Example: kappa = 1e6 gives sigma = 1e-3; +/-3 sigma holds ~99.73%
      sigma = 1e-3
      expect(@engine.arc_probability(-3.0 * sigma, 6.0 * sigma, 0.0, 1e6)).to be_within(1e-9).of(0.997300190641155)
      expect(@engine.cdf(0.0, 0.0, 1e6)).to be_within(1e-12).of(0.5)
      expect(@engine.quantile(0.975, 0.0, 1e6)).to be_within(1e-12).of(0.001959964543249)
      expect(@engine.quantile(0.025, 0.0, 1e6)).to be_within(1e-12).of(-0.001959964543249)
    end

    it 'is invariant when every angle is shifted by whole periods' do
      [-2, 1, 5].each do |k|
        shift = k * TWO_PI
        expect(@engine.pdf(1.234 + shift, 0.5 + shift, 2.5)).to be_within(1e-13).of(@engine.pdf(1.234, 0.5, 2.5))
        expect(@engine.cdf(2.0 + shift, 1.0 + shift, 2.5, 0.5 + shift)).to be_within(1e-12).of(@engine.cdf(2.0, 1.0, 2.5, 0.5))
        expect(@engine.arc_probability(6.1 + shift, 0.3, 0.2 + shift, 3.0)).to be_within(1e-12).of(@engine.arc_probability(6.1, 0.3, 0.2, 3.0))
        expect(@engine.quantile(0.3, 0.7 + shift, 1.5, 0.1 + shift)).to be_within(1e-9).of(@engine.quantile(0.3, 0.7, 1.5, 0.1) + shift)
      end
    end

    it '.rng with a specified seed should be reproducible' do
      seed = Random.new_seed.modulo 100000007
      gen_a = @engine.rng(1.0, 5.0, seed)
      gen_b = @engine.rng(1.0, 5.0, seed)
      expect(gen_a.call).to eq(gen_b.call)
    end

    it '.rng generates angles around the mean direction' do
      gen = @engine.rng(1.0, 5.0, 12345)
      samples = Array.new(10_000) { gen.call }
      # Circular mean of the samples should be close to mu
      sin_mean = samples.sum { |a| Math.sin(a) } / samples.size
      cos_mean = samples.sum { |a| Math.cos(a) } / samples.size
      expect(Math.atan2(sin_mean, cos_mean)).to be_within(0.05).of(1.0)
    end

    it '.p_value should return same result as .quantile' do
      expect(@engine.p_value(0.42, 0.7, 2.0, 1.1)).to eq(@engine.quantile(0.42, 0.7, 2.0, 1.1))
    end
  end

  describe 'singleton' do
    before do
      @engine = Distribution::Vonmises
    end
    it_should_behave_like 'vonmises engine'
  end

  describe Distribution::Vonmises::Ruby_ do
    before do
      @engine = Distribution::Vonmises::Ruby_
    end
    it_should_behave_like 'vonmises engine'
  end
end
