require File.expand_path(File.dirname(__FILE__) + '/spec_helper.rb')

describe Distribution::PhaseType do
  # Hand-computed reference: Erlang-2 with rate 2 (two exponential
  # phases in series). S(t) = exp(-2t)(1 + 2t), f(t) = 4t exp(-2t),
  # E[T^k] = (k + 1)! / 2^k.
  describe 'series Erlang-2 model' do
    before do
      @erlang = Distribution::PhaseType.new([1, 0], [[-2, 2], [0, -2]])
    end

    it 'should have no mass at zero' do
      expect(@erlang.zero_mass).to be_within(1e-12).of(0.0)
    end

    it 'should return the hand-computed survival function' do
      [0.1, 0.5, 1.0, 2.3].each do |t|
        expect(@erlang.survival(t)).to be_within(1e-12).of(Math.exp(-2 * t) * (1 + 2 * t))
      end
      expect(@erlang.survival(0)).to eq 1.0
    end

    it 'should return the hand-computed density' do
      [0.1, 0.5, 1.0, 2.3].each do |t|
        expect(@erlang.pdf(t)).to be_within(1e-12).of(4 * t * Math.exp(-2 * t))
      end
    end

    it 'should return the hand-computed cdf' do
      [0.1, 0.5, 1.0, 2.3].each do |t|
        expect(@erlang.cdf(t)).to be_within(1e-12).of(1 - Math.exp(-2 * t) * (1 + 2 * t))
      end
    end

    it 'should return the hand-computed moments' do
      expect(@erlang.moment(1)).to be_within(1e-12).of(1.0)
      expect(@erlang.moment(2)).to be_within(1e-12).of(1.5)
      expect(@erlang.moment(3)).to be_within(1e-12).of(3.0)
      expect(@erlang.mean).to be_within(1e-12).of(1.0)
      expect(@erlang.variance).to be_within(1e-12).of(0.5)
    end

    it 'should return the hand-computed residual lifetime distribution' do
      # At t = 0.5 the state probabilities are proportional to
      # [exp(-1), exp(-1)], so the residual starts in each phase
      # with probability 1/2 and S_residual(s) = exp(-2s)(1 + s).
      residual = @erlang.residual_lifetime(0.5)
      expect(residual.alpha[0]).to be_within(1e-12).of(0.5)
      expect(residual.alpha[1]).to be_within(1e-12).of(0.5)
      expect(residual.zero_mass).to be_within(1e-12).of(0.0)
      [0.1, 0.3, 1.2].each do |s|
        expect(residual.survival(s)).to be_within(1e-12).of(Math.exp(-2 * s) * (1 + s))
      end
      # Residual mean: integral of its survival function = 3/4.
      expect(residual.mean).to be_within(1e-12).of(0.75)
    end
  end

  # Hand-computed reference: hyperexponential, branching at time zero.
  # alpha = [0.3, 0.7], T = diag(-1, -4).
  # S(t) = 0.3 exp(-t) + 0.7 exp(-4t), f(t) = 0.3 exp(-t) + 2.8 exp(-4t).
  describe 'two-phase model branching at initialization' do
    before do
      @hyper = Distribution::PhaseType.new([0.3, 0.7], [[-1, 0], [0, -4]])
    end

    it 'should return the hand-computed survival function and density' do
      [0.2, 0.7, 1.5].each do |t|
        expect(@hyper.survival(t)).to be_within(1e-12).of(0.3 * Math.exp(-t) + 0.7 * Math.exp(-4 * t))
        expect(@hyper.pdf(t)).to be_within(1e-12).of(0.3 * Math.exp(-t) + 2.8 * Math.exp(-4 * t))
      end
    end

    it 'should return the hand-computed moments' do
      expect(@hyper.mean).to be_within(1e-12).of(0.3 + 0.7 / 4.0)
      expect(@hyper.moment(2)).to be_within(1e-12).of(2 * 0.3 + 2 * 0.7 / 16.0)
      expect(@hyper.variance).to be_within(1e-12).of(0.6875 - 0.475**2)
    end

    it 'should reweight the branches in the residual lifetime' do
      t = 0.7
      residual = @hyper.residual_lifetime(t)
      total = 0.3 * Math.exp(-t) + 0.7 * Math.exp(-4 * t)
      expect(residual.alpha[0]).to be_within(1e-12).of(0.3 * Math.exp(-t) / total)
      expect(residual.alpha[1]).to be_within(1e-12).of(0.7 * Math.exp(-4 * t) / total)
    end
  end

  # Hand-computed reference: Coxian model, branching after the first
  # phase. alpha = [1, 0], T = [[-3, 2], [0, -2]].
  # S(t) = 2 exp(-2t) - exp(-3t), f(t) = 4 exp(-2t) - 3 exp(-3t),
  # E[T] = 2/3, E[T^2] = 7/9, Var(T) = 1/3.
  describe 'two-phase model branching after the first phase' do
    before do
      @coxian = Distribution::PhaseType.new([1, 0], [[-3, 2], [0, -2]])
    end

    it 'should return the hand-computed survival function and density' do
      [0.2, 0.7, 1.5].each do |t|
        expect(@coxian.survival(t)).to be_within(1e-12).of(2 * Math.exp(-2 * t) - Math.exp(-3 * t))
        expect(@coxian.pdf(t)).to be_within(1e-12).of(4 * Math.exp(-2 * t) - 3 * Math.exp(-3 * t))
      end
    end

    it 'should return the hand-computed moments' do
      expect(@coxian.mean).to be_within(1e-12).of(2.0 / 3)
      expect(@coxian.moment(2)).to be_within(1e-12).of(7.0 / 9)
      expect(@coxian.variance).to be_within(1e-12).of(1.0 / 3)
    end
  end

  describe 'mass at zero' do
    before do
      @defective = Distribution::PhaseType.new([0.6, 0.2], [[-2, 2], [0, -2]])
    end

    it 'should be reported separately instead of being normalized away' do
      expect(@defective.zero_mass).to be_within(1e-12).of(0.2)
      expect(@defective.survival(0)).to be_within(1e-12).of(0.8)
      expect(@defective.cdf(0)).to be_within(1e-12).of(0.2)
    end

    it 'should keep the continuous part unnormalized' do
      t = 0.4
      expect(@defective.survival(t)).to be_within(1e-12).of(0.6 * Math.exp(-2 * t) * (1 + 2 * t) + 0.2 * Math.exp(-2 * t))
      expect(@defective.mean).to be_within(1e-12).of(0.7)
    end
  end

  describe 'model validation' do
    it 'should refuse a generator that never absorbs' do
      # Closed communicating class: rows sum to zero, generator singular.
      expect do
        Distribution::PhaseType.new([1, 0], [[-1, 1], [1, -1]])
      end.to raise_error(Distribution::PhaseType::NotAbsorbingError)
    end

    it 'should refuse a single-state model with zero exit rate' do
      expect do
        Distribution::PhaseType.new([1], [[0]])
      end.to raise_error(Distribution::PhaseType::InvalidModelError)
    end

    it 'should refuse a non-substochastic generator' do
      expect do
        Distribution::PhaseType.new([1, 0], [[-1, 2], [0, -1]])
      end.to raise_error(Distribution::PhaseType::InvalidModelError)
    end

    it 'should refuse negative off-diagonal rates' do
      expect do
        Distribution::PhaseType.new([1, 0], [[-1, -0.5], [0, -1]])
      end.to raise_error(Distribution::PhaseType::InvalidModelError)
    end

    it 'should refuse a non-square generator' do
      expect do
        Distribution::PhaseType.new([1, 0], [[-1, 1, 0], [0, -1, 1]])
      end.to raise_error(Distribution::PhaseType::InvalidModelError)
    end

    it 'should refuse an initial vector summing above one' do
      expect do
        Distribution::PhaseType.new([0.7, 0.7], [[-2, 2], [0, -2]])
      end.to raise_error(Distribution::PhaseType::InvalidModelError)
    end

    it 'should refuse a mismatched initial vector' do
      expect do
        Distribution::PhaseType.new([1, 0, 0], [[-2, 2], [0, -2]])
      end.to raise_error(Distribution::PhaseType::InvalidModelError)
    end
  end
end
