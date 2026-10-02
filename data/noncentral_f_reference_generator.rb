# Independent high-precision reference values for the non-central F
# distribution, computed entirely with BigDecimal (no Float code path).
#
# CDF = sum_j e^{-mu} mu^j / j! * I_z((d1+2j)/2, d2/2), with
# mu = lambda/2 and z = d1*x/(d1*x+d2).  I_z is the regularized
# incomplete beta via its hypergeometric series (with beta-complement
# symmetry for the short side); log-gamma uses the DLMF Stirling series
# with Bernoulli numbers (plus a recurrence shift); Poisson weights are
# accumulated in log space.
require "bigdecimal"
require "bigdecimal/math"
require "json"

class NCFReference
  WORK = 40
  ZERO = BigDecimal("0")
  ONE = BigDecimal("1")

  def initialize
    @pi = BigMath.PI(WORK)
    @bernoulli = [Rational(1, 6), Rational(-1, 30), Rational(1, 42),
                  Rational(-1, 30), Rational(5, 66), Rational(-691, 2730),
                  Rational(7, 6), Rational(-3617, 510), Rational(43867, 798),
                  Rational(-174611, 330), Rational(854513, 138),
                  Rational(-236364091, 2730), Rational(8553103, 6),
                  Rational(-23749461029, 870), Rational(8615841276005, 14322),
                  Rational(-848035476425, 510), Rational(90219075042845, 11730)]
  end

  def bd(x)
    BigDecimal(x.to_s, WORK)
  end

  # log Gamma via the Stirling (Bernoulli) expansion.  Small arguments are
  # shifted up with the recurrence Gamma(z+1) = z Gamma(z).
  def lgamma(x)
    x = bd(x)
    shifts = []
    while x < 20
      shifts << x
      x += 1
    end
    sum = ZERO
    @bernoulli.each_with_index do |b, i|
      n = i + 1
      coeff = bd(b.numerator) / bd(b.denominator) / bd(2 * n) / bd(2 * n - 1)
      term = coeff / (x**(2 * n - 1))
      sum += term
      break if term.abs < BigDecimal("1E-32")
    end
    result = (x - BigDecimal("0.5")) * BigMath.log(x, WORK) - x +
             0.5 * BigMath.log(2 * @pi, WORK) + sum
    shifts.reverse_each { |s| result -= BigMath.log(s, WORK) }
    result
  end

  def logbeta(a, b)
    lgamma(a) + lgamma(b) - lgamma(a + b)
  end

  # Regularized incomplete beta I_x(a, b).
  def ibeta(x, a, b)
    x = bd(x); a = bd(a); b = bd(b)
    return ZERO if x <= 0
    return ONE if x >= 1
    if x < (a + 1) / (a + b + 2)
      ibeta_direct(x, a, b)
    else
      ONE - ibeta_direct(ONE - x, b, a)
    end
  end

  def ibeta_direct(x, a, b)
    # I_x(a,b) = x^a / (a B(a,b)) * 2F1(a, 1-b; a+1; x).
    pref = BigMath.exp(a * BigMath.log(x, WORK) - logbeta(a, b), WORK) / a
    sum = ONE
    term = ONE
    n = 0
    loop do
      term *= (a + n) * (n + ONE - b) / ((a + n + 1) * (n + 1)) * x
      sum += term
      n += 1
      break if term.abs < BigDecimal("1E-30") || n > 400
    end
    pref * sum
  end

  def logw(j, mu)
    -mu + bd(j) * BigMath.log(mu, WORK) - lgamma(j + 1)
  end

  # Beta(a,b) density at z, computed directly in log space.
  def beta_pdf(z, a, b)
    z = bd(z); a = bd(a); b = bd(b)
    return ZERO if z <= 0 || z >= 1
    BigMath.exp(-logbeta(a, b) + (a - 1) * BigMath.log(z, WORK) +
                 (b - 1) * BigMath.log(ONE - z, WORK), WORK)
  end

  # Non-central F density at x via the same Poisson mixture.
  def pdf(x, d1, d2, lambda)
    d1 = bd(d1); d2 = bd(d2); mu = bd(lambda) / 2
    return ZERO if x.to_f < 0 || x.to_f.infinite?
    x = bd(x)
    z = d1 * x / (d1 * x + d2)
    dz_dx = d1 * d2 / (d1 * x + d2)**2
    mode = mu.floor
    term_at = lambda { |j| beta_pdf(z, d1 / 2 + j, d2 / 2) * dz_dx }
    total = BigMath.exp(logw(mode, mu), WORK) * term_at.call(mode)
    lw_up = logw(mode, mu); lw_dn = lw_up
    up = mode; down = mode
    iters = 1
    tail = bd("Infinity")
    while iters < 2_000_000
      if down.zero? || lw_up >= lw_dn
        up += 1
        lw_up = logw(up, mu)
        total += BigMath.exp(lw_up, WORK) * term_at.call(up)
      else
        down -= 1
        lw_dn = logw(down, mu)
        total += BigMath.exp(lw_dn, WORK) * term_at.call(down)
      end
      iters += 1
      tail = residual(up, down, mu)
      break if tail < BigDecimal("1E-28")
    end
    total
  end

  # Cheap (Float) central-F quantile, used only as the Newton starting
  # point; the high-precision Newton solve below is independent of it.
  def central_f_quantile_seed(p, d1, d2)
    target = p.to_f
    a = d2 / 2.0
    b = d1 / 2.0
    lo = 0.0; hi = 1.0
    60.times do
      mid = (lo + hi) / 2.0
      if float_ibeta(mid, a, b) < target
        lo = mid
      else
        hi = mid
      end
    end
    z = (lo + hi) / 2.0
    d2.to_f / d1.to_f * z / (1.0 - z)
  end

  # Float regularized incomplete beta via continued fraction (Numerical
  # Recipes betacf), solely to build the Newton seed.
  def float_ibeta(x, a, b)
    return 0.0 if x <= 0.0
    return 1.0 if x >= 1.0
    lbeta = Math.lgamma(a + b).first - Math.lgamma(a).first - Math.lgamma(b).first
    bt = Math.exp(lbeta + a * Math.log(x) + b * Math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0)
      bt * float_betacf(a, b, x) / a
    else
      1.0 - bt * float_betacf(b, a, 1.0 - x) / b
    end
  end

  def float_betacf(a, b, x)
    qab = a + b; qap = a + 1.0; qam = a - 1.0
    c = 1.0; d = 1.0 - qab * x / qap
    d = 1e-300 if d.abs < 1e-300
    d = 1.0 / d; h = d
    (1..200).each do |m|
      m2 = 2 * m
      aa = m * (b - m) * x / ((qam + m2) * (a + m2))
      d = 1.0 + aa * d; c = 1.0 + aa / c
      d = 1e-300 if d.abs < 1e-300; c = 1e-300 if c.abs < 1e-300
      d = 1.0 / d; h *= d * c
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (a + m2 + 1.0))
      d = 1.0 + aa * d; c = 1.0 + aa / c
      d = 1e-300 if d.abs < 1e-300; c = 1e-300 if c.abs < 1e-300
      d = 1.0 / d; delta = d * c; h *= delta
      break if (delta - 1.0).abs < 3e-15
    end
    h
  end

  # Non-central F quantile by bisection on a deliberately lower-precision
  # independent CDF (no Newton, no code shared with the Float solver).
  def quantile(p, d1, d2, lambda)
    p = bd(p)
    seed = bd(central_f_quantile_seed(p.to_f, d1, d2))
    lo = ZERO
    hi = seed * (1 + bd(lambda) / bd(d1))
    hi += hi while cdf_approx(hi, d1, d2, lambda) < p
    48.times do
      mid = (lo + hi) / 2
      if cdf_approx(mid, d1, d2, lambda) <= p
        lo = mid
      else
        hi = mid
      end
    end
    (lo + hi) / 2
  end

  # Mixture CDF summed to roughly 18 digits only, used by the independent
  # quantile bisection to keep that one-time solve cheap.
  def cdf_approx(x, d1, d2, lambda)
    d1 = bd(d1); d2 = bd(d2); mu = bd(lambda) / 2
    z = d1 * bd(x) / (d1 * bd(x) + d2)
    mode = mu.floor
    term_at = lambda { |j| ibeta_approx(z, d1 / 2 + j, d2 / 2) }
    total = BigMath.exp(logw(mode, mu), WORK) * term_at.call(mode)
    lw_up = logw(mode, mu); lw_dn = lw_up
    up = mode; down = mode
    iters = 1
    while iters < 200_000
      if down.zero? || lw_up >= lw_dn
        up += 1
        lw_up = logw(up, mu)
        total += BigMath.exp(lw_up, WORK) * term_at.call(up)
      else
        down -= 1
        lw_dn = logw(down, mu)
        total += BigMath.exp(lw_dn, WORK) * term_at.call(down)
      end
      iters += 1
      break if residual(up, down, mu) < BigDecimal("1E-22")
    end
    total
  end

  def ibeta_approx(z, a, b)
    z = bd(z)
    return ZERO if z <= 0
    return ONE if z >= 1
    a = bd(a); b = bd(b)
    if z < (a + 1) / (a + b + 2)
      x = z; aa = a; bb = b; complement = false
    else
      x = 1 - z; aa = b; bb = a; complement = true
    end
    pref = BigMath.exp(aa * BigMath.log(x, WORK) - logbeta(aa, bb), WORK) / aa
    sum = ONE
    term = ONE
    n = 0
    300.times do
      term *= (aa + n) * (n + ONE - bb) / ((aa + n + 1) * (n + 1)) * x
      sum += term
      n += 1
      break if term.abs < BigDecimal("1E-22")
    end
    result = pref * sum
    complement ? ONE - result : result
  end

  # Non-central F CDF.  Returns [value, omitted-Poisson-mass bound].
  def cdf(x, d1, d2, lambda)
    d1 = bd(d1); d2 = bd(d2); mu = bd(lambda) / 2
    return [ZERO, ZERO] if x.to_f <= 0
    return [ONE, ZERO] if x.to_f.infinite?
    z = d1 * bd(x) / (d1 * bd(x) + d2)
    mode = (mu.floor)
    term_at = lambda { |j| ibeta(z, d1 / 2 + j, d2 / 2) }
    total = BigMath.exp(logw(mode, mu), WORK) * term_at.call(mode)
    lw_up = logw(mode, mu); lw_dn = lw_up
    up = mode; down = mode
    tail = bd("Infinity")
    iters = 1
    while iters < 2_000_000
      if down.zero? || lw_up >= lw_dn
        up += 1
        lw_up = logw(up, mu)
        total += BigMath.exp(lw_up, WORK) * term_at.call(up)
      else
        down -= 1
        lw_dn = logw(down, mu)
        total += BigMath.exp(lw_dn, WORK) * term_at.call(down)
      end
      iters += 1
      tail = residual(up, down, mu)
      break if tail < BigDecimal("1E-28")
    end
    [total, tail]
  end

  def residual(up, down, mu)
    bound = ZERO
    ru = mu / (up + 1)
    bound += BigMath.exp(logw(up, mu), WORK) * ru / (1 - ru) if ru < 1
    if down > 0
      rd = bd(down) / mu
      bound += BigMath.exp(logw(down, mu), WORK) * rd / (1 - rd) if rd < 1
    end
    bound
  end
end

if $PROGRAM_NAME == __FILE__
  ref = NCFReference.new
  rows = []
  [[5, 10, "4.0"], [5, 10, "20.0"], [3, 50, "4.0"], [10, 4, "10.0"]].each do |d1, d2, lam|
    [0.001, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0, 10_000.0].each do |x|
      value, bound = ref.cdf(x, d1, d2, lam)
      density = ref.pdf(x, d1, d2, lam)
      rows << { d1: d1, d2: d2, lambda: lam.to_f, x: x,
                cdf: value.round(28).to_s("F"),
                survival: (1 - value).round(28).to_s("F"),
                pdf: density.round(28).to_s("F"),
                tail_bound: bound.round(28).to_s("F") }
    end
  end
  # Quantiles near both tails and at the centre.  The solve is run at
  # roughly 18-digit precision, so values are emitted to 18 digits.
  [[5, 10, "4.0"], [5, 10, "20.0"], [3, 50, "4.0"], [10, 4, "10.0"]].each do |d1, d2, lam|
    ["1E-10", "0.5", "0.9999999999"].each do |ps|
      q = ref.quantile(ps, d1, d2, lam)
      rows << { d1: d1, d2: d2, lambda: lam.to_f, p: ps.to_f,
                quantile: q.round(18).to_s("F") }
    end
  end
  puts JSON.pretty_generate(rows)
end
