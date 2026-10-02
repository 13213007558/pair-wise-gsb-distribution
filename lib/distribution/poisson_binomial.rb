require 'distribution/poisson_binomial/ruby'

module Distribution
  # Calculate statisticals for the Poisson-binomial distribution:
  # the distribution of the number of successes in +n+ independent
  # Bernoulli trials with individual success probabilities
  # (e.g. independent sensors, each with its own failure probability).
  #
  # All methods take an Array of probabilities as parameter:
  #   Distribution::PoissonBinomial.sf(3, [0.1, 0.2, 0.05, 0.5])
  #   => P(at least 3 failures)
  #
  # Probabilities of 0 and 1 are treated as deterministic events.
  # If all probabilities are equal to +p+, results coincide with
  # Distribution::Binomial with parameters (n, p).
  #
  # == Numerical notes and input limits
  #
  # * The pmf is computed by O(n^2) convolution and renormalized so
  #   the total probability mass is exactly 1.
  # * Accumulated roundoff grows slowly with n (roughly n * 2^-53);
  #   practical input sizes go up to n ~ 10^4.
  # * +sf+ (right tail) sums the pmf from the top down, so tiny tail
  #   probabilities remain accurate where 1 - cdf would underflow.
  # * Individual masses below ~1e-308 underflow to 0.0.
  module PoissonBinomial
    SHORTHAND = 'poibino'

    extend Distributable
    create_distribution_methods

    ##
    # :singleton-method: pmf_array(probs)
    # Returns the full pmf as an Array of Floats:
    # entry +i+ is P(X = i). The array sums to 1.

    ##
    # :singleton-method: pdf(k, probs)
    # Returns P(X = k) for independent Bernoulli trials
    # with probabilities +probs+

    ##
    # :singleton-method: cdf(k, probs)
    # Returns P(X <= k) for independent Bernoulli trials
    # with probabilities +probs+

    ##
    # :singleton-method: sf(k, probs)
    # Returns the right tail P(X >= k), computed without
    # 1 - cdf cancellation, for independent Bernoulli trials
    # with probabilities +probs+

    ##
    # :singleton-method: p_value(qn, probs)
    # Returns the quantile: smallest k with P(X <= k) >= +qn+
  end
end
