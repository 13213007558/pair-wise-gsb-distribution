require 'distribution/wallenius_hypergeometric/ruby'

module Distribution
  # Wallenius' noncentral hypergeometric distribution (two colors).
  #
  # Balls of two classes are drawn one at a time, without
  # replacement, from an urn. At each draw a class is picked with
  # probability proportional to (remaining balls of the class)
  # times (weight of the class). The distribution models the number
  # of drawn balls of the interesting class and, unlike the ordinary
  # hypergeometric distribution, can express that one class is more
  # (or less) likely to be picked than the other.
  #
  # NOTE: this is *not* Fisher's noncentral hypergeometric
  # distribution. Fisher's model describes independent binomial
  # sampling conditioned on the total, while Wallenius' model
  # describes a sequential urn process. The two models differ
  # whenever the odds differ from 1.
  #
  # == Usage
  #   Distribution::WalleniusHypergeometric.pdf(2, 5, 4, 12, 3.0)
  #   Distribution::WalleniusHypergeometric.cdf(2, 5, 4, 12, 3.0)
  #   Distribution::WalleniusHypergeometric.mean(5, 4, 12, 3.0)
  #   Distribution::WalleniusHypergeometric.rng(5, 4, 12, 3.0)
  #
  # The point probabilities need a numerical integration. The
  # +*_with_status+ variants return a Result struct whose +status+
  # is :ok when the integral converged and :not_converged otherwise:
  #   r = Distribution::WalleniusHypergeometric.pdf_with_status(2, 5, 4, 12, 3.0)
  #   r.value     # probability
  #   r.status    # :ok or :not_converged
  module WalleniusHypergeometric
    SHORTHAND = 'whypg'
    extend Distributable

    create_distribution_methods

    ##
    # :singleton-method: pdf(k,m,n,total,odds)
    # Probability of drawing exactly +k+ balls of the interesting
    # class when drawing +n+ balls without replacement from an urn
    # with +m+ interesting balls out of +total+, where the
    # interesting class has relative weight +odds+ (w1/w2).

    ##
    # :singleton-method: cdf(k,m,n,total,odds)
    # Probability of drawing +k+ or fewer balls of the interesting
    # class. Same parameters as pdf.

    ##
    # :singleton-method: mean(m,n,total,odds)
    # Expected number of drawn balls of the interesting class.

    ##
    # :singleton-method: rng(m,n,total,odds,seed)
    # Sequential (ball-by-ball) sample from the distribution.
  end
end
