require 'distribution/wallenius_hypergeometric/ruby'

module Distribution
  # Two-color Wallenius noncentral hypergeometric distribution.
  #
  # This model describes sequential draws without replacement when each of
  # two colors may have a different draw weight.  Remaining items are sampled
  # in proportion to their current weights, so removing an item changes the
  # next draw probability.
  #
  # This is the Wallenius, or sequentially-dependent, noncentral hypergeometric
  # distribution.  It is not Fisher's noncentral hypergeometric model.  Fisher's
  # model is an independent/parallel weighting model and has a different mass
  # function.
  #
  # == Parameters
  # +k+:: number of drawn items from the first color
  # +m1+:: initial number of items of the first color
  # +m2+:: initial number of items of the second color
  # +n+:: number of items drawn without replacement
  # +odds+:: per-item draw-weight ratio of the first color to the second
  #
  # == Example
  #   Distribution::WalleniusNoncentralHypergeometric.pdf(2, 4, 3, 5, 2.0)
  module WalleniusNoncentralHypergeometric
    SHORTHAND = 'wnhypg'
    extend Distributable

    create_distribution_methods

    alias_method :pmf, :pdf if singleton_methods.include? :pdf
  end

  WalleniusHypergeometric = WalleniusNoncentralHypergeometric
end
