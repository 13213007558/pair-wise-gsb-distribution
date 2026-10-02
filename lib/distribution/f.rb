require 'distribution/f/ruby'
require 'distribution/f/gsl'
require 'distribution/f/statistics2'
require 'distribution/f/java'
module Distribution
  # Calculate cdf and inverse cdf for F Distribution.
  #
  module F
    SHORTHAND = 'fdist'
    extend Distributable
    create_distribution_methods

    class << self
      def pdf(x, k1, k2, ncp = nil, **options)
        return Distribution::NoncentralF::Ruby_.pdf(x, k1, k2, ncp, **options) unless ncp.nil?
        engine_for(:pdf).pdf(x, k1, k2)
      end

      def cdf(x, k1, k2, ncp = nil, **options)
        return Distribution::NoncentralF::Ruby_.cdf(x, k1, k2, ncp, **options) unless ncp.nil?
        engine_for(:cdf).cdf(x, k1, k2)
      end

      def survival(x, k1, k2, ncp = nil, **options)
        return Distribution::NoncentralF::Ruby_.survival(x, k1, k2, ncp, **options) unless ncp.nil?
        Distribution::F::Ruby_.q(x, k1, k2)
      end

      alias_method :q, :survival

      def quantile(probability, k1, k2, ncp = nil, **options)
        return Distribution::NoncentralF::Ruby_.quantile(probability, k1, k2, ncp, **options) unless ncp.nil?
        engine_for(:quantile).quantile(probability, k1, k2)
      end

      alias_method :p_value, :quantile
      alias_method :inverse_cdf, :quantile

      private

      def engine_for(method)
        engine_name = Distribution.libraries_order.find do |name|
          const_get(name).respond_to?(method)
        end
        const_get(engine_name || :Ruby_)
      end
    end

    ##
    # :singleton-method: pdf(x,k1,k2)
    # Returns the PDF of F distribution
    # with +k1+ and +k2+ degrees of freedom over [0, +x+]

    ##
    # :singleton-method: p_value(qn, k1, k2)
    # Return the P-value of the corresponding integral +qn+ with
    # +k1+ and +k2+ degrees of freedom

    ##
    # :singleton-method: cdf(x,k1,k2)
    # Returns the integral of F distribution
    # with +k1+ and +k2+ degrees of freedom over [0, +x+]
  end
end
