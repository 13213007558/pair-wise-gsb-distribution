require 'distribution/stable/ruby'

module Distribution
  # Nolan's S_0 parameterization is the public four-parameter interface.
  # Use Stable.s1_parameters and Stable.from_s1_parameters when working
  # with S_1 inputs.
  module Stable
    SHORTHAND = 'stab'
    extend Distributable
    create_distribution_methods

    class IntegrationError < StandardError
      attr_reader :diagnostics

      def initialize(message, diagnostics)
        super(message)
        @diagnostics = diagnostics
      end
    end
  end
end
