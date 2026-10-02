module Distribution
  # Estimates rectangle probabilities of a finite-dimensional multivariate
  # normal distribution, i.e.
  #
  #   P(lower[i] <= X[i] <= upper[i], i = 0...m)  for  X ~ N(mean, sigma)
  #
  # using Alan Genz's separation-of-variables Monte Carlo algorithm
  # ("Numerical Computation of Multivariate Normal Probabilities",
  # J. Comput. Graph. Statist., 1992), implemented in pure Ruby.
  #
  # == Usage
  #    result = Distribution::NormalMultivariate.cdf(
  #      [0, 0], [Float::INFINITY, Float::INFINITY], # first quadrant
  #      [0, 0],                                     # zero mean
  #      [[1, 0], [0, 1]],                           # covariance
  #      random: Random.new(42)
  #    )
  #    result.estimate # => ~0.25
  #    result.error    # => statistical error estimate (see notes)
  #    result.achieved # => true if error target was met within budget
  #
  # == Notes on the error estimate
  # The returned probability is a Monte Carlo estimate, never an exact
  # value. +error+ is +alpha+ times the sample standard error of the
  # estimate, a *statistical* error estimate: it quantifies the random
  # fluctuation of the estimator itself and does not bound the true
  # integration error deterministically. Supplying a seeded random source
  # only makes a run reproducible; it does not turn the estimate into an
  # exact probability.
  module NormalMultivariate
    # Result of a NormalMultivariate.cdf call.
    #
    # * +estimate+:   Monte Carlo estimate of the rectangle probability.
    # * +error+:      statistical error estimate (+alpha+ standard errors
    #                 of the Monte Carlo mean).
    # * +achieved+:   true if the final error estimate met the requested
    #                 +epsilon+ target within the +max_iterations+ budget.
    # * +iterations+: number of Monte Carlo samples actually used.
    Result = Struct.new(:estimate, :error, :achieved, :iterations)

    class << self
      # Estimates the multivariate normal rectangle probability
      # P(lower <= X <= upper) for X ~ N(mean, sigma).
      #
      # == Arguments
      # * +lower+: array of lower bounds, one per dimension. An entry may
      #   be +-Float::INFINITY+ (or +nil+) for an unbounded side.
      # * +upper+: array of upper bounds. An entry may be
      #   ++Float::INFINITY+ (or +nil+) for an unbounded side.
      # * +mean+: array of means, one per dimension.
      # * +sigma+: covariance matrix as an array of +m+ rows of +m+
      #   numerics. Must be symmetric positive definite; an
      #   +ArgumentError+ is raised otherwise.
      # * +epsilon+: target for the statistical error estimate.
      # * +alpha+: number of standard errors used for the error estimate.
      # * +max_iterations+: Monte Carlo sample budget. Sampling always
      #   uses at least 100 iterations (or the whole budget, if smaller)
      #   so the error estimate is based on a meaningful sample size.
      # * +random+: explicit random source; must respond to +rand+
      #   returning uniform values in [0, 1). Pass a seeded
      #   +Random.new(seed)+ for reproducible runs.
      #
      # Returns a Result struct; see its documentation. The estimate is
      # statistical in nature: only in degenerate cases (e.g. a fully
      # unbounded rectangle, or a one-dimensional problem) is it exact.
      def cdf(lower, upper, mean, sigma, epsilon: 0.001, alpha: 2.5,
              max_iterations: 100_000, random: Random.new)
        m = lower.size
        unless upper.size == m && mean.size == m && sigma.size == m &&
               sigma.all? { |row| row.size == m }
          fail ArgumentError,
               'lower, upper, mean and sigma must describe the same ' +               'number of dimensions'
        end
        fail ArgumentError, 'dimension must be at least 1' if m < 1

        a = Array.new(m) do |i|
          bound = lower[i]
          bound.nil? ? -Float::INFINITY : bound.to_f - mean[i]
        end
        b = Array.new(m) do |i|
          bound = upper[i]
          bound.nil? ? Float::INFINITY : bound.to_f - mean[i]
        end
        (0...m).each do |i|
          fail ArgumentError, "lower[#{i}] exceeds upper[#{i}]" if a[i] > b[i]
        end

        c = cholesky(sigma)

        d = Array.new(m)
        e = Array.new(m)
        f = Array.new(m)
        y = Array.new(m)
        d[0] = a[0] == -Float::INFINITY ? 0.0 : uPhi(a[0] / c[0][0])
        e[0] = b[0] == Float::INFINITY ? 1.0 : uPhi(b[0] / c[0][0])
        f[0] = e[0] - d[0]

        intsum = 0.0
        varsum = 0.0
        n = 0
        error = Float::INFINITY
        # A minimum number of samples is enforced: with fewer samples the
        # variance of the estimator (and thus the error estimate) is not
        # meaningful -- with a single sample it is always zero.
        min_iterations = [100, max_iterations].min
        while n < min_iterations || (error > epsilon && n < max_iterations)
          (1...m).each do |i|
            w = random.rand
            y[i - 1] = iPhi(d[i - 1] + w * (e[i - 1] - d[i - 1]))
            sumc = 0.0
            (0...i).each { |j| sumc += c[i][j] * y[j] }
            d[i] = a[i] == -Float::INFINITY ? 0.0 : uPhi((a[i] - sumc) / c[i][i])
            e[i] = b[i] == Float::INFINITY ? 1.0 : uPhi((b[i] - sumc) / c[i][i])
            f[i] = (e[i] - d[i]) * f[i - 1]
          end
          sample = f[m - 1]
          intsum += sample
          varsum += sample * sample
          n += 1
          estimate = intsum / n
          variance = varsum / n - estimate * estimate
          # Values at the level of floating-point noise are not real
          # variance (e.g. every sample is identical in degenerate cases).
          variance = 0.0 if variance < Float::EPSILON * varsum
          error = alpha * Math.sqrt(variance / n)
        end

        Result.new(intsum / n, error, error <= epsilon, n)
      end

      # Inverse of the standard univariate normal cdf. :nodoc:
      def iPhi(pr)
        Distribution::Normal.p_value(pr)
      end

      # Standard univariate normal cdf. :nodoc:
      def uPhi(x)
        Distribution::Normal.cdf(x)
      end

      private

      # Lower-triangular Cholesky factor of a symmetric positive definite
      # matrix given as an array of rows. Raises ArgumentError unless the
      # matrix is positive definite.
      def cholesky(sigma)
        m = sigma.size
        c = Array.new(m) { Array.new(m, 0.0) }
        (0...m).each do |i|
          (0..i).each do |j|
            sum = sigma[i][j].to_f
            (0...j).each { |k| sum -= c[i][k] * c[j][k] }
            if i == j
              if sum <= 0.0
                fail ArgumentError,
                     'covariance matrix is not positive definite'
              end
              c[i][j] = Math.sqrt(sum)
            else
              c[i][j] = sum / c[j][j]
            end
          end
        end
        c
      end
    end
  end
end
