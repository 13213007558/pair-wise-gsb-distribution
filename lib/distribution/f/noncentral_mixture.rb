# Non-central F distributions are Poisson mixtures of central F (beta)
# distributions.  If X follows F(d1, d2; lambda), then with mu = lambda/2
# and z = d1*x/(d1*x + d2):
#
#   P(X <= x) = sum over j of w_j * I_z((d1 + 2j)/2, d2/2),
#
# where w_j = exp(-mu) mu^j / j!.  The density is the same mixture with
# the regularized betas replaced by the corresponding beta densities.
#
# Directly forming the weights overflows or underflows for large lambda,
# so every quantity here is built from log weights and the series is
# expanded outward from the Poisson mode.  The omitted part of each side
# is bounded rigorously: a geometric bound for the density terms, whose
# ratio in j is explicit and strictly decreasing, and a Poisson tail or
# Chernoff bound for probabilities.  This lets callers detect a result
# that did not reach the requested tolerance instead of silently
# accepting a truncated sum.
module Distribution
  module F
    # @api private
    module NoncentralMixture
      class << self
        # Sum term(j) * w_j over j = 0..Infinity, expanding outward from
        # the Poisson mode floor(mu).
        #
        # When up_ratio/down_ratio are nil the terms are probabilities in
        # [0, 1] (CDF and survival), so the residual is bounded by omitted
        # Poisson mass.  Otherwise they give ratios of consecutive weighted
        # density terms and a geometric bound is used.
        #
        # Returns [value, error_bound, iterations, converged].
        def sum(mu, mode, budget, tolerance,
                term, up_ratio = nil, down_ratio = nil)
          term0 = term.call(mode).to_f
          return [term0, 0.0, 1, true] if mu == 0.0

          ln_mu = Math.log(mu)
          log_w0 = log_weight(mode, mu, ln_mu)
          value = term0 > 0.0 ? Math.exp(log_w0) * term0 : 0.0

          log_w_up = log_w0
          log_w_dn = log_w0
          log_p_up = term0 > 0.0 ? log_w0 + Math.log(term0) : -Float::INFINITY
          log_p_dn = log_p_up

          up = mode
          down = mode
          iterations = 1
          bound = Float::INFINITY

          while iterations < budget
            # Expand the side with the larger current weight first, so the
            # partial sum stays closest to its limit.
            if down.zero? || log_w_up >= log_w_dn
              up += 1
              iterations += 1
              log_w_up += ln_mu - Math.log(up)
              tv = term.call(up).to_f
              if tv > 0.0
                log_p_up = log_w_up + Math.log(tv)
                value += Math.exp(log_p_up)
              end
            else
              down -= 1
              iterations += 1
              log_w_dn -= ln_mu - Math.log(down + 1)
              tv = term.call(down).to_f
              if tv > 0.0
                log_p_dn = log_w_dn + Math.log(tv)
                value += Math.exp(log_p_dn)
              else
                log_p_dn = -Float::INFINITY
              end
            end

            bound = residual_bound(
              up, down, mu, log_w_up, log_w_dn,
              log_p_up, log_p_dn, up_ratio, down_ratio
            )
            bound += rounding_error(iterations, value)

            break if bound <= tolerance
          end

          [value, bound, iterations, bound <= tolerance]
        end

        # log of the Poisson(mu) probability mass at j.
        def log_weight(j, mu, ln_mu = Math.log(mu))
          -mu + j * ln_mu - Math.lgamma(j + 1.0).first
        end

        private

        def rounding_error(iterations, value)
          6.0 * Float::EPSILON * iterations * [value.abs, Float::EPSILON].max
        end

        def residual_bound(up, down, mu, log_w_up, log_w_dn,
                           log_p_up, log_p_dn, up_ratio, down_ratio)
          if up_ratio.nil?
            probability_tail(up, down, mu, log_w_up, log_w_dn)
          else
            value_tail(up, down, log_p_up, log_p_dn, up_ratio, down_ratio)
          end
        end

        # Omitted Poisson mass above up and below down.
        def probability_tail(up, down, mu, log_w_up, log_w_dn)
          bound = 0.0

          r_up = mu / (up + 1.0)
          bound += if r_up < 1.0
                     Math.exp(log_w_up) * r_up / (1.0 - r_up)
                   else
                     chernoff(up + 1.0, mu)
                   end

          if down > 0
            r_dn = down.to_f / mu
            bound += if r_dn < 1.0
                       Math.exp(log_w_dn) * r_dn / (1.0 - r_dn)
                     else
                       chernoff(down.to_f, mu)
                     end
          end

          bound
        end

        # Geometric bound for weighted density terms beyond the included
        # range [down, up].  The ratios decrease as the series moves away
        # from its mode, so the ratio at the first omitted index dominates
        # every later one.  Bounding from the first omitted term (rather
        # than the last included one) keeps the bound rigorous.
        def value_tail(up, down, log_p_up, log_p_dn, up_ratio, down_ratio)
          # up_ratio.call(k) is p_{k+1}/p_k.  The last included index is
          # `up`; the first omitted term p_{up+1} = p_up * r(up) and every
          # later ratio is <= r(up), so the tail <= p_{up+1}/(1-r(up)).
          r_up = up_ratio.call(up).to_f
          return Float::INFINITY unless r_up > 0.0 && r_up < 1.0 &&
                                        log_p_up > -Float::INFINITY
          first_omitted_up = Math.exp(log_p_up) * r_up
          bound = first_omitted_up / (1.0 - r_up)

          if down > 0
            r_dn = down_ratio.call(down).to_f
            return Float::INFINITY unless r_dn > 0.0 && r_dn < 1.0 &&
                                          log_p_dn > -Float::INFINITY
            first_omitted_dn = Math.exp(log_p_dn) * r_dn
            bound += first_omitted_dn / (1.0 - r_dn)
          end

          bound
        end

        # Chernoff bound for a Poisson tail starting at k on the far side
        # of the mean: exp(-mu + k(1 + log(mu/k))).
        def chernoff(k, mu)
          return Float::INFINITY if k <= 0.0
          Math.exp(-mu + k * (1.0 + Math.log(mu / k)))
        end
      end
    end
  end
end
