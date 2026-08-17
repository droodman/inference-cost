# Stochastic frontier model with a one-sided truncated-normal inefficiency term u,
# combined with the Papke-Wooldridge fractional logit response for output in [0, 1].
#
#   ln y*_i = x_i'beta - u_i,      u_i ~ N+(mu_i, sigma_u^2)  (normal truncated below at 0)
#   mu_i    = w_i'delta
#   l_i(u)  = y_i * z_i - log(1 + exp(z_i)),   z_i = ln y*_i = x_i'beta - u
#
# (This is the y_i*ln(y*_i) - log(1+y*_i) form after substituting y*_i = exp(z_i);
#  the z-space version above is the numerically stable one to compute with.)
#
# u is latent, so the observation likelihood integrates it out:
#   L_i = int_0^Inf exp{l_i(u)} f_N+(u; mu_i, sigma_u^2) du
# The logistic link breaks the usual Aigner-Lovell-Schmidt closed-form convolution,
# so this is evaluated by quadrature. Substituting p = F_u(u) (the truncated-normal
# CDF) turns the integral into int_0^1 exp{l_i(u(p))} dp, which a FIXED Gauss-Legendre
# rule on (0,1) can integrate with the same nodes for every observation -- this keeps
# the likelihood surface smooth as parameters (and hence mu_i, sigma_u) move during
# optimization, unlike per-observation adaptive quadrature.

suppressMessages(library(maxLik))

## ---- stable elementary pieces ----------------------------------------------

log1pexp <- function(z) {
  # numerically stable log(1 + exp(z)) for any sign/magnitude of z
  out <- z
  small <- z <= 0
  out[small]  <- log1p(exp(z[small]))
  out[!small] <- z[!small] + log1p(exp(-z[!small]))
  out
}

gauss_legendre_01 <- function(n) {
  # n-point Gauss-Legendre nodes/weights on (0,1)
  gl <- pracma::gaussLegendre(n, 0, 1)
  list(nodes = gl$x, weights = gl$w)
}

## ---- likelihood -------------------------------------------------------------

make_frontier_loglik <- function(y, Xb, Xmu, quad) {
  K <- length(quad$nodes)
  logw <- log(quad$weights)
  log1m_p <- log1p(-quad$nodes)   # log(1 - p_k), node-dependent only
  skeleton <- list(beta = numeric(ncol(Xb)), delta = numeric(ncol(Xmu)), log_sigma_u = 0)

  function(theta) {
    p <- relist(theta, skeleton)
    beta    <- p$beta
    delta   <- p$delta
    sigma_u <- exp(p$log_sigma_u)

    eta <- as.vector(Xb %*% beta)     # x_i'beta
    mu  <- as.vector(Xmu %*% delta)   # pre-truncation mean of u_i

    # Quantile-transform quadrature nodes p_k into u-space via
    #   u(p) = mu + sigma_u * qnorm(Phi_c + (1-Phi_c)*p),   Phi_c = pnorm(-mu/sigma_u)
    # A large fraction of observations here have mu/sigma_u far from 0 (mu_i very
    # negative relative to sigma_u => Phi_c essentially 1: u is squeezed hard against
    # its lower truncation bound), so evaluating this directly loses precision or
    # returns qnorm(1) = Inf. Route it through the exact identity q(x) = -q(1-x) and
    # 1-x = (1-Phi_c)*(1-p) = Phi(mu/sigma_u)*(1-p), which R's log.p=TRUE machinery
    # keeps accurate arbitrarily far into either tail -- no clipping required.
    log_1mPhi_c <- pnorm(mu, sd = sigma_u, log.p = TRUE)           # log(1 - Phi_c)
    log_1mx <- outer(log_1mPhi_c, log1m_p, "+")                    # n x K: log(1-x_{i,k})
    Q <- -qnorm(log_1mx, log.p = TRUE)                             # n x K
    U <- mu + sigma_u * Q                              # n x K
    Z <- eta - U                                       # n x K: z_i(u) at each node
    H <- y * Z - log1pexp(Z)                           # n x K: l_i(u) at each node

    M <- H + matrix(logw, nrow(H), K, byrow = TRUE)
    ll <- matrixStats::rowLogSumExps(M)   # per-observation log L_i (vector -> BHHH-ready)

    # Analytic score, returned as the "gradient" attribute maxLik looks for.
    # This differentiates the QUADRATURE SUM computed above -- nodes p_k and
    # weights w_k are constants -- so gradient and value are exactly consistent
    # (differentiating the underlying integral instead would give a slightly
    # different object at finite K, desynchronizing the optimizer's line search
    # and convergence test from the objective it is actually climbing).
    #
    #   dl_i/dtheta = sum_k q_ik * dH_ik/dtheta,  q_ik = softmax_k(M_ik)
    #   dH/dZ = y_i - G(Z),   dZ = -dU,   U = mu + sigma_u*Q
    #   dU/dmu = 1 - R,   dU/dlog_sigma_u = sigma_u*Q + mu*R = U + mu*(R - 1)
    # where R = (1-p_k)*phi(mu/sigma_u)/phi(Q) is a ratio of normal densities.
    # Both phi() terms underflow separately at the |mu/sigma_u| ~ 10-100 seen in
    # this data, so form it in logs; and both formulas need R - 1 with R -> 1 in
    # exactly that heavy-truncation regime, so take it via expm1 rather than
    # subtracting (otherwise the leading digits cancel where most of the data is).
    q  <- exp(M - ll)                                  # n x K quadrature posterior
    r  <- q * (y - plogis(Z))                          # n x K: q_ik * dH_ik/dZ_ik
    logR <- outer(-0.5 * (mu / sigma_u)^2, log1m_p, "+") + 0.5 * Q^2
    Rm1 <- expm1(logR)                                 # n x K: R - 1

    attr(ll, "gradient") <- cbind(
      Xb  * rowSums(r),                                # dl/dbeta
      Xmu * rowSums(r * Rm1),                          # dl/ddelta  (= -dU/dmu)
      log_sigma_u = -rowSums(r * (U + mu * Rm1))       # dl/dlog_sigma_u
    )
    ll
  }
}

## ---- fitting wrapper ---------------------------------------------------------

fit_fractional_frontier <- function(formula_beta, formula_mu, data,
                                     n_quad = 64, start = NULL, method = "BFGS", ...) {
  vars <- unique(c(all.vars(formula_beta), all.vars(formula_mu)))
  data <- data[complete.cases(data[vars]), , drop = FALSE]

  mf <- model.frame(formula_beta, data)
  y  <- model.response(mf)
  if (any(y < 0 | y > 1)) stop("response must lie in [0, 1]")
  Xb  <- model.matrix(formula_beta, data)
  Xmu <- model.matrix(formula_mu, data)

  quad <- gauss_legendre_01(n_quad)
  ll <- make_frontier_loglik(y, Xb, Xmu, quad)

  if (is.null(start)) {
    b0 <- coef(glm(formula_beta, data = data, family = quasibinomial(link = "logit")))
    d0 <- c(0, rep(0, ncol(Xmu) - 1))
    names(b0) <- paste0("beta_", colnames(Xb))
    names(d0) <- paste0("delta_", colnames(Xmu))
    start <- c(b0, d0, log_sigma_u = log(0.3))
  }

  res <- maxLik(ll, start = start, method = method, ...)
  attr(res, "Xb_names")   <- colnames(Xb)
  attr(res, "Xmu_names")  <- colnames(Xmu)
  attr(res, "loglik_fn")  <- ll
  res
}

## ---- robust (sandwich) standard errors ---------------------------------------
#
# The Bernoulli-type log-likelihood here is a QUASI-likelihood in the Papke-Wooldridge
# sense: a working likelihood for the conditional mean, not the true density of a
# fractional y. As they note for the plain fractional logit case, that calls for
# robust (sandwich) standard errors rather than the naive inverse-Hessian ones.

vcov_robust <- function(fit) {
  ll <- attr(fit, "loglik_fn")
  if (is.null(ll)) stop("fit has no attached loglik_fn (did you fit it with fit_fractional_frontier?)")
  theta_hat <- coef(fit)

  scores <- attr(ll(theta_hat), "gradient")     # n x p per-observation score contributions
  meat   <- crossprod(scores)                   # p x p, sum_i s_i s_i'
  bread  <- solve(-hessian(fit))

  bread %*% meat %*% bread
}

summary_robust <- function(fit) {
  lmtest::coeftest(fit, vcov = vcov_robust(fit))[activePar(fit), , drop = FALSE]
}

## ---- fit separately per group (e.g. per benchmark) ---------------------------

fit_fractional_frontier_by_group <- function(formula_beta, formula_mu, data, group, ...) {
  group_var <- data[[group]]
  levels <- sort(unique(group_var[!is.na(group_var)]))

  fits <- setNames(vector("list", length(levels)), as.character(levels))
  for (lev in levels) {
    fits[[as.character(lev)]] <- fit_fractional_frontier(
      formula_beta, formula_mu, data[group_var == lev, , drop = FALSE], ...
    )
  }
  fits
}

# sigma_u at the reference date (tc = 0). With a time-varying scale the other
# logsig_* terms describe how it moves from there, so the intercept alone is the
# comparable single number.
sigma_u_hat <- function(fit) {
  cf <- coef(fit)
  i <- match("logsig_(Intercept)", names(cf))
  if (is.na(i)) i <- length(cf)     # older single-scale fits
  unname(exp(cf[i]))
}

## ---- example usage ------------------------------------------------------------

if (sys.nframe() == 0) {
  source(if (file.exists("src/paths.R")) "src/paths.R" else "paths.R")
  src_source("prepare_data.R")   # single data-prep path; no .dta is read
  d <- build_runs()
  d <- d[stats::complete.cases(d[c("acc", "lncost", "releasedate")]), ]
  d$t <- as.numeric(d$releasedate - as.Date("2023-01-01")) / 365.25

  fits <- fit_fractional_frontier_by_group(
    acc ~ lncost + t,
    # ~ t,   # dropped: t in both formula_beta and formula_mu looked like it was
    ~ 1,     # driving the sigma_u boundary/identification problem
    data = d,
    group = "benchmark",
    # delta_(Intercept) was weakly identified even as a bare constant (moves in
    # lockstep with log_sigma_u); fixing it at 0 collapses mu to 0 for everyone,
    # i.e. the classical half-normal special case of the same likelihood
    fixed = "delta_(Intercept)"
  )

  for (bench in names(fits)) {
    cat(sprintf("\n==== benchmark: %s ====\n", bench))
    print(summary_robust(fits[[bench]]))
    cat(sprintf("sigma_u (implied): %.4f\n", sigma_u_hat(fits[[bench]])))
  }
}
