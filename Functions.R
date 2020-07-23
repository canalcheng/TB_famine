
# apcglmkfit
# The apclinkfit() function in APCG1 package was adapted to handle generalized linear regression
apcglmkfit <- function(r, header = F, n.risk = NA, fam = "loglin", Plot = F, pylim = c(0, 0), Scale = 1e-05, c1 = 1, c2 = 1, cc = 1, lindex = 0, lvalue = 0, amin = 1, pmin = 1, cmin = 1, n.interval = 1, stderrplot = T, pcex = 1, p0 = 1, k = 1, agapyr = 1, lambda = 0, pgapyr = 1, cgapyr = 1){
  
  # header and design matrix
  if (header == T) {
    age = floor(as.numeric(dimnames(r)[[1]]))
    period = floor(as.numeric(dimnames(r)[[2]]))
    cohort0 = period[1] - age[length(age)]
  }
  
  a <- nrow(r)
  p <- ncol(r)
  y <- c(t(r))
  x <- APCG1:::apckmat(a, p, p0 = p0, k = k)
  n.coh = ncol(x) - a - p + 2
  if (cmin > 1000) 
    cohort0 = cmin
  cohort = cohort0 + (1:n.coh - 1) * (age[3] - age[2])
  if (header == F) {
    age = amin + agapyr * (1:a - 1)
    period = pmin + pgapyr * (1:p - 1)
    cohort = cmin + cgapyr * (1:(ncol(x) - a - p + 1))
  } 
  
  x1 = x[, 1:a]
  x2 = x[, -1:-a]
  v <- eigen(t(x) %*% x)$vector
  x <- x %*% v[, -(ncol(v))]
  
  # calculate case count
  rr = c(t(round(r*n.risk/100000, 0)))
  
  if (fam == "loglin") {
    fit = glm(rr ~ x - 1 + offset(log(c(t(n.risk)))), family = poisson(link = log))
  } else if (fam == "qlik") {
    fit = glm(rr ~ x - 1 + offset(log(c(t(n.risk)))), family = quasipoisson(link = log))
  } else
  {
    stop("Wrong family of model!")
  }
  
  # calculate variance
  b = as.vector(v %*% as.matrix(c(fit$coef, 0)))
  bvar.1 = diag(summary.glm(fit)$coef[, 2]) %*% summary.glm(fit, corr = T)$corr %*% diag(summary.glm(fit)$coef[, 2])
  bvarall = v %*% cbind(rbind(bvar.1, 0), 0) %*% t(v)
  b.int = b[1]
  b.int.var = bvarall[1, 1]
  b = b[-1]
  
  a.se = sqrt(matrix(-1, 1, a - 1) %*% bvarall[2:a, 2:a] %*% matrix(-1, a - 1, 1))
  p.se = sqrt(matrix(-1, 1, p - 1) %*% bvarall[2:p - 1 + a, 2:p - 1 + a] %*% matrix(-1, p - 1, 1))
  c.se = sqrt(matrix(-1, 1, n.coh - 1) %*% bvarall[-1:-(p - 1 + a), -1:-(p - 1 + a)] %*% matrix(-1, n.coh - 1, 1))
  
  b.se = rbind(c(b.int, sqrt(b.int.var), b.int.t <- b.int/sqrt(b.int.var), 2 * (1 - pt(abs(b.int.t), df = fit$df.resid))), 
               cbind(b[1:(a - 1)], sqrt(diag(bvarall)[2:a]), tt <- b[1:(a - 1)]/sqrt(diag(bvarall)[2:a]), 2 * (1 - pt(abs(tt), df = fit$df.resid))), 
               c(-sum(b[1:(a - 1)]), a.se, tt <- -sum(b[1:(a - 1)])/a.se, 2 * (1 - pt(abs(tt), df = fit$df.resid))),
               
               cbind(b[a:(a + p - 2)], sqrt(diag(bvarall)[1 + a:(a + p - 2)]), tt <- b[a:(a + p - 2)]/sqrt(diag(bvarall)[1 + a:(a + p - 2)]), 2 * (1 - pt(abs(tt), df = fit$df.resid))),
               c(-sum(b[a:(a + p - 2)]), p.se, tt <- -sum(b[a:(a + p - 2)])/p.se, 2 * (1 - pt(abs(tt), df = fit$df.resid))),
               
               cbind(b[-1:-(a + p - 2)], sqrt(diag(bvarall)[-1:-(a + p - 1)]), tt <- b[-1:-(a + p - 2)]/sqrt(diag(bvarall)[-1:-(a + p - 1)]), 2 * (1 - pt(abs(tt), df = fit$df.resid))), 
               c(-sum(b[-1:-(a + p - 2)]), c.se, tt <- -sum(b[-1:-(a + p - 2)])/c.se, 2 * (1 - pt(abs(tt), df = fit$df.resid))))
  
  
  
  out.PCA = b.se
  
  
  diag.mu = diag(fit$fitted)
  b.var.delta = solve(t(x2) %*% diag.mu %*% x2) %*% (t(x2) %*% diag.mu %*% x1) %*% bvarall[1:a, 1:a] %*% t(x1) %*% diag.mu %*% x2 %*% solve(t(x2) %*% diag.mu %*% x2)
  bse.delta = sqrt(diag(b.var.delta))
  
  b.last.se.delta.p = sqrt(c(matrix(1, 1, p - 1) %*% b.var.delta[1:(p - 1), 1:(p - 1)] %*% matrix(1, p - 1, 1)))
  
  b.last.se.delta.c = sqrt(c(matrix(1, 1, n.coh - 1) %*% b.var.delta[-1:-(p - 1), -1:-(p - 1)] %*% matrix(1, n.coh - 1, 1)))
  
  p.delta = cbind(est <- c(b[a:(a + p - 2)], -sum(b[a:(a + p - 2)])), bse.est <- c(bse.delta[1:(p - 1)], b.last.se.delta.p), tt <- est/bse.est, 2 * (1 - pt(abs(tt), fit$df.resid)))
  
  
  c.delta = cbind(est <- c(b[-1:-(a + p - 2)], -sum(b[-1:-(a + p - 2)])), bse.est <- c(bse.delta[-1:-(p - 1)],  b.last.se.delta.c), tt <- est/bse.est, 2 * (1 - pt(abs(tt), fit$df.resid)))
  
  b.se.delta = rbind(c(b.int, sqrt(b.int.var), b.int.t <- b.int/sqrt(b.int.var), 2 * (1 - pt(abs(b.int.t), df = fit$df.resid))),
                     cbind(b[1:(a - 1)], sqrt(diag(bvarall)[2:a]), tt <- b[1:(a - 1)]/sqrt(diag(bvarall)[2:a]), 2 * (1 - pt(abs(tt), df = fit$df.resid))), 
                     c(-sum(b[1:(a - 1)]), a.se, tt <- -sum(b[1:(a - 1)])/a.se, 2 * (1 - pt(abs(tt), df = fit$df.resid))), p.delta, c.delta)
  
  b.label = dimnames(b.se)[[1]]
  b.label[1] = "Intercept"
  
  for (i in 1:a) b.label[i + 1] = paste("Age", as.character(age[i]))
  for (i in (a + 1):(a + p)) b.label[i + 1] = paste("Period",  as.character(period[i - a]))
  for (i in (a + p + 1):(a + p + 1 + n.coh - 1)) b.label[i + 1] = paste("Cohort", as.character(cohort[i - a - p]))
  
  dimnames(b.se)[[1]] = b.label
  
  out = vector("list")
  out$model = summary(fit)$call
  
  dev = matrix(c(summary(fit)$deviance, summary(fit)$df.resid), 
               1, 2)
  dimnames(dev)[[2]] = c("Deviance", "DF")
  out$deviance = dev
  out$pearson.chisq = sum((resid(fit, type = "pearson"))^2)
  p.val = 1 - pchisq(dev[1], dev[2])
  names(p.val) = "p.value"
  out$p.val = p.val
  if (fam == "qlik") 
    out$dispersion = summary(fit)$dispersion
  
  par.delta = b.se.delta
  dimnames(par.delta)[[1]] = dimnames(b.se)[[1]]
  dimnames(par.delta)[[2]] = dimnames(summary(fit)$coef)[[2]]
  out$parameter = par.delta
  
  out$a.vcov = bvarall[2:a, 2:a]
  out$p.vcov = b.var.delta[1:(p - 1), 1:(p - 1)]
  out$c.vcov = b.var.delta[-1:-(p - 1), -1:-(p - 1)]
  
  
  return(out)
}

# used by the apcglmkfit function internally
apclinkfit_run_update <- function (r, c1 = 1, c2 = 1, cc = 1, apc = 0, lambda = 0, gamma = 2, Plot = F, CIplot = T, transform = "log", panelplot = T, amin = 1, pmin = 1, cmin = 1, agapyr = 1, pgapyr = 1, cgapyr = 1, p0 = 0, k = 1, header = F, resid.out = F, ModelDiag = F, Diaglim = c(0, 0), Diag.na = c(0, 0), pplim = c(0, 0)) 
{
  if (header == T) {
    age = floor(as.numeric(dimnames(r)[[1]]))
    period = floor(as.numeric(dimnames(r)[[2]]))
    cohort0 = period[1] - age[length(age)]
    if (age[2] - age[1] != age[3] - age[2]) 
      age[1] = 2 * age[2] - age[3]
    if (period[2] - period[1] != period[3] - period[2]) 
      period[1] = 2 * period[2] - period[3]
  }
  a <- nrow(r)
  p <- ncol(r)
  if (transform == "log") 
    r = log(r)
  y <- c(t(r))
  x <- APCG1:::apckmat(a, p, p0 = p0, k = k)
  n.coh = ncol(x) - a - p + 2
  if (cmin > 1000) 
    cohort0 = cmin
  cohort = cohort0 + (1:n.coh - 1) * (age[3] - age[2])
  if (header == F) {
    age = amin + agapyr * (1:a - 1)
    period = pmin + pgapyr * (1:p - 1)
    cohort = cmin + cgapyr * (1:(ncol(x) - a - p + 1))
  }
  if (apc == 0) {
    v <- eigen(t(x) %*% x)$vector
    xx <- x %*% v[, -(ncol(v))]
    if (prod(Diag.na) > 0) {
      if (length(as.vector(Diag.na)) == 2) 
        r[Diag.na[1], Diag.na[2]] = -Inf
      else for (i in 1:(length(as.vector(Diag.na))/2)) r[Diag.na[i, 
                                                                 1], Diag.na[i, 2]] = -Inf
    }
    r = as.vector(t(r))
    r.na = (r == -Inf)
    y = r[!r.na]
    x.no.na = xx[!r.na, ]
    xx = x.no.na
    fit <- lm(y ~ xx - 1)
    b <- as.vector(v %*% as.matrix(c(fit$coef, 0)))
    aa <- b[1:(a - 1) + 1]
    pp <- b[2:p + a - 1]
    coh <- b[-1:-(a + p - 1)]
    n.coh = length(coh) + 1
    bvar.red <- diag(summary.lm(fit)$coef[, 2]) %*% summary.lm(fit, 
                                                               corr = T)$corr %*% diag(summary.lm(fit)$coef[, 2])
    bvar <- v %*% (cbind(rbind(bvar.red, 0), 0)) %*% t(v)
    alast.se <- sqrt(t(as.matrix(rep(-1, a - 1))) %*% bvar[2:a, 
                                                           2:a] %*% as.matrix(rep(-1, a - 1)))
    plast.se <- sqrt(t(as.matrix(rep(-1, p - 1))) %*% bvar[2:p - 
                                                             1 + a, 2:p - 1 + a] %*% as.matrix(rep(-1, p - 1)))
    clast.se <- sqrt(t(as.matrix(rep(-1, n.coh - 1))) %*% 
                       bvar[-1:-(p - 1 + a), -1:-(p - 1 + a)] %*% as.matrix(rep(-1, 
                                                                                n.coh - 1)))
    b.para <- cbind(c(b[1], aa, -sum(aa), pp, -sum(pp), coh, 
                      -sum(coh)), c(sqrt(diag(bvar[1:a, 1:a])), alast.se, 
                                    sqrt(diag(bvar[2:p + a - 1, 2:p + a - 1])), plast.se, 
                                    sqrt(diag(bvar[-1:-(a + p - 1), -1:-(a + p - 1)])), 
                                    clast.se))
    b.para.full = cbind(b.para, Tval <- b.para[, 1]/b.para[, 
                                                           2], 2 * (1 - pt(abs(Tval), summary(fit)$df)))
    anames = paste("Age ", as.character(age))
    pnames = paste("Period ", as.character(period))
    cnames = paste("Cohort ", as.character(cohort))
    dimnames(b.para.full)[[1]] = c("intercept", anames, pnames, 
                                   cnames)
    dimnames(b.para.full)[[2]] = c("Parameter", "Standard.Error", 
                                   "t", "P-value")
  }
  if (apc > 0) {
    xc <- constr(x[, -1], a, p, c1, c2, cc, apc)
    fit <- lm(y ~ xc)
    sol$inter <- fit$coef[1]
    if (apc == 1) {
      aa <- fit$coef[2:(a - 1)]
      aa <- c(aa[1:(c2 - 1)], cc * aa[c1], aa[-1:-(c2 - 
                                                     1)])
      pp <- fit$coef[a:(a + p - 2)]
      coh <- fit$coef[-1:-(a + p - 2)]
    }
    if (apc == 2) {
      aa <- fit$coef[2:a]
      pp <- fit$coef[(a + 1):(a + p - 2)]
      pp <- c(pp[1:(c2 - 1)], cc * pp[c1], pp[-1:-(c2 - 
                                                     1)])
      coh <- fit$coef[-1:-(a + p - 2)]
    }
    if (apc == 3) {
      aa <- fit$coef[2:a]
      pp <- fit$coef[(a + 1):(a + p - 1)]
      coh <- fit$coef[(a + p):length(fit$coef)]
      coh <- c(coh[1:(c2 - 1)], cc * coh[c1], coh[-1:-(c2 - 
                                                         1)])
    }
  }
  if (CIplot == T) {
    alim = range(b.para[1:a + 1, 1] + 1.96 * b.para[1:a + 
                                                      1, 2], b.para[1:a + 1, 1] - 1.96 * b.para[1:a + 1, 
                                                                                                2])
    plim = range(b.para[1:p + a + 1, 1] + 1.96 * b.para[1:p + 
                                                          a + 1, 2], b.para[1:p + a + 1, 1] - 1.96 * b.para[1:p + 
                                                                                                              a + 1, 2])
    clim = range(b.para[-1:-(a + p + 1), 1] + 1.96 * b.para[-1:-(a + 
                                                                   p + 1), 2], b.para[-1:-(a + p + 1), 1] - 1.96 * b.para[-1:-(a + 
                                                                                                                                 p + 1), 2])
  }
  if (CIplot == F) {
    alim = range(b.para[1:a + 1, 1])
    plim = range(b.para[1:p + a + 1, 1])
    clim = range(b.para[-1:-(a + p + 1), 1])
  }
  if (sum(abs(pplim)) > 0) {
    alim = plim = clim = pplim
  }
  if (panelplot == T) 
    par(mfrow = c(2, 2))
  if (Plot == T) {
    plot(age, c(aa, -sum(aa)), ylab = "effect", type = "l", 
         ylim = alim)
    points(age, c(aa, -sum(aa)))
    title(main = "Age effect")
    if (CIplot == T) {
      lines(age, b.para[1:a + 1, 1] - 1.96 * b.para[1:a + 
                                                      1, 2], lty = 2)
      lines(age, b.para[1:a + 1, 1] + 1.96 * b.para[1:a + 
                                                      1, 2], lty = 2)
    }
    plot(period, c(pp, -sum(pp)), ylab = "effect", type = "l", 
         ylim = plim)
    points(period, c(pp, -sum(pp)))
    title(main = "Period effect")
    if (CIplot == T) {
      lines(period, b.para[1:p + a + 1, 1] - 1.96 * b.para[1:p + 
                                                             a + 1, 2], lty = 2)
      lines(period, b.para[1:p + a + 1, 1] + 1.96 * b.para[1:p + 
                                                             a + 1, 2], lty = 2)
    }
    if (panelplot == T) 
      par(mfrow = c(2, 1), new = T, fig = c(0, 1, 0, 0.56), 
          cex = 0.85)
    plot(cohort, c(coh, -sum(coh)), ylab = "effect", type = "l", 
         ylim = clim)
    points(cohort, c(coh, -sum(coh)))
    title(main = "Cohort effect")
    if (CIplot == T) {
      lines(cohort, b.para[-1:-(p + a + 1), 1] - 1.96 * 
              b.para[-1:-(p + a + 1), 2], lty = 2)
      lines(cohort, b.para[-1:-(p + a + 1), 1] + 1.96 * 
              b.para[-1:-(p + a + 1), 2], lty = 2)
    }
    par(mfrow = c(1, 1))
  }
  if (ModelDiag == T) {
    plot(fit$fitted, fit$residual, xlab = "fitted value", 
         ylab = "residual")
    abline(a = 0, b = 0, col = 2)
    title(main = "Residual plot for linear model diagnostics")
    if (var(Diaglim) > 0) {
      abline(a = Diaglim[1], b = 0, lty = 2, col = 4)
      abline(a = Diaglim[2], b = 0, lty = 2, col = 4)
      obs.id = c((1:length(y))[fit$resid >= max(Diaglim)], 
                 (1:length(y))[fit$resid <= min(Diaglim)])
      n.obs.id = length(obs.id)
      outlier.obs = matrix(NA, n.obs.id, 4)
      for (i in 1:n.obs.id) {
        outlier.obs[i, 1] = fit$resid[obs.id[i]]
        outlier.obs[i, 2] = y[obs.id[i]]
        outlier.obs[i, 3] = 1 + floor((obs.id[i] - 1)/p)
        outlier.obs[i, 4] = obs.id[i] - p * floor(obs.id[i]/p)
        if (outlier.obs[i, 4] == 0) 
          outlier.obs[i, 4] = p
      }
    }
  }
  Rsquare = c(summary(fit)$r.squared, summary(fit)$adj.r.squared, 
              ncol(v) - 1)
  names(Rsquare) = c("Rsquare", "Adjusted.Rsquare", "Model df")
  Varcomp = c(summary(fit)$sig, summary(fit)$df[2])
  names(Varcomp) = c("Variance Component", "Residual df")
  sol <- vector("list")
  sol$model = fit$call
  sol$scale = transform
  sol$Rsquare = Rsquare
  sol$varcomp_df = Varcomp
  sol$parameters = b.para.full
  if (ModelDiag == T && var(Diaglim) > 0) 
    sol$outlier = outlier.obs
  sol$fitted = fit$fitted
  sol$residual = fit$residual
  return(sol)
}

# used by the apcglmkfit function internally
apclinkfit_update <- function (r, Plot = F, CIplot = T, transform = "log", amin = 1,  pmin = 1, cmin = 1, agapyr = 1, lambda = 0, pgapyr = 1, cgapyr = 1,  p0 = 0, k = 1, header = F, ModelDiag = F, pplim = c(0, 0)) 
{
  out <- apclinkfit_run_update(r = r, apc = 0, lambda = 0, Plot = Plot, CIplot = CIplot, transform = transform, amin = amin,pmin = pmin, cmin = cmin, agapyr = agapyr, pgapyr = pgapyr, cgapyr = cgapyr, p0 = p0, k = k, header = header, resid.out = F, 
ModelDiag = ModelDiag, Diaglim = c(0, 0), pplim = pplim)
  return(out)
}



#' Format APC result
#'
#' @param apc.result apc mode
#'
#' @return a data frame with 4 cols: Type, X, parameter, and standard deviation
#' @export
#'
#' @examples
reshape.APC.result <- function(apc.result)
{
  formatted.pre <- data.frame(apc.result$parameter)
  
  formatted <- data.frame(format = row.names(formatted.pre)) %>% separate(format, into = c("Type","X")) %>% mutate(parameter = formatted.pre$Estimate, sd = formatted.pre$Std..Error)
  
  formatted
}




#' Get ggplot color
#'
#' @param n The number of colors to get
#'
#' @return a string of length n representing the colors
#' @export
#'
#' @examples
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


province.cohort.smooth <- function(result.mf, F1.start = 1957, F1.end = 1963, F2.start = 1978, F2.end = 1990, method = c("GAM","loess"))
{
  result.smooth <- NULL
  
  for(i in 1:2)
  {
    current.result <- result.mf[result.mf$Sex == unique(result.mf$Sex)[i],]
    current.cohort <- current.result[current.result$Type == "Cohort", c("X","parameter", "sd")]
    colnames(current.cohort) <- c("Cohort","Cohort.effect.e","Cohort.std")
    current.cohort$Sex <-  unique(result.mf$Sex)[i]
    
    # use span from expectation
    # optimal df for GAM
    current.cohort.rm <- current.cohort[!current.cohort$Cohort %in% c(seq(F1.start, F1.end,3), seq(F2.start, F2.end,3)),]
    current.cohort.rm <- current.cohort.rm[c(-1,-2, -nrow(current.cohort.rm),-(nrow(current.cohort.rm)-1)),]
    # current.cohort.rm <- current.cohort.rm[c(-1, -nrow(current.cohort.rm)),]
    
    if(method == "GAM")
    {
      current.gam <- gam(Cohort.effect.e ~ s(Cohort), data = current.cohort.rm)
      current.cohort$cohort.predict <- predict(current.gam, newdata = data.frame(Cohort = current.cohort$Cohort))
    }
    
    if(method == "loess")
    {
      current.cohort.loess <- loess.as(current.cohort.rm$Cohort, current.cohort.rm$Cohort.effect.e, criterion = "gcv")
      current.cohort$cohort.predict <- predict(current.cohort.loess, newdata = data.frame(x = current.cohort$Cohort))
    }
    
    result.smooth <- rbind(result.smooth, current.cohort)
  }
  
  result.smooth
}


province.RR <- function(result.mf, result.mf.vcov, F1.start = 1957, F1.end = 1963, F2.start = 1978, F2.end = 1990, F2 = 1981, n.rep = 100, method = c("GAM","loess", "linear"), Sichuan.pop.long, Sichuan.case.long, ncores = 7)
{
  result.RR <- NULL
  
  writeLines("","parameterset.txt")
  cl <- makeCluster(ncores, type = "SOCK")
  registerDoSNOW(cl)
  
  for(i in 1:2)
  {
    sink("parameterset.txt",append=TRUE)
    
    current.result <- result.mf[result.mf$Sex == unique(result.mf$Sex)[i],]
    current.age <- current.result[current.result$Type == "Age", c("X","parameter", "sd")]
    colnames(current.age) <- c("Age","Age.effect.e","Age.std")
    current.period <- current.result[current.result$Type == "Period", c("X","parameter", "sd")]
    colnames(current.period) <- c("DiagnoseYear","Period.effect.e","Period.std")
    current.cohort <- current.result[current.result$Type == "Cohort", c("X","parameter", "sd")]
    colnames(current.cohort) <- c("Cohort","Cohort.effect.e","Cohort.std")
    
    # use span from expectation
    # optimal df for GAM
    current.cohort.rm <- current.cohort[!current.cohort$Cohort %in% c(seq(F1.start, F1.end,3), seq(F2.start, F2.end,3)),]
    current.cohort.rm <- current.cohort.rm[c(-1, -nrow(current.cohort.rm)),]
    
    if(method == "GAM")
    {
      current.gam <- gam(Cohort.effect.e ~ s(Cohort), data = current.cohort.rm)
      optimal.sp <- current.gam$sp
    }
    
    if(method == "loess")
    {
      current.cohort.loess <- loess.as(current.cohort.rm$Cohort, current.cohort.rm$Cohort.effect.e, criterion = "gcv")
      optimal.span <- current.cohort.loess$pars$span
    }
    
    result <- foreach(j = 1:n.rep, .combine = "rbind", .packages = c("fANCOVA","tidyverse", "splines", "mgcv")) %dopar%
    {
      cat("i = ",i, ", j = ",j,"\n")
      current.RR <- data.frame(rep = j, Sex = unique(result.mf$Sex)[i], F1.RR = rep(NA, n.rep), F2.RR = rep(NA, n.rep), F1.avertcase = rep(NA, n.rep), F2.avertcase = rep(NA, n.rep))
      
      age.first <- mvnfast::rmvn(n = 1, mu = current.age$Age.effect.e[-length(current.age$Age.effect.e)], sigma = result.mf.vcov[[i]]$a.vcov)[1,]
      age.lastone <- -sum(age.first)
      current.age$Age.effect <- c(age.first, age.lastone)
      
      period.first <- mvnfast::rmvn(n = 1, mu = current.period$Period.effect.e[-length(current.period$Period.effect.e)], sigma = result.mf.vcov[[i]]$p.vcov)[1,]
      period.lastone <- -sum(period.first)
      current.period$Period.effect <- c(period.first, period.lastone)
      
      cohort.first <- mgcv::rmvn(n = 1, mu = current.cohort$Cohort.effect.e[-length(current.cohort$Cohort.effect.e)], V = result.mf.vcov[[i]]$c.vcov)
      cohort.lastone <- -sum(cohort.first)
      current.cohort$Cohort.effect <- c(cohort.first, cohort.lastone)
      
      current.cohort.rm <- current.cohort[!current.cohort$Cohort %in% c(seq(F1.start, F1.end,3), seq(F2.start, F2.end,3)),]
      current.cohort.rm <- current.cohort.rm[c(-1,-nrow(current.cohort.rm)),]
      
      if(method == "GAM")
      {
        current.cohort.gam <- gam(Cohort.effect ~ s(Cohort), data = current.cohort.rm, sp = optimal.sp)
        # plot(current.cohort.gam)
        current.cohort.gam.predict <- predict(current.cohort.gam, newdata = current.cohort, se.fit = TRUE)
        current.cohort$cohort.predict <- current.cohort.gam.predict$fit
        current.cohort$cohort.predict.sd <- current.cohort.gam.predict$se.fit
      }
      
      if(method == "loess")
      {
        current.cohort.loess <- loess.as(current.cohort.rm$Cohort, current.cohort.rm$Cohort.effect, user.span = optimal.span)
        current.cohort.loess.predict <- predict(current.cohort.loess, newdata = data.frame(x = current.cohort$Cohort), se = TRUE)
        current.cohort$cohort.predict <- current.cohort.loess.predict$fit
        current.cohort$cohort.predict.sd <- current.cohort.loess.predict$se.fit
      }
      
      if(method == "linear")
      {
        current.cohort$cohort.predict <- current.cohort$Cohort.effect
        # F1 interpolate
        F1.dat <- data.frame(x = c(F1.start - 3, F1.end + 3), y = c(current.cohort$Cohort.effect[current.cohort$Cohort == F1.start - 3], current.cohort$Cohort.effect[current.cohort$Cohort == F1.end + 3]))
        F1.lm <- lm(y ~ x, F1.dat)
        current.cohort$cohort.predict[current.cohort$Cohort %in% seq(F1.start, F1.end, 3)] <- predict(F1.lm, newdata = data.frame(x = seq(F1.start, F1.end, 3)))
        # F2 interpolate
        F2.dat <- data.frame(x = c(F2.start - 3, F2.end + 3), y = c(current.cohort$Cohort.effect[current.cohort$Cohort == F2.start - 3], current.cohort$Cohort.effect[current.cohort$Cohort == F2.end + 3]))
        F2.lm <- lm(y ~ x, F2.dat)
        current.cohort$cohort.predict[current.cohort$Cohort %in% seq(F2.start, F2.end, 3)] <- predict(F2.lm, newdata = data.frame(x = seq(F2.start, F2.end, 3)))
        
        current.cohort$cohort.predict.sd <- 0
      }
      
      cohort.n <- (F2.end - F2.start)/3 + 1
      
      famine.cohort.F1 <- data.frame(Age = rep(seq(46,58,3), each = 3), DiagnoseYear = 2005:2019, Cohort = (F1.start + F1.end)/2)[-15,]
      famine.cohort.F1$Age <- c(2006,2006,2006, 2009,2009,2009, 2012,2012,2012, 2015,2015,2015,2018,2018) - (F1.start + F1.end)/2
      
      # one cohort F2
      famine.cohort.F2 <- data.frame(DiagnoseYear = 2005:2018, Cohort = F2, period.collapse = c(2006,2006,2006, 2009,2009,2009, 2012,2012,2012, 2015,2015,2015,2018,2018))
      famine.cohort.F2$Age <- famine.cohort.F2$period.collapse - famine.cohort.F2$Cohort
      
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.age, by = "Age")
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.period, by = "DiagnoseYear")
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.cohort, by = "Cohort")
      famine.cohort.F1$intercept <- current.result$parameter[current.result$Type == "Intercept"]
      famine.cohort.F1$intercept.std <- current.result$sd[current.result$Type == "Intercept"]
      
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.age, by = "Age")
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.period, by = "DiagnoseYear")
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.cohort, by = "Cohort")
      famine.cohort.F2$intercept <- current.result$parameter[current.result$Type == "Intercept"]
      famine.cohort.F2$intercept.std <- current.result$sd[current.result$Type == "Intercept"]
      
      # population 
      current.pop <- Sichuan.pop.long[Sichuan.pop.long$Sex == unique(result.mf$Sex)[i],]
      
      # observed case
      current.case <- Sichuan.case.long[Sichuan.case.long$Sex == unique(result.mf$Sex)[i],]
      
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.pop[,c("DiagnoseYear","Age","Population")], by = c("Age", "DiagnoseYear"))
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.case[,c("DiagnoseYear","Age","CaseCount")], by = c("Age", "DiagnoseYear"))
      
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.pop[,c("DiagnoseYear","Age","Population")], by = c("DiagnoseYear", "Age"))
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.case[,c("DiagnoseYear","Age","CaseCount")], by = c("DiagnoseYear", "Age"))
      
      famine.cohort.F1.sample <- matrix(NA, nrow(famine.cohort.F1), n.rep)
      famine.cohort.F2.sample <- matrix(NA, nrow(famine.cohort.F2), n.rep)
      
      for(k in 1:n.rep)
      {
        famine.cohort.F1.sample[,k] <- rnorm(nrow(famine.cohort.F1), mean = famine.cohort.F1$cohort.predict, sd = famine.cohort.F1$cohort.predict.sd)
        famine.cohort.F2.sample[,k] <- rnorm(nrow(famine.cohort.F2), mean = famine.cohort.F2$cohort.predict, sd = famine.cohort.F2$cohort.predict.sd)
      }
      
      current.RR$F1.RR <- apply(famine.cohort.F1.sample, 2, function(x) sum(famine.cohort.F1$CaseCount)/sum(exp(famine.cohort.F1$Age.effect + famine.cohort.F1$Period.effect + x + famine.cohort.F1$intercept)*famine.cohort.F1$Population))  
      current.RR$F1.avertcase <- apply(famine.cohort.F1.sample, 2, function(x) sum(famine.cohort.F1$CaseCount)-round(sum(exp(famine.cohort.F1$Age.effect + famine.cohort.F1$Period.effect + x + famine.cohort.F1$intercept)*famine.cohort.F1$Population)))  
      current.RR$F2.RR <- apply(famine.cohort.F2.sample, 2, function(x) sum(famine.cohort.F2$Case)/sum(exp(famine.cohort.F2$Age.effect + famine.cohort.F2$Period.effect + x + famine.cohort.F2$intercept)*famine.cohort.F2$Population)) 
      current.RR$F2.avertcase <- apply(famine.cohort.F2.sample, 2, function(x) sum(famine.cohort.F2$Case)- round(sum(exp(famine.cohort.F2$Age.effect + famine.cohort.F2$Period.effect + x + famine.cohort.F2$intercept)*famine.cohort.F2$Population))) 
      
      current.RR
    }
    
    result.RR <- rbind(result.RR, result)
    
    sink()
  }
  
  stopCluster(cl)
  
  sink()
  
  result.RR
}



pref.RR <- function(result.mf, result.mf.vcov, F1.start = 1957, F1.end = 1963, F2.start = rep(1978,21), F2.end = rep(1990,21), F2 = rep(1981,21), n.rep = 100, method = c("GAM","loess"), Sichuan.pop.long, Sichuan.case.long, collapse = FALSE, ncore = 7)
{
  result.mf$X <- as.numeric(result.mf$X)
  
  result.RR <- NULL
  
  writeLines("","parameterset.txt")
  cl <- makeCluster(ncore, type = "SOCK")
  registerDoSNOW(cl)
  
  for(i in 1:length(unique(result.mf$Pref)))
  {
    sink("parameterset.txt",append=TRUE)
    
    current.result <- result.mf[result.mf$Pref == unique(result.mf$Pref)[i],]
    current.age <- current.result[current.result$Type == "Age", c("X","parameter", "sd")]
    colnames(current.age) <- c("Age","Age.effect.e","Age.std")
    current.period <- current.result[current.result$Type == "Period", c("X","parameter", "sd")]
    colnames(current.period) <- c("DiagnoseYear","Period.effect.e","Period.std")
    current.cohort <- current.result[current.result$Type == "Cohort", c("X","parameter", "sd")]
    colnames(current.cohort) <- c("Cohort","Cohort.effect.e","Cohort.std")
    
    
    # use span from expectation
    # optimal df for GAM
    current.cohort.rm <- current.cohort[!current.cohort$Cohort %in% c(seq(F1.start, F1.end,3), seq(F2.start[i], F2.end[i],3)),]
    current.cohort.rm <- current.cohort.rm[c(-1,-2, -nrow(current.cohort.rm),-(nrow(current.cohort.rm)-1)),]
    
    if(method == "GAM")
    {
      current.gam <- gam(Cohort.effect.e ~ s(Cohort), data = current.cohort.rm)
      optimal.sp <- current.gam$sp
    }
    
    if(method == "loess")
    {
      current.cohort.loess <- loess.as(current.cohort.rm$Cohort, current.cohort.rm$Cohort.effect.e, criterion = "gcv")
      optimal.span <- current.cohort.loess$pars$span
    }
    
    result <- foreach(j = 1:n.rep, .combine = "rbind", .packages = c("fANCOVA","tidyverse", "splines", "mgcv")) %dopar%
    {
      cat("i = ",i, ", j = ",j,"\n")
      current.RR <- data.frame(rep = j, Pref = unique(result.mf$Pref)[i], F1.RR = rep(NA, n.rep), F2.RR = rep(NA, n.rep))
      
      age.first <- mvnfast::rmvn(n = 1, mu = current.age$Age.effect.e[-length(current.age$Age.effect.e)], sigma = result.mf.vcov[[i]]$a.vcov)[1,]
      age.lastone <- -sum(age.first)
      current.age$Age.effect <- c(age.first, age.lastone)
      
      period.first <- mvnfast::rmvn(n = 1, mu = current.period$Period.effect.e[-length(current.period$Period.effect.e)], sigma = result.mf.vcov[[i]]$p.vcov)[1,]
      period.lastone <- -sum(period.first)
      current.period$Period.effect <- c(period.first, period.lastone)
      
      cohort.first <- mgcv::rmvn(n = 1, mu = current.cohort$Cohort.effect.e[-length(current.cohort$Cohort.effect.e)], V = result.mf.vcov[[i]]$c.vcov)
      cohort.lastone <- -sum(cohort.first)
      current.cohort$Cohort.effect <- c(cohort.first, cohort.lastone)
      
      current.cohort.rm <- current.cohort[!current.cohort$Cohort %in% c(seq(F1.start, F1.end,3), seq(F2.start[i], F2.end[i],3)),]
      current.cohort.rm <- current.cohort.rm[c(-1, -nrow(current.cohort.rm)),]
      
      if(method == "GAM")
      {
        current.cohort.gam <- gam(Cohort.effect ~ s(Cohort), data = current.cohort.rm, sp = optimal.sp)
        # plot(current.cohort.gam)
        current.cohort.gam.predict <- predict(current.cohort.gam, newdata = current.cohort, se.fit = TRUE)
        current.cohort$cohort.predict <- current.cohort.gam.predict$fit
        current.cohort$cohort.predict.sd <- current.cohort.gam.predict$se.fit
      }
      
      if(method == "loess")
      {
        current.cohort.loess <- loess.as(current.cohort.rm$Cohort, current.cohort.rm$Cohort.effect, user.span = optimal.span)
        current.cohort.loess.predict <- predict(current.cohort.loess, newdata = data.frame(x = current.cohort$Cohort), se = TRUE)
        current.cohort$cohort.predict <- current.cohort.loess.predict$fit
        current.cohort$cohort.predict.sd <- current.cohort.loess.predict$se.fit
      }
      
      cohort.n <- (F2.end[i] - F2.start[i])/3 + 1
      
      if(collapse)
      {
        famine.cohort.F1 <- data.frame(Age = seq(46,58,3), DiagnoseYear = seq(2006,2018,3), Cohort = 1960)
      } else
      {
        famine.cohort.F1 <- data.frame(Age = rep(seq(46,58,3), each = 3), DiagnoseYear = 2005:2019, Cohort = 1960)[-15,]
      }
      
      if(collapse)
      {
        famine.cohort.F2 <- data.frame(DiagnoseYear = seq(2006,2018,3), Cohort = F2[i])
        famine.cohort.F2$Age <- famine.cohort.F2$DiagnoseYear - famine.cohort.F2$Cohort
      } else
      {
        famine.cohort.F2 <- data.frame(DiagnoseYear = 2005:2018, Cohort = F2[i], period.collapse = c(2006,2006,2006, 2009,2009,2009, 2012,2012,2012, 2015,2015,2015,2018,2018))
        famine.cohort.F2$Age <- famine.cohort.F2$period.collapse - famine.cohort.F2$Cohort
      }
      
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.age, by = "Age")
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.period, by = "DiagnoseYear")
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.cohort, by = "Cohort")
      famine.cohort.F1$intercept <- current.result$parameter[current.result$Type == "Intercept"]
      famine.cohort.F1$intercept.std <- current.result$sd[current.result$Type == "Intercept"]
      
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.age, by = "Age")
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.period, by = "DiagnoseYear")
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.cohort, by = "Cohort")
      famine.cohort.F2$intercept <- current.result$parameter[current.result$Type == "Intercept"]
      famine.cohort.F2$intercept.std <- current.result$sd[current.result$Type == "Intercept"]
      
      # population 
      current.pop <- Sichuan.pop.long[Sichuan.pop.long$Pref == unique(result.mf$Pref)[i],]
      
      # observed case
      current.case <- Sichuan.case.long[Sichuan.case.long$Pref == unique(result.mf$Pref)[i],]
      
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.pop[,c("DiagnoseYear","Age","Population")], by = c("Age", "DiagnoseYear"))
      famine.cohort.F1 <- left_join(famine.cohort.F1, current.case[,c("DiagnoseYear","Age","CaseCount")], by = c("Age", "DiagnoseYear"))
      
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.pop[,c("DiagnoseYear","Age","Population")], by = c("DiagnoseYear", "Age"))
      famine.cohort.F2 <- left_join(famine.cohort.F2, current.case[,c("DiagnoseYear","Age","CaseCount")], by = c("DiagnoseYear", "Age"))
      
      famine.cohort.F1.sample <- matrix(NA, nrow(famine.cohort.F1), n.rep)
      famine.cohort.F2.sample <- matrix(NA, nrow(famine.cohort.F2), n.rep)
      
      for(k in 1:n.rep)
      {
        famine.cohort.F1.sample[,k] <- rnorm(nrow(famine.cohort.F1), mean = famine.cohort.F1$cohort.predict, sd = famine.cohort.F1$cohort.predict.sd)
        famine.cohort.F2.sample[,k] <- rnorm(nrow(famine.cohort.F2), mean = famine.cohort.F2$cohort.predict, sd = famine.cohort.F2$cohort.predict.sd)
      }
      
      current.RR$F1.RR <- apply(famine.cohort.F1.sample, 2, function(x) sum(famine.cohort.F1$CaseCount)/sum(exp(famine.cohort.F1$Age.effect + famine.cohort.F1$Period.effect + x + famine.cohort.F1$intercept)*famine.cohort.F1$Population))  
      current.RR$F2.RR <- apply(famine.cohort.F2.sample, 2, function(x) sum(famine.cohort.F2$Case)/sum(exp(famine.cohort.F2$Age.effect + famine.cohort.F2$Period.effect + x + famine.cohort.F2$intercept)*famine.cohort.F2$Population)) 
      
      current.RR
    }
    
    result.RR <- rbind(result.RR, result)
    
    sink()
  }
  
  stopCluster(cl)
  
  sink()
  
  result.RR
}



pref.cohort.smooth <- function(result.mf, F1.start = 1957, F1.end = 1963, F2.start = rep(1978,21), F2.end = rep(1990,21), method = c("GAM","loess"))
{
  result.smooth <- NULL
  
  result.mf$X <- as.numeric(result.mf$X)
  
  for(i in 1:21)
  {
    current.result <- result.mf[result.mf$Pref == unique(result.mf$Pref)[i],]
    current.cohort <- current.result[current.result$Type == "Cohort", c("X","parameter", "sd")]
    colnames(current.cohort) <- c("Cohort","Cohort.effect.e","Cohort.std")
    current.cohort$Pref <-  unique(result.mf$Pref)[i]
    
    # use span from expectation
    # optimal df for GAM
    current.cohort.rm <- current.cohort[!current.cohort$Cohort %in% c(seq(F1.start, F1.end,3), seq(F2.start[i], F2.end[i],3)),]
    current.cohort.rm <- current.cohort.rm[c(-1,-2, -nrow(current.cohort.rm),-(nrow(current.cohort.rm)-1)),]
    
    if(method == "GAM")
    {
      current.gam <- gam(Cohort.effect.e ~ s(Cohort), data = current.cohort.rm)
      current.cohort$cohort.predict <- predict(current.gam, newdata = data.frame(Cohort = current.cohort$Cohort))
    }
    
    if(method == "loess")
    {
      current.cohort.loess <- loess.as(current.cohort.rm$Cohort, current.cohort.rm$Cohort.effect.e, criterion = "gcv")
      current.cohort$cohort.predict <- predict(current.cohort.loess, newdata = data.frame(x = current.cohort$Cohort))
    }
    
    result.smooth <- rbind(result.smooth, current.cohort)
  }
  
  result.smooth
}
