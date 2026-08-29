library(data.table)
library(fixest)

d <- fread("analysis_work/common_matched.csv")
setDT(d)

sep_safe <- function(x, yvar = "bank_musd") {
  repeat {
    n0 <- nrow(x)
    for (g in c("pair_id", "op_year", "dp_year")) {
      keep <- x[, .(ysum = sum(get(yvar), na.rm = TRUE)), by = g][ysum > 0, get(g)]
      x <- x[get(g) %in% keep]
    }
    if (nrow(x) == n0) break
  }
  x
}

fit_pair <- function(x) {
  x <- sep_safe(copy(x), "bank_musd")
  fb <- fepois(bank_musd ~ ipd_z | pair_id + op_year + dp_year,
               data=x, vcov=~undirected_pair, fixef.rm="perfect_fit", notes=FALSE)
  ft <- fepois(trade_musd ~ ipd_z | pair_id + op_year + dp_year,
               data=x, vcov=~undirected_pair, fixef.rm="perfect_fit", notes=FALSE)
  data.table(
    beta_bank=coef(fb)[["ipd_z"]], se_bank=se(fb)[["ipd_z"]], n_bank=nobs(fb),
    beta_trade=coef(ft)[["ipd_z"]], se_trade=se(ft)[["ipd_z"]], n_trade=nobs(ft),
    n_raw=nrow(x), pairs=uniqueN(x$pair_id), clusters=uniqueN(x$undirected_pair)
  )
}

fb <- fepois(bank_musd ~ ipd_z | pair_id + op_year + dp_year,
             data=d, vcov=~undirected_pair, fixef.rm="perfect_fit", notes=FALSE)
ft <- fepois(trade_musd ~ ipd_z | pair_id + op_year + dp_year,
             data=d, vcov=~undirected_pair, fixef.rm="perfect_fit", notes=FALSE)
main <- data.table(outcome=c("Interbank","Trade"),
                   beta=c(coef(fb)[["ipd_z"]],coef(ft)[["ipd_z"]]),
                   se=c(se(fb)[["ipd_z"]],se(ft)[["ipd_z"]]),
                   nobs=c(nobs(fb),nobs(ft)))
main[, pct := 100*(exp(beta)-1)]
main[, lo95_pct := 100*(exp(beta-1.96*se)-1)]
main[, hi95_pct := 100*(exp(beta+1.96*se)-1)]
fwrite(main,"analysis_work/main_estimates.csv")

roll <- rbindlist(lapply(1987:2020,function(e) {
  x <- d[year >= e-9 & year <= e]
  if(nrow(x)<500) return(NULL)
  z <- tryCatch(fit_pair(x),error=function(err) NULL)
  if(is.null(z)) return(NULL)
  z[, `:=`(start=e-9,end=e)]
  z
}),fill=TRUE)
for(v in c("bank","trade")) {
  set(roll, j=paste0("pct_",v), value=100*(exp(roll[[paste0("beta_",v)]])-1))
  set(roll, j=paste0("lo95_",v), value=100*(exp(roll[[paste0("beta_",v)]]-1.96*roll[[paste0("se_",v)]])-1))
  set(roll, j=paste0("hi95_",v), value=100*(exp(roll[[paste0("beta_",v)]]+1.96*roll[[paste0("se_",v)]])-1))
}
fwrite(roll,"analysis_work/rolling_estimates.csv")

periods <- list(
  "1978-1989"=c(1978,1989),
  "1990-2001"=c(1990,2001),
  "2002-2008"=c(2002,2008),
  "2009-2013"=c(2009,2013),
  "2014-2020"=c(2014,2020)
)
per <- rbindlist(lapply(names(periods),function(nm) {
  rr <- periods[[nm]]; x <- d[year>=rr[1]&year<=rr[2]]
  if(nrow(x)<300) return(NULL)
  z <- tryCatch(fit_pair(x),error=function(err) NULL)
  if(is.null(z)) return(NULL)
  z[, period:=nm]
  z
}),fill=TRUE)
for(v in c("bank","trade")) set(per, j=paste0("pct_",v), value=100*(exp(per[[paste0("beta_",v)]])-1))
fwrite(per,"analysis_work/period_estimates.csv")

direct <- d[!is.na(reported_claim_musd)]
if(nrow(direct)>500) {
  direct[, bank_direct := reported_claim_musd]
  direct <- sep_safe(direct, "bank_direct")
  fd <- fepois(bank_direct ~ ipd_z | pair_id + op_year + dp_year,
               data=direct,vcov=~undirected_pair,fixef.rm="perfect_fit",notes=FALSE)
  rob <- data.table(spec="Reported claims only", beta=coef(fd)[["ipd_z"]],
                    se=se(fd)[["ipd_z"]],nobs=nobs(fd))
  rob[,pct:=100*(exp(beta)-1)]
  fwrite(rob,"analysis_work/direct_robustness.csv")
}

capture.output(summary(fb), file="analysis_work/model_bank.txt")
capture.output(summary(ft), file="analysis_work/model_trade.txt")
print(main)
