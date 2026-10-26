## District estimates from the Fay-Herriot (EBLUP) and latent-field (TMB) fits of the
## province model without covariates, one panel per country (Lesotho at admin-1).
## Reads:  output/{country}/fits/district_predictions.rds and
##         province_calibration_conditional.rds (not in the release)
## Writes: figures/fig_eblup_vs_tmb_scatter.pdf
## Usage:  Rscript R/fig_eblup_vs_tmb_scatter.R (R/eblup_vs_latentfield_agreement.R reuses the fit)
suppressMessages({library(TMB); library(dplyr); library(ggplot2)})
source("R/utils.R")
if (!file.exists("tmb/mu_free.so")) compile("tmb/mu_free.cpp")
dyn.load(dynlib("tmb/mu_free"))
cts <- c("Mozambique","Zimbabwe","Malawi","Lesotho")

fit_country <- function(cty){
  c <- tolower(cty)
  dpf <- if (c=="lesotho") "district_predictions_admin1_backup.rds" else "district_predictions.rds"
  dp <- harmonise_register_units(as.data.frame(readRDS(file.path("output",c,"fits",dpf))), c) %>% filter(is.finite(z_d), is.finite(v_d), v_d>0, n_positive>0) %>%
    select(-any_of(c("Z_hat_p","psi_p","lambda_p","var_lambda_p")))   # drop NA + zero-positive hospital reporting-units
  pcc <- as.data.frame(readRDS(file.path("output",c,"fits","province_calibration_conditional.rds"))) %>%
    select(province, Z_hat_p, psi_p) %>% filter(is.finite(Z_hat_p), psi_p>0)
  d <- dp %>% inner_join(pcc, by="province")
  us <- sort(unique(d$province)); uidx <- match(d$province, us)
  su <- d %>% distinct(province,.keep_all=TRUE) %>% arrange(match(province,us))
  Zhat <- as.numeric(su$Z_hat_p); psi <- as.numeric(su$psi_p)
  X <- unname(model.matrix(~factor(province)-1, data=d)); storage.mode(X) <- "double"   # province dummies = lambda_p
  zz <- as.numeric(d$z_d); vv <- as.numeric(d$v_d)
  ## TMB (latent-field, marginal ML)
  obj <- MakeADFun(list(z=zz,v=vv,prov_d=as.integer(uidx-1),Xlam=X,Zhat=Zhat,psi=psi,w=rep(0,nrow(d)),agg=0L),
                   list(eta_p=Zhat,beta=rep(0,ncol(X)),log_sigma_u=log(0.3),u=rep(0,nrow(d))),
                   random="u", DLL="mu_free", silent=TRUE)
  nlminb(obj$par,obj$fn,obj$gr,control=list(iter.max=4000,eval.max=6000))
  pl <- obj$env$parList(par=obj$env$last.par.best)   # eta_p, beta, log_sigma_u, u
  eta_tmb <- pl$eta_p[uidx] + pl$u
  ## EBLUP (Fay-Herriot): GLS for lambda, ML for sigma_u^2, BLUP for u
  D <- d$z_d - Zhat[uidx]
  blk <- function(a) lapply(seq_along(us), function(i){ix<-which(uidx==i);n<-length(ix);solve(a*diag(n)+diag(d$v_d[ix],n)+psi[i]*matrix(1,n,n))})
  gls <- function(a){S<-blk(a);XX<-matrix(0,ncol(X),ncol(X));XD<-numeric(ncol(X))
    for(i in seq_along(us)){ix<-which(uidx==i);Xp<-X[ix,,drop=FALSE];XX<-XX+t(Xp)%*%S[[i]]%*%Xp;XD<-XD+as.numeric(t(Xp)%*%S[[i]]%*%D[ix])};solve(XX,XD)}
  nr <- function(la){a<-exp(la);S<-blk(a);b<-gls(a);v<-0
    for(i in seq_along(us)){ix<-which(uidx==i);n<-length(ix);Sp<-a*diag(n)+diag(d$v_d[ix],n)+psi[i]*matrix(1,n,n)
      r<-D[ix]-X[ix,,drop=FALSE]%*%b;v<-v+as.numeric(determinant(Sp,log=TRUE)$mod)+as.numeric(t(r)%*%S[[i]]%*%r)};0.5*v}
  a <- exp(optimize(nr,c(-9,3))$min); b <- gls(a); S <- blk(a); eta_eb <- numeric(nrow(d))
  for(i in seq_along(us)){ix<-which(uidx==i);r<-D[ix]-X[ix,,drop=FALSE]%*%b;u<-as.numeric(a*S[[i]]%*%r);eta_eb[ix]<-Zhat[i]+u}
  data.frame(country=cty, eblup=plogis(eta_eb), tmb=plogis(eta_tmb))
}

df <- bind_rows(lapply(cts, function(k) tryCatch(fit_country(k), error=function(e){cat(k,"ERR:",conditionMessage(e),"\n");NULL})))
df$country <- factor(df$country, levels=cts)
rr <- df %>% group_by(country) %>% summarise(r=cor(eblup,tmb), .groups="drop")
g <- ggplot(df, aes(100*eblup, 100*tmb)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey60") +
  geom_point(colour="#0072B2", alpha=.55, size=1.6) +
  facet_wrap(~country, scales="free", nrow=1) +
  labs(x="Fay--Herriot (EBLUP) undiagnosed prevalence (%)", y="latent-field (TMB) (%)") +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), plot.background=element_rect(fill="white",colour=NA),
        panel.background=element_rect(fill="white",colour=NA))
ggsave("figures/fig_eblup_vs_tmb_scatter.pdf", g, width=12, height=3.4)
print(rr); cat("wrote fig_eblup_vs_tmb_scatter.pdf (", nrow(df), "districts )\n")
