#!/usr/bin/env Rscript
# ELSA-style discovery: single-PGS delta-R2; multivariable elastic net with
# 10x10 repeated cross-validation, 10 lambda candidates, 1,000 bootstrap
# resamples; nonlinear longitudinal GAMMs (>=3 observations/person).
# Run on authorised participant-level data only; do not commit data/results.

packages <- c('dplyr', 'caret', 'glmnet', 'gamm4', 'mgcv')
not_installed <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(not_installed)) stop('Install missing packages: ', paste(not_installed, collapse=', '))

cfg <- list(
  frailty_file = 'data/frailty_long.csv', # id,age,sex,wave,frailty
  pgs_file     = 'data/pgs.csv',          # id,PGS_1,...,PGS_26
  pcs_file     = 'data/pcs.csv',          # id,PC_1,...,PC_10
  output_dir   = 'results',
  pgs_prefix   = 'PGS_',
  pc_prefix    = 'PC_',
  n_pcs        = 10,             # set to 4 for LBC-style replication
  min_age      = 65,
  n_pgs_expected = 26,        # change if reproducing a different analysis
  sqrt_frailty = TRUE,       # ELSA method; set FALSE if already transformed
  standardise_frailty = TRUE,# standardises transformed FI in this cohort
  alpha        = 0.5,        # elastic net mixing parameter
  seed         = 123,
  cv_folds     = 10,
  cv_repeats   = 10,
  n_lambdas    = 10,
  n_boot       = 1000,
  minimum_long_visits = 3,
  spline_k     = 5,
  # NULL: model trajectories for all PGSs. Otherwise a vector of PGS column names.
  longitudinal_pgs = NULL
)

read_csv <- function(path) {
  if (!file.exists(path)) stop('Missing input file: ', path)
  read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
}
fi <- read_csv(cfg$frailty_file)
pgs <- read_csv(cfg$pgs_file)
pcs <- read_csv(cfg$pcs_file)
for (nm in c('id','age','sex','wave','frailty')) if (!nm %in% names(fi)) stop('Missing frailty field: ',nm)
for (dat in list(pgs,pcs)) if (!'id' %in% names(dat)) stop('Each input must have id')
if (anyDuplicated(pgs$id) || anyDuplicated(pcs$id)) stop('PGS/PC tables need one row per person')
if (anyDuplicated(fi[c('id','wave')])) stop('Frailty table has duplicate id-wave rows')
pc_names <- paste0(cfg$pc_prefix, seq_len(cfg$n_pcs))
if (!all(pc_names %in% names(pcs))) stop('Missing required PC columns')
pgs_names <- grep(paste0('^', cfg$pgs_prefix), names(pgs), value=TRUE)
if (length(pgs_names) != cfg$n_pgs_expected) stop('Expected ',cfg$n_pgs_expected,' PGSs; found ',length(pgs_names))
if (anyDuplicated(c(names(fi),names(pgs)[-1],names(pcs)[-1]))) stop('Non-ID column names overlap across inputs')
long <- fi %>%
  dplyr::inner_join(pgs[,c('id',pgs_names)],by='id') %>%
  dplyr::inner_join(pcs[,c('id',pc_names)],by='id') %>%
  dplyr::filter(!is.na(age),!is.na(frailty),age>=cfg$min_age)
long$sex <- factor(long$sex)
if (any(long$frailty < 0, na.rm=TRUE)) stop('Frailty index must be nonnegative')
if (cfg$sqrt_frailty) long$frailty <- sqrt(long$frailty)
if (cfg$standardise_frailty) {
  # Derive mean/sd from unique participants at baseline rather than all visits.
  baseline_pre <- long %>% arrange(id,age,wave) %>% group_by(id) %>% slice(1) %>% ungroup()
  y_mu <- mean(baseline_pre$frailty); y_sd <- sd(baseline_pre$frailty)
  if (!is.finite(y_sd) || y_sd==0) stop('Outcome has zero variance')
  long$frailty <- (long$frailty-y_mu)/y_sd
}
base <- long %>% arrange(id,age,wave) %>% group_by(id) %>% slice(1) %>% ungroup()
covs <- c('age','sex',pc_names)
features <- c(covs,pgs_names)
# To ensure all analyses use identical baseline participants, use complete cases.
# The original study's missing-data rule should be checked against this choice.
base_cc <- base[complete.cases(base[,c('frailty',features)]),]
if (nrow(base_cc)<cfg$cv_folds+1) stop('Too few complete baseline observations')
dir.create(cfg$output_dir,showWarnings=FALSE,recursive=TRUE)
write.csv(data.frame(n_baseline=nrow(base_cc),n_longitudinal=n_distinct(long$id),
                     n_longitudinal_records=nrow(long)),
          file.path(cfg$output_dir,'sample_sizes.csv'),row.names=FALSE)

# 1. SINGLE-PGS MODELS: incremental R2 over age, sex and ancestry PCs.
null_fit <- lm(reformulate(covs,response='frailty'), data=base_cc)
null_r2 <- summary(null_fit)$r.squared
single <- lapply(pgs_names,function(score){
  mod <- lm(reformulate(c(covs,score),response='frailty'),data=base_cc)
  sm <- summary(mod)
  row <- sm$coefficients[score,,drop=FALSE]
  data.frame(PGS=score,beta=row[1,'Estimate'],se=row[1,'Std. Error'],
             p=row[1,'Pr(>|t|)'],R2_full=sm$r.squared,R2_null=null_r2,
             delta_R2=sm$r.squared-null_r2)
}) %>% dplyr::bind_rows()
single$p_fdr <- p.adjust(single$p,method='BH')
write.csv(single,file.path(cfg$output_dir,'single_pgs_models.csv'),row.names=FALSE)

# 2. MULTI-PGS ELASTIC NET: 10-fold / 10-repeat caret CV; 10 lambdas.
# glmnet standardizes X internally when fitting, including in CV refits.
# Covariates are unpenalised; PGSs are penalised, unlike a penalised-all model.
make_x <- function(data) model.matrix(reformulate(features),data=data)[,-1,drop=FALSE]
X <- make_x(base_cc)
y <- base_cc$frailty
if (any(!is.finite(X))) stop('Non-finite predictors in model matrix')
# Match columns generated by model.matrix (e.g. sex category indicators).
penalty <- as.numeric(colnames(X) %in% pgs_names)
if (sum(penalty)!=length(pgs_names)) stop('PGS names not correctly encoded in model matrix')
# glmnet lambda scale is response-dependent; get sensible candidates from a pilot fit.
pilot <- glmnet::glmnet(X,y,alpha=cfg$alpha,penalty.factor=penalty,standardize=TRUE)
lambda_grid <- exp(seq(log(max(pilot$lambda)),log(min(pilot$lambda)),length.out=cfg$n_lambdas))
set.seed(cfg$seed)
ctrl <- caret::trainControl(method='repeatedcv',number=cfg$cv_folds,
                           repeats=cfg$cv_repeats,savePredictions='final')
fit <- caret::train(x=X,y=y,method='glmnet',metric='RMSE',maximize=FALSE,
                    trControl=ctrl,tuneGrid=expand.grid(alpha=cfg$alpha,lambda=lambda_grid),
                    penalty.factor=penalty,standardize=TRUE)
best_lambda <- fit$bestTune$lambda
write.csv(fit$results,file.path(cfg$output_dir,'elastic_net_cv_results.csv'),row.names=FALSE)
write.csv(data.frame(alpha=cfg$alpha,best_lambda=best_lambda,
                     RMSE_CV=min(fit$results$RMSE)),
          file.path(cfg$output_dir,'elastic_net_best_tuning.csv'),row.names=FALSE)
final <- glmnet::glmnet(X,y,alpha=cfg$alpha,lambda=best_lambda,
                         penalty.factor=penalty,standardize=TRUE)
betas <- as.matrix(coef(final))[,1]
coef_table <- data.frame(predictor=names(betas),beta=as.numeric(betas),
                         nonzero=as.numeric(betas)!=0)
write.csv(coef_table,file.path(cfg$output_dir,'elastic_net_coefficients.csv'),row.names=FALSE)

# 3. BOOTSTRAP COEFFICIENT STABILITY at the CV-selected lambda.
# Each resample refits the model; lambda is held fixed. Not a formal post-selection p-value.
set.seed(cfg$seed)
boot <- matrix(NA_real_,nrow=cfg$n_boot,ncol=length(betas),
               dimnames=list(NULL,names(betas)))
for (b in seq_len(cfg$n_boot)) {
  rows <- sample(seq_len(nrow(X)),replace=TRUE)
  fit_b <- glmnet::glmnet(X[rows,,drop=FALSE],y[rows],alpha=cfg$alpha,
                           lambda=best_lambda,penalty.factor=penalty,standardize=TRUE)
  boot[b,] <- as.matrix(coef(fit_b))[,1]
}
bootstrap_summary <- data.frame(
  predictor=colnames(boot),
  mean_beta=colMeans(boot),
  q025=apply(boot,2,quantile,probs=.025),
  q975=apply(boot,2,quantile,probs=.975),
  selection_frequency=colMeans(boot!=0)
)
write.csv(bootstrap_summary,file.path(cfg$output_dir,'elastic_net_bootstrap.csv'),row.names=FALSE)

# 4. LONGITUDINAL GAMMs: nonlinear age trend and age-by-PGS smooth interaction.
# Fit separately for each score, adjusting for sex and ancestry PCs.
long_cc <- long[complete.cases(long[,c('id','age','frailty',covs,pgs_names)]),]
visit_count <- table(long_cc$id)
long_cc <- long_cc[long_cc$id %in% names(visit_count[visit_count>=cfg$minimum_long_visits]),]
long_cc$id <- factor(long_cc$id)
if (n_distinct(long_cc$id)<10) stop('Too few people with >=3 complete frailty observations')
traj_pgs <- if(is.null(cfg$longitudinal_pgs)) pgs_names else cfg$longitudinal_pgs
if (!all(traj_pgs %in% pgs_names)) stop('Unknown score in longitudinal_pgs')
# Centering score makes the main age smooth correspond to mean genetic liability.
long_results <- list()
for (score in traj_pgs) {
  d <- long_cc
  s <- sd(d[[score]])
  if (!is.finite(s) || s==0) next
  d$score_z <- (d[[score]]-mean(d[[score]]))/s
  # Numeric 'by' smooth models the age-varying effect of PGS; include main effect.
  formula_gam <- as.formula(paste0(
    'frailty ~ s(age, k=',cfg$spline_k,') + score_z + s(age, by=score_z, k=',
    cfg$spline_k,') + ',paste(c('sex',pc_names),collapse=' + ')))
  m <- tryCatch(gamm4::gamm4(formula_gam,random=~(1|id),data=d),
                error=function(e){warning('GAMM failed for ',score,': ',conditionMessage(e));NULL})
  if (is.null(m)) next
  tab <- summary(m$gam)$s.table
  idx <- grep('s\\(age\\):score_z',rownames(tab))
  if(length(idx)!=1) {
    warning('Cannot locate interaction smooth for ',score)
    next
  }
  long_results[[score]] <- data.frame(
    PGS=score,n_participants=n_distinct(d$id),n_observations=nrow(d),
    edf_interaction=tab[idx,'edf'],p_interaction=tab[idx,'p-value'])
}
if(length(long_results)) {
  trajectories <- bind_rows(long_results)
  trajectories$p_fdr <- p.adjust(trajectories$p_interaction,method='BH')
  write.csv(trajectories,file.path(cfg$output_dir,'longitudinal_gamm.csv'),row.names=FALSE)
}
cat('Analysis completed; outputs in ',cfg$output_dir,'\n',sep='')
