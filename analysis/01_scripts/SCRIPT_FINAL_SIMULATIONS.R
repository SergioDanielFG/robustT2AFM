# =====================================================================
# SCRIPT_FINAL_SIMULATIONS.R
# ---------------------------------------------------------------------
# WHAT IT COMPUTES
#   The whole Monte Carlo simulation study, with a single logic and a
#   single seed scheme, producing four tables:
#     - false-alarm calibration: the ARL0 each method's own analytic
#       limit produces under control
#     - detection power: the true positive rate of both methods,
#       equalised to the same ARL0, for four mean shifts and two
#       contamination levels
#     - design sensitivity: the same demanding scenario varying h, J
#       and I one at a time
#     - contamination sensitivity: the same scenario raising the share
#       of contaminated Phase 1 batches to 30% and 40%
#
# DESIGN DECISIONS THAT AFFECT THE VALIDITY OF THE RESULT
#   - The two metrics answer different questions and are computed
#     differently. In the calibration each method keeps its OWN analytic
#     limit, because the question is whether the limit that would be
#     used in practice meets the false-alarm rate it promises. In the
#     detection study both methods are first brought to the same
#     false-alarm rate through the empirical 1-alpha quantile of their
#     own statistic over an independent set of N_CAL in-control
#     batches, so that a difference in TPR reflects the quality of the
#     estimator and not a stricter threshold.
#   - That calibration set is the same for both methods within a
#     replicate and is regenerated in every replicate.
#   - The equicorrelation structure Sigma_jj' = rho, Sigma_jj = 1
#     privileges no particular relation between variables.
#   - The number of contaminated batches is FIXED at 6 in the power and
#     design-sensitivity tables; K is NOT varied, which is a
#     methodological decision: the effect of K on the limit is measured
#     separately in 15_arl0_convergence_K_v3.R.
#   - Each table uses its own seed offsets (SEED_BASE*100 for the
#     calibration, *200/*250/*300 for the power study, *400/*450/*500
#     for the design sensitivity and *800/*850/*900 for the
#     contamination sensitivity), and each replicate seeds its own
#     stream, so every table is reproducible on its own and the ones
#     that share offsets are comparable replicate by replicate.
#
# ANCHORS
#   None. This script is the ORIGIN of the published tables: the
#   anchors used by the other scripts of the folder are the values it
#   produces, not values it checks.
#
# SUCCESS CRITERIA
#   None fixed in advance. This is the descriptive campaign; the
#   pre-registered criteria belong to the scripts that test specific
#   claims against these tables.
#
# OUTPUT
#   02_resultados/table_3_2_calibration.csv
#   02_resultados/table_3_3_power.csv
#   02_resultados/table_3_4_sensitivity.csv
#   02_resultados/table_3_4b_contamination.csv
# =====================================================================
library(robustT2AFM); library(MASS)
dir.create("02_resultados", showWarnings = FALSE)

# --- Global parameters (set once, for all the tables) ---
N_REP <- 2000
N_CAL <- 5000
K1<-30; K2<-100; I<-20; J<-4
RHO<-0.6; H_MCD<-0.67; ALPHA<-0.001
OB<-6; OR<-0.20; OS<-4
SEED_BASE<-2026
VARS<-paste0("Var",1:J)

make_equi <- function(p,rho){ S<-matrix(rho,p,p); diag(S)<-1; S }
Sigma_EQ  <- make_equi(J,RHO)
t2_manual <- function(Xb,mu,Si,n){ d<-colMeans(Xb)-mu; as.numeric(n*t(d)%*%Si%*%d) }
se <- function(x) sd(x)/sqrt(length(x))
t0 <- Sys.time()

# =====================================================================
# CALIBRATION OF THE ARL0 (each method with its own analytic UCL)
# =====================================================================
cat("== Table: false-alarm calibration ==\n")
far <- matrix(NA,N_REP,2,dimnames=list(NULL,c("Rob","Hot")))
for(rep in seq_len(N_REP)){
  set.seed(SEED_BASE*100+rep)
  sim<-simulate_batch_process(K1=K1,K2=K2,I=I,J=J,rho=RHO,Sigma=Sigma_EQ,
                              outlier_batches_F1=OB,outlier_rate=OR,outlier_shift=OS,prop_ooc_F2=0)
  f1<-subset(sim,Phase=="Phase 1"); f2<-subset(sim,Phase=="Phase 2")
  calR<-suppressMessages(calibrate_afm_mcd(f1,VARS,mcd_alpha=H_MCD))
  calH<-hotelling_classical_calibrate(f1,VARS)
  uR<-ucl_F_adjusted(calR,I=I,alpha=ALPHA)$UCL
  uH<-hotelling_classical_ucl(K=calH$n_batches,I=I,J=J,alpha=ALPHA,phase="II")$UCL
  far[rep,"Rob"]<-mean(monitor_afm_mcd(f2,calR,VARS)$T2 > uR)
  far[rep,"Hot"]<-mean(hotelling_classical_monitor(f2,calH,VARS)$T2 > uH)
  if(rep%%500==0) cat("  rep",rep,"\n")
}
fm<-colMeans(far); sf<-apply(far,2,se)
tab32<-data.frame(Limit=c("Robusto F (m*)","Clasico Montgomery"),
                  FAR=round(fm,6), SE_FAR=signif(sf,3), ARL0=round(1/fm,1))
print(tab32,row.names=FALSE)
write.csv(tab32,"02_resultados/table_3_2_calibration.csv",row.names=FALSE)

# =====================================================================
# DETECTION POWER (equalised by empirical quantile)
# =====================================================================
cat("== Table: detection power ==\n")
res33<-data.frame()
for(sh in c(0.5,1.0,1.5,2.0)) for(ob in c(0,6)){
  tprH<-numeric(N_REP); tprR<-numeric(N_REP)
  for(rep in seq_len(N_REP)){
    set.seed(SEED_BASE*200+rep)
    sim<-simulate_batch_process(K1=K1,K2=1,I=I,J=J,rho=RHO,Sigma=Sigma_EQ,
                                outlier_batches_F1=ob,outlier_rate=OR,outlier_shift=OS,prop_ooc_F2=0)
    f1<-subset(sim,Phase=="Phase 1")
    calR<-suppressMessages(calibrate_afm_mcd(f1,VARS,mcd_alpha=H_MCD))
    calH<-hotelling_classical_calibrate(f1,VARS)
    SiR<-solve(calR$Sw); muR<-calR$mu_r; SiH<-solve(calH$Sp); muH<-calH$mu_global
    set.seed(SEED_BASE*250+rep)
    f1c<-lapply(seq_len(N_CAL),function(k){X<-mvrnorm(I,rep(0,J),Sigma_EQ);colnames(X)<-VARS;X})
    qR<-quantile(sapply(f1c,t2_manual,muR,SiR,I),1-ALPHA)
    qH<-quantile(sapply(f1c,t2_manual,muH,SiH,I),1-ALPHA)
    set.seed(SEED_BASE*300+rep)
    f2<-lapply(seq_len(K2),function(k){X<-mvrnorm(I,rep(sh,J),Sigma_EQ);colnames(X)<-VARS;X})
    tprR[rep]<-mean(sapply(f2,t2_manual,muR,SiR,I)>qR)
    tprH[rep]<-mean(sapply(f2,t2_manual,muH,SiH,I)>qH)
  }
  res33<-rbind(res33,data.frame(contam_batches=ob,shift_F2=sh,
                                TPR_Hot=round(mean(tprH),4),SE_Hot=round(se(tprH),4),
                                TPR_Rob=round(mean(tprR),4),SE_Rob=round(se(tprR),4),
                                Advantage=round(mean(tprR)-mean(tprH),4)))
  cat(sprintf("  contam=%d shift=%.1f Adv=%+.3f\n",ob,sh,mean(tprR)-mean(tprH)))
}
print(res33,row.names=FALSE)
write.csv(res33,"02_resultados/table_3_3_power.csv",row.names=FALSE)

# =====================================================================
# SENSITIVITY to h, J and I (6 contaminated batches FIXED, K not varied)
# =====================================================================
cat("== Table: design sensitivity h/J/I ==\n")
CFG<-data.frame(config=c("base","h=0.75","J=2","J=6","I=15"),
                K=c(30,30,30,30,30), I=c(20,20,20,20,15), J=c(4,4,2,6,4),
                h=c(0.67,0.75,0.67,0.67,0.67), stringsAsFactors=FALSE)
res34<-data.frame()
for(r in seq_len(nrow(CFG))){
  Kc<-CFG$K[r]; Ic<-CFG$I[r]; Jc<-CFG$J[r]; hc<-CFG$h[r]
  vv<-paste0("Var",1:Jc); Sg<-make_equi(Jc,RHO)
  tprH<-numeric(N_REP); tprR<-numeric(N_REP)
  for(rep in seq_len(N_REP)){
    set.seed(SEED_BASE*400+rep)
    sim<-simulate_batch_process(K1=Kc,K2=1,I=Ic,J=Jc,rho=RHO,Sigma=Sg,
                                outlier_batches_F1=6,outlier_rate=OR,outlier_shift=OS,prop_ooc_F2=0)
    f1<-subset(sim,Phase=="Phase 1")
    calR<-suppressMessages(calibrate_afm_mcd(f1,vv,mcd_alpha=hc))
    calH<-hotelling_classical_calibrate(f1,vv)
    SiR<-solve(calR$Sw); muR<-calR$mu_r; SiH<-solve(calH$Sp); muH<-calH$mu_global
    set.seed(SEED_BASE*450+rep)
    f1c<-lapply(seq_len(N_CAL),function(k){X<-mvrnorm(Ic,rep(0,Jc),Sg);colnames(X)<-vv;X})
    qR<-quantile(sapply(f1c,t2_manual,muR,SiR,Ic),1-ALPHA)
    qH<-quantile(sapply(f1c,t2_manual,muH,SiH,Ic),1-ALPHA)
    set.seed(SEED_BASE*500+rep)
    f2<-lapply(seq_len(K2),function(k){X<-mvrnorm(Ic,rep(1.0,Jc),Sg);colnames(X)<-vv;X})
    tprR[rep]<-mean(sapply(f2,t2_manual,muR,SiR,Ic)>qR)
    tprH[rep]<-mean(sapply(f2,t2_manual,muH,SiH,Ic)>qH)
  }
  res34<-rbind(res34,data.frame(config=CFG$config[r],K=Kc,I=Ic,J=Jc,h=hc,
                                TPR_Hot=round(mean(tprH),4),SE_Hot=round(se(tprH),4),
                                TPR_Rob=round(mean(tprR),4),SE_Rob=round(se(tprR),4),
                                Advantage=round(mean(tprR)-mean(tprH),4)))
  cat(sprintf("  %-7s Adv=%+.3f\n",CFG$config[r],mean(tprR)-mean(tprH)))
}
print(res34,row.names=FALSE)
write.csv(res34,"02_resultados/table_3_4_sensitivity.csv",row.names=FALSE)

# =====================================================================
# SENSITIVITY to the contamination level (30%, 40%)
# =====================================================================
cat("== Table: contamination sensitivity 30/40 ==\n")
resB<-data.frame()
for(ob in c(9,12)){
  tprH<-numeric(N_REP); tprR<-numeric(N_REP)
  for(rep in seq_len(N_REP)){
    set.seed(SEED_BASE*800+rep)
    sim<-simulate_batch_process(K1=30,K2=1,I=20,J=4,rho=RHO,Sigma=Sigma_EQ,
                                outlier_batches_F1=ob,outlier_rate=OR,outlier_shift=OS,prop_ooc_F2=0)
    f1<-subset(sim,Phase=="Phase 1")
    calR<-suppressMessages(calibrate_afm_mcd(f1,VARS,mcd_alpha=H_MCD))
    calH<-hotelling_classical_calibrate(f1,VARS)
    SiR<-solve(calR$Sw); muR<-calR$mu_r; SiH<-solve(calH$Sp); muH<-calH$mu_global
    set.seed(SEED_BASE*850+rep)
    f1c<-lapply(seq_len(N_CAL),function(k){X<-mvrnorm(20,rep(0,4),Sigma_EQ);colnames(X)<-VARS;X})
    qR<-quantile(sapply(f1c,t2_manual,muR,SiR,20),1-ALPHA)
    qH<-quantile(sapply(f1c,t2_manual,muH,SiH,20),1-ALPHA)
    set.seed(SEED_BASE*900+rep)
    f2<-lapply(seq_len(K2),function(k){X<-mvrnorm(20,rep(1.0,4),Sigma_EQ);colnames(X)<-VARS;X})
    tprR[rep]<-mean(sapply(f2,t2_manual,muR,SiR,20)>qR)
    tprH[rep]<-mean(sapply(f2,t2_manual,muH,SiH,20)>qH)
  }
  resB<-rbind(resB,data.frame(config=paste0("contam=",round(100*ob/30),"%"),
                              contam_pct=round(100*ob/30),
                              TPR_Hot=round(mean(tprH),4),SE_Hot=round(se(tprH),4),
                              TPR_Rob=round(mean(tprR),4),SE_Rob=round(se(tprR),4),
                              Advantage=round(mean(tprR)-mean(tprH),4)))
  cat(sprintf("  contam=%d%% Adv=%+.3f\n",round(100*ob/30),mean(tprR)-mean(tprH)))
}
print(resB,row.names=FALSE)
write.csv(resB,"02_resultados/table_3_4b_contamination.csv",row.names=FALSE)

cat("\n===== CAMPAIGN COMPLETE |",round(difftime(Sys.time(),t0,units="mins"),1),"min =====\n")
cat("Check 02_resultados/ for the four CSV files.\n")
