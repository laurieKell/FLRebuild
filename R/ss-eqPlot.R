I #' Build equilibrium and production-curve inputs from SS output
#'
#' Process Stock Synthesis output to calculate yield and surplus production
#' series for equilibrium-curve plotting.
#'
#' @param object A list returned by `r4ss::SS_output()`, or a list of character
#'   paths to Stock Synthesis run directories.
#' @param maxY Multiplier used to set the upper y-limit proxy for triangle data.
#'
#' @return A list with elements `tseries`, `curve`, `refpts`, `triangle`, and
#'   `derived`. For list input containing paths, each data.frame includes an
#'   identifying column `id`.
#' @export
setMethod("curveSS", signature(object="list"), function(object, maxY = 1.5) {
  # If input is a list of character paths, dispatch per element and bind with ID.
  if (length(object) > 0 && all(vapply(object, function(x) is.character(x) && length(x) == 1, logical(1)))) {
    ids <- names(object)
    if (is.null(ids) || any(ids == "")) ids <- as.character(seq_along(object))

    runs <- lapply(seq_along(object), function(i) {
      out <- curveSS(object[[i]], maxY = maxY)
      lapply(out, function(df) {
        df <- as.data.frame(df)
        tryIt(cbind(id = ids[i], df, stringsAsFactors = FALSE))
      })
    })

    parts <- names(runs[[1]])
    combined <- lapply(parts, function(p) do.call(rbind, lapply(runs, `[[`, p)))
    names(combined) <- parts
    return(combined)
  }

  eqlYield  =object[["equil_yield"]]
  timeseries=object[["timeseries"]]
  rfs=c("SSB_unfished","Totbio_unfished","SmryBio_unfished","Recr_unfished","SSB_Btgt","SPR_Btgt","annF_Btgt",  
        "Dead_Catch_Btgt","SSB_SPR", "annF_SPR","Dead_Catch_SPR","SSB_MSY","SPR_MSY","annF_MSY",   
        "Dead_Catch_MSY", "Ret_Catch_MSY","B_MSY/SSB_unfished")
  dq=object[["derived_quants"]][unique(object$derived_quants$Label)[rfs],]
  
  if (!("Tot_Catch"%in%names(eqlYield))&("Catch"%in%names(eqlYield)))
    names(eqlYield)[names(eqlYield)=="Catch"]="Tot_Catch"

  if (any(tolower(names(timeseries))%in%"yr"))
    names(timeseries)[tolower(names(timeseries))=="yr"]="year"

  # Filter and prepare data for plotting
  ts <- timeseries %>%
    dplyr::filter(!Era %in% c("VIRG", "FORE")) %>%
    dplyr::mutate(catch_total = rowSums(dplyr::select(., dplyr::starts_with("dead(B)")), na.rm = TRUE)) %>%
    dplyr::group_by(year, Seas) %>%
    dplyr::summarise(
      sum_bio_all = sum(Bio_all, na.rm = TRUE),
      sum_spawn_bio = sum(SpawnBio, na.rm = TRUE),
      sum_catch_total = sum(catch_total, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      mean_bio_all = mean(sum_bio_all),
      mean_spawn_bio = mean(sum_spawn_bio),
      catch_total = sum(sum_catch_total),
      .groups = "drop"
    )
  names(ts)=c("year","biomass","ssb","yield")
  
  nyears=nrow(ts)
  
  ts$sp   =c(NA, ts$biomass[-1] - ts$biomass[-nyears] + ts$yield[-nyears])
  eql     =dplyr::transmute(eqlYield, yield=Tot_Catch,ssb=SSB)
  rfs     =dplyr::transmute(subset(eqlYield, Tot_Catch==max(Tot_Catch)), msy=Tot_Catch,bmsy=SSB)[1,]
  maxY    =signif(max(c(ts$yield,eql$yield))*maxY,1)
  triangle=data.frame(x=c(rfs$bmsy, rfs$bmsy, rfs$bmsy*maxY/rfs$msy, rfs$bmsy),
                      y=c(rfs$msy,      maxY,                  maxY, rfs$msy))
  vBiomass=vBio(object)
  
  prodFun<-splinefun(x=eql$ssb, y =eql$yield, method = "natural")
  
  ts=tryIt(cbind(ts,pf=prodFun(ts$ssb)))
  #ts=tryIt(cbind((ts,pe=c(with(ts, (ssb[-1]-tail(ssb,-1)-tail(yield,-1)+tail(pf,-1))/tail(ssb,-1)),NA))
  ts=tryIt(cbind(ts,pe =c(with(ts, log(ssb[-1]/tail(ssb-yield+pf,-1))),NA),
                    pe2=c(with(ts,     (tail(ssb-yield+pf,-1)-ssb[-1])/tail(ssb-yield+pf,-1)), NA)))
  
  labs =object$parameters[grep("Rec",object$parameters$Label),2]
  years=as.integer(sub(".*_(\\d{4})$", "\\1", labs))
  ts   =merge(ts,data.frame(year=years,
                      recDevs=object$parameters[grep("Rec",object$parameters$Label),3]))

  return(list(tseries=merge(ts,vBiomass),curve=eql,refpts=rfs,triangle=triangle,derived=dq))
  })

#' Build equilibrium and production-curve inputs from SS output directory
#'
#' Read Stock Synthesis output from a directory and return the same structure as
#' `curveSS` for a list object.
#'
#' @param object Character path (or vector of paths) to Stock Synthesis run
#'   directories.
#' @param maxY Multiplier used to set the upper y-limit proxy for triangle data.
#' @param covar Logical; passed to `r4ss::SS_output`.
#'
#' @return A list with elements `tseries`, `curve`, `refpts`, `triangle`, and
#'   `derived`.
#' @export
setMethod("curveSS", signature(object = "character"), function(object, maxY = 1.5, covar = FALSE) {
  if (length(object) > 1) {
    ids <- names(object)
    if (is.null(ids) || any(ids == "")) ids <- as.character(seq_along(object))

    runs <- lapply(seq_along(object), function(i) {
      out <- curveSS(object[[i]], maxY = maxY, covar = covar)
      lapply(out, function(df) {
        df <- as.data.frame(df)
        tryIt(cbind(id = object[i], df, stringsAsFactors = FALSE))
      })
    })

    parts <- names(runs[[1]])
    combined <- lapply(parts, function(p) do.call(rbind, lapply(runs, `[[`, p)))
    names(combined) <- parts
    return(combined)
  }

  ss_rep <- r4ss::SS_output(dir = object, covar = covar, verbose = FALSE, printstats = FALSE)
  curveSS(ss_rep, maxY = maxY)
})



logistic <- plogis  # maps to standard logistic function
logit <- qlogis     # maps to inverse logistic function

#' Estimate vulnerable biomass from SS report components
#'
#' Computes sex-specific vulnerable biomass using selectivity parameters and
#' summary biomass time series from Stock Synthesis output.
#'
#' @param object A list returned by `r4ss::SS_output()`.
#'
#' @return A data frame with columns `year`, `female`, and `male`.
#' @export
setMethod("vBio", signature(object = "list"), function(object) {
  # Get selectivity parameters
  sel_F_peak  =unlist(c(subset(object$parameters, Label=="SzSel_Fem_Peak_Fishery_2(2)","Value")))
  sel_F_ascend=unlist(c(subset(object$parameters, Label=="SzSel_Fem_Ascend_Fishery_2(2)","Value")))
  sel_M_peak  =unlist(c(subset(object$parameters, Label=="SzSel_Male_Peak_Fishery_3(3)","Value")))
  sel_M_ascend=unlist(c(subset(object$parameters, Label=="SzSel_Male_Ascend_Fishery_3(3)","Value")))
  
  # Calculate selectivity using logistic transform
  calc_sel<-function(len, peak, ascend) 
    logistic((len - peak)/exp(ascend))
  
  # Get length bins from data
  len_bins=as.numeric(gsub("\\D", "", names(object$sizeselex)[-(1:4)]))
  len_bins[is.na(len_bins)]=0
  
  # Calculate sex-specific selectivity
  sel_F=calc_sel(len_bins, sel_F_peak, sel_F_ascend)
  sel_M=calc_sel(len_bins, sel_M_peak, sel_M_ascend)
  
  # Extract and scale biomass (assuming SmryBio is in log space)
  bio_F=object$timeseries$`SmryBio_SX:1_GP:1`
  bio_M=object$timeseries$`SmryBio_SX:2_GP:1`
  
  # Calculate vulnerable biomass
  rtn=data.frame(
    year   = object$timeseries$Yr,
    female = bio_F * sum(sel_F),
    male   =NA)
  
  if(!is.null(bio_M)&!is.null(sel_M))
    rtn$male = bio_M * sum(sel_M)

  return(rtn)
})


if (FALSE){
    library(remotes)
    remotes::install_github("flr/FLRebuild")
    
    
    library(FLRebuild)
    library(r4ss)
    library(tidyverse)
    
    pfs=curveSS("C:/active/flrpapers/PE/SS3/alb-natl")
      #c("C:/active/flrpapers/PE/SS3/sma/Scenario-3"))
    
    p1=ggplot(pfs$tseries)+
      geom_line(aes(year,log(yield/pf)))+
      geom_point(aes(year,log(yield/pf)))+
      geom_point(aes(year,log(sp/pf)),col="red")+
      geom_line( aes(year,log(sp/pf)),col="red")+
      geom_line(aes( year,recDevs),col="blue")+
      geom_point(aes(year,recDevs),col="blue")+
      geom_hline(aes(yintercept=0),linetype=2)+
      theme_minimal()+
      xlab("Year")+ylab("Process Error")+
      labs(title="Process Error")

    p2=ggplot(pfs$curve)+
      geom_line( aes(ssb,yield))+
      geom_point(aes(ssb,yield),fill="grey",shape=21,data=pfs$tseries)+
      geom_path(aes(ssb,yield), lwd=0.25,data=pfs$tseries)+
      theme_minimal()+
      xlab("SSB")+ylab("Yield")+
      labs(title="Yield v SSB")+
      scale_y_continuous(limits=c(NA,9100))
    
    p3=ggplot(pfs$curve)+
      geom_line( aes(ssb,yield))+
      geom_point(aes(ssb,sp),fill="red",shape=21,data=pfs$tseries)+
      geom_path(aes(ssb,sp), col="red",lwd=0.25,data=pfs$tseries)+
      theme_minimal()+
      xlab("SSB")+ylab("Surplus Production")+
      labs(title="Surplus Production v SSB")+
      scale_y_continuous(limits=c(NA,9100))
   
    (p2 | p3) / p1}



ss3Objects<-function(replist,fls,seas=1) {
  
  ## natage → N
  nat=replist$natage
  if (is.null(nat) || !nrow(nat)) stop("replist$natage not found")
  
  if (!"Beg/Mid" %in% names(nat)) stop("'Beg/Mid' column not found in natage")
  nat=nat[nat$Seas == seas & nat[["Beg/Mid"]] == "B" & nat$Era == "TIME",
          , drop = FALSE]
  
  age_cols=as.character(0:100)
  age_cols=age_cols[age_cols %in% names(nat)]
  if (!length(age_cols)) stop("No age columns found in natage")
  
  ages =as.numeric(age_cols)
  years=sort(unique(nat$Yr))
  
  N=matrix(NA_real_, nrow = length(ages), ncol = length(years),
           dimnames = list(age = ages, year = years))
  for (j in seq_along(years)) {
    yr=years[j]
    rows=nat[nat$Yr == yr, age_cols, drop = FALSE]
    N[, j]=colSums(rows, na.rm = TRUE)
  }
  
  ## wtatage → W
  wt=replist$wtatage
  if (is.null(wt) || !nrow(wt)) stop("replist$wtatage not found")
  wt=wt[wt$seas == seas, , drop = FALSE]
  
  W=matrix(NA_real_, nrow = length(ages), ncol = length(years),
           dimnames = list(age = ages, year = years))
  for (j in seq_along(years)) {
    yr=years[j]
    rows=wt[wt$year == yr, age_cols, drop = FALSE]
    if (!nrow(rows)) {
      W[, j]=NA_real_
    } else {
      W[, j]=colMeans(rows, na.rm = TRUE)
    }
  }
  
  ## timeseries → SSB
  ts=replist$timeseries
  if (is.null(ts) || !nrow(ts)) stop("replist$timeseries not found")
  ts_sub=ts[ts$Era == "TIME" & ts$Seas == seas,
            c("Yr", "SpawnBio"), drop = FALSE]
  names(ts_sub)=c("year", "SSB")
  
  ## m(fls) → M
  M_flq=m(fls)  # FLQuant: age × year × unit × season × area × iter
  M_slice=M_flq[as.character(ages),
                as.character(years),
                1, 1, 1, 1, drop = TRUE]
  M_mat=as.matrix(M_slice)
  
  list(
    N     = N,
    M     = M_mat,
    W     = W,
    SSB   = ts_sub,
    ages  = ages,
    years = years)}

ss3qSSB<-function(W, Mat, M, theta = 0.5, plusgroup = TRUE) {
  ages =seq_len(nrow(W))
  years=seq_len(ncol(W))
  
  q=matrix(0, nrow = nrow(W), ncol = ncol(W), dimnames = dimnames(W))
  
  for (ia in ages) {
    for (jt in years) {
      q_at=0
      surv=1
      
      for (ik in ia:nrow(W)) {
        ku=jt + (ik - ia)
        if (ku > ncol(W)) break
        
        if (ik > ia) {
          surv=surv * exp(-M[ik - 1, ku - 1])
        }
        
        q_at=q_at + surv * W[ik, ku] * Mat[ik, ku] * exp(-theta * M[ik, ku])
      }
      
      if (plusgroup) {
        last_u=jt + (nrow(W) - ia)
        if (last_u <= ncol(W)) {
          q_at=q_at + surv *
            (W[nrow(W), last_u] *
               Mat[nrow(W), last_u] *
               exp(-theta * M[nrow(W), last_u]) /
               (1 - exp(-M[nrow(W), last_u])))
        }
      }
      
      q[ia, jt]=q_at}}
  
  q}

ss3qStationary<-function(W, Mat, M, theta = 0.5, plusgroup = TRUE) {
  Wa  =rowMeans(W,   na.rm = TRUE)
  Mata=rowMeans(Mat, na.rm = TRUE)
  Ma  =rowMeans(M,   na.rm = TRUE)
  
  A=length(Wa)
  q_age=numeric(A)
  
  for (a in seq_len(A)) {
    surv=1
    q_at=0
    
    for (k in a:A) {
      if (k > a) {
        surv=surv * exp(-Ma[k - 1])
      }
      
      q_at=q_at + surv * Wa[k] * Mata[k] * exp(-theta * Ma[k])
    }
    
    if (plusgroup) {
      q_at=q_at + surv *
        (Wa[A] * Mata[A] * exp(-theta * Ma[A]) / (1 - exp(-Ma[A])))
    }
    
    q_age[a]=q_at}
  
  q_mat=matrix(
    rep(q_age, times = ncol(W)),
    nrow = A,
    ncol = ncol(W),
    dimnames = dimnames(W))
  
  q_mat}

ss3EqvlCatch<-function(CN, W, Mat, M,
                       theta = 0.5,
                       plusgroup = TRUE,
                       method = c("finite", "stationary")) {
  method=match.arg(method)
  
  stopifnot(
    is.matrix(CN), is.matrix(W), is.matrix(Mat), is.matrix(M),
    all(dim(CN)  == dim(W)),
    all(dim(Mat) == dim(W)),
    all(dim(M)   == dim(W)))
  
  q=switch(
    method,
    finite = ss3qSSB(
      W = W, Mat = Mat, M = M,
      theta = theta, plusgroup = plusgroup),
    stationary = ss3qStationary(
      W = W, Mat = Mat, M = M,
      theta = theta, plusgroup = plusgroup))
  
  C_ssb_at=CN * q
  C_ssb=colSums(C_ssb_at, na.rm = TRUE)
  
  list(
    C_ssb = C_ssb,
    C_ssb_at = C_ssb_at,
    q = q,
    method = method)}

ss3VB<-function(N, W, V, Z = NULL, timing = c("start", "mid"), fraction = 0.5) {
  timing=match.arg(timing)
  if (timing == "start") {
    return(colSums(N * W * V, na.rm = TRUE))}
  
  if (is.null(Z)) stop("Z required for mid-year vulnerable biomass")
  colSums(N * exp(-fraction * Z) * W * V, na.rm = TRUE)}

getSP<-function(object){
  
  replist=SS_output(object,verbose=FALSE,covar=FALSE,printstats=FALSE)
  fls    =ss3om:::readFLSss3(object)
  obj    =ss3Objects(replist, fls)
  
  ages =obj$ages
  years=obj$years
  
  ## 1. Maturity, selectivity, catch-at-age from FLR on SS3 grid
  
  Mat=mat(fls)
  Mat=Mat[as.character(ages),as.character(years),1, 1, 1, 1, drop = TRUE]
  Mat=as.matrix(Mat)
  
  # Selectivity proxy: harvest as F-at-age pattern
  Sel=harvest(fls)
  Sel=Sel[as.character(ages),as.character(years),1, 1, 1, 1, drop = TRUE]
  
  V=as.matrix(Sel)
  
  # Catch numbers-at-age
  Ctc=catch.n(fls)  # FLQuant age × year
  Ctc=Ctc[as.character(ages),as.character(years),1, 1, 1, 1, drop = TRUE]
  Ctc=as.matrix(Ctc)
  
  
  ## 2. SSB-equivalent catch and vulnerable biomass
  # Choose q method: "finite" or "stationary"
  q_method="stationary"
  
  ssbEqv=ss3EqvlCatch(CN       = Ctc,
                      W        = obj$W,
                      Mat      = Mat,
                      M        = obj$M,
                      theta    = 0.5,
                      plusgroup= TRUE,
                      method   = q_method)
  
  VB=ss3VB(N      = obj$N,
           W      = obj$W,
           V      = V,
           Z      = NULL,
           timing = "start")
  
  
  
  B_df =obj$SSB                # columns: year, SSB
  C_ssb=ssbEqv$C_ssb           # named vector or same order
  
  yrs=intersect(B_df$year, as.numeric(names(C_ssb)))
  B  =B_df$SSB[match(yrs, B_df$year)]
  C_t=C_ssb[match(yrs, as.numeric(names(C_ssb)))]
  
  # P_t = B_{t+1} - B_t + C_t
  B_t  =B[-length(B)]
  B_tp1=B[-1]
  C_t_ =C_t[-length(C_t)]
  
  P_ssb=B_tp1 - B_t + C_t_
  
  data.frame(
    year = yrs[-length(yrs)],
    ssb  = B_t,
    vb   = VB[-length(VB)],
    sp2  = P_ssb)}


# The plot shows that the stock’s realised yield and surplus production at a given SSB vary widely over time, and that the deviations are strongly structured, not simple white‑noise process error. 
#  
# ## Yield and surplus production vs SSB
#  
# + The yield v SSB panel shows the time series looping around the equilibrium curve, with periods where yield at a given SSB is much higher or lower than expected, indicating transient effects, e.g. due to recruitment or exploitation history. 
# + The surplus production v SSB panel shows a similar looping and vertical spread, showing that productivity at a given SSB changes substantially through time, inconsistent with a constant productivity surplus production model.
#     
# ## Time series of process error
#  
# - The lower panel shows multi‑year runs of positive and negative process error revealing low‑frequency structure or regimes rather than independent deviations around zero. 
#     
# These patterns means estimates of MSY, and stock status are sensitive to the time window used, and management procedures that assume stationary productivity and IID process error are likely over‑confident and potentially biased. 
# 
# In Stock Synthesis the blue series (recruitment deviates) is the primary source of the low‑frequency structure seen in the black and red “process error” series, and strong correlations between them indicate that much of the apparent process error is actually recruitment‑driven rather than generic productivity noise. [sedarweb](https://sedarweb.org/documents/sedar-39-rd08-appendix-a-technical-description-of-the-stock-synthesis-assessment-program/)
# 
# ## Role of recruitment deviations in SS
# 
# - Recruitment deviations in SS are annual log‑scale random effects around the stock–recruit curve that allow recruitment to vary over time while keeping the long‑term mean consistent with the specified relationship via bias‑adjustment. 
# - Blocks of positive rec devs imply extended periods of above‑average recruitment at given SSB, while negative blocks imply the opposite; these directly alter biomass trajectories and hence yield and surplus production at a given SSB. 
# 
# ## Correlation with process error
# 
# - In the plot the blue rec devs track many of the multi‑year swings in the black (log(yield/pf)) and red (log(sp/pf)) series, implying strong correlation between recruitment deviations and the surplus‑production “process error.” 
# - When such correlation exists, process error is not independent white noise around a fixed production curve; instead, recruitment deviations act as a compensatory, systematic driver of productivity, so models or HCRs that assume stationary productivity and IID process error will tend to underestimate uncertainty and can be biased, especially under extreme productivity scenarios. 
# 
# Yes—because SS3 represents most process variability through recruitment deviations, they can easily soak up model misspecification, and that has important consequences for inference and management advice. [iccat](https://www.iccat.int/Documents/CVSP/CV070_2014/n_5/CV070052069.pdf)
# 
# ## How rec devs absorb misspecification
# 
# - In SS3, recruitment deviations are flexible time‑varying random effects, while many other processes (growth, maturity, natural mortality, stock–recruit form, often selectivity) are fixed or strongly constrained, so unexplained variation in the data tends to be pushed into the rec devs.
# 
# - When key biological or fishery processes are mis‑specified, recruitment deviations can develop trends, autocorrelation, or “boom–bust” patterns that compensate for errors in productivity, selectivity, catch history, or data weighting rather than reflecting true recruitment variability. 
# 
# 
# ## Consequences for assessment outputs
# 
# - If rec devs are absorbing misspecification, the implied productivity and reference points (MSY, \(B_{\text{MSY}}\), \(F_{\text{MSY}}\)) can be badly biased, because the model attributes persistent surplus or deficit of biomass to recruitment instead of to wrong life‑history or selectivity assumptions. 
# 
# - This can produce extreme high or low productivity scenarios that still fit the data but diverge strongly from deterministic or alternative models, leading to conflicting advice across model configurations.
# 
# ## Consequences for diagnostics and management
# 
# - Non‑random patterns in rec devs (trends, strong autocorrelation, long runs of one sign) are therefore a diagnostic red flag that the model is using recruitment as a **sink** for unmodelled time‑varying processes and latent misspecification. 
# 
# - If such patterns are ignored and a single SS3 run is used for advice, uncertainty is understated and HCRs tuned to that model can be over‑optimistic or overly pessimistic, so ensemble modelling, alternative structural hypotheses, and explicit process‑error diagnostics on rec devs are needed before relying on the assessment for management.
# 
# 
#     
#     baSS=SS_output("C:/active/flrpapers/PE/SS3/bass-north")
#     
#     dim(baSS$parameters[grep("Rec",baSS$parameters$Label),])
#     
#     
#     bass$tseries=transform(bass$tseries, pe2=sp-pf)
#     bass$tseries=transform(bass$tseries, 
#                            recDevs=baSS$parameters[grep("Rec",baSS$parameters$Label),3])
#     ggplot(bass$tseries)+geom_point(aes(pe,recDevs))
#     
#     ccf(head(bass$tseries$pe,-1),head(bass$tseries$recDevs,-1))
# 
#     
#     sma=curveSS(c("C:/active/flrpapers/PE/SS3/sma/Scenario-1",
#                   "C:/active/flrpapers/PE/SS3/sma/Scenario-2",
#                   "C:/active/flrpapers/PE/SS3/sma/Scenario-3"))
#     
#     ggplot(sma$ts)+geom_line(aes(year,pe,col=id))
#     
#     ggplot(sma$curve)+
#       geom_line(aes(ssb,yield,col=id))+
#       geom_line(aes(ssb,yield,col=id),data=sma$tseries)
#     
#     her=curveSS(c("C:/active/flrpapers/PE/SS3/sma/herring"))
#     
#     ggplot(her$tseries)+
#       geom_line(aes(year,pe))
#     
#     
#     ggplot(her$curve)+
#       geom_line(aes(ssb,yield,col=.id))+
#       geom_path(aes(ssb,yield,col=.id),data=her$tseries)
# }
# 
# 
# # Ignoring PE diagnostics increases the risk that advice is based on a structurally wrong model with over‑confident uncertainty, which can drive biased catch recommendations and delayed corrective action. [edepot.wur](https://edepot.wur.nl/556274)
# # 
# # ## Hidden mis‑specification and biased advice
# # 
# # - Without checking PE for trends, autocorrelation, or regime shifts, time‑varying productivity, selectivity, or catch reporting changes can be misinterpreted as stable dynamics. [repository.library.noaa](https://repository.library.noaa.gov/view/noaa/52235/noaa_52235_DS1.pdf)
# # - This can bias reference points (e.g. \(F_{\text{MSY}}\), \(B_{\text{MSY}}\)) and stock‑status estimates, so catches are set too high when productivity has declined or too low when it has increased, undermining risk equivalence across stocks. [onlinelibrary.wiley](https://onlinelibrary.wiley.com/doi/10.1111/faf.12550)
# # 
# # ## Over‑confident uncertainty and poor prediction skill
# # 
# # - Standard likelihoods assume independent, homoscedastic residuals; if PE has low‑frequency structure, variance and autocorrelation are underestimated, so confidence intervals around biomass and forecasts are too narrow. [spo.nmfs.noaa](https://spo.nmfs.noaa.gov/sites/default/files/TMSPO240.pdf)
# # - Hindcasting studies show that models with unexamined structured PE often have poor out‑of‑sample prediction skill, leading to optimistic rebuilding trajectories and underestimation of the probability of falling below limit reference points. [sciencedirect](https://www.sciencedirect.com/science/article/abs/pii/S0165783616301540)
# # 
# # ## Consequences for management outcomes
# # 
# # - **Higher probability of overfishing and limit breaches**: Biased forecasts combined with narrow uncertainty bands reduce the apparent risk of exceeding \(F_{\text{MSY}}\) or dropping below \(B_{\text{lim}}\), encouraging TACs that are less precautionary than intended. [edepot.wur](https://edepot.wur.nl/556274)
# # - **Slow detection of problems**: If retrospective patterns or PE signals are ignored, managers may only react after strong declines are evident in surveys or catches, by which time rebuilding is longer and more costly. [academic.oup](https://academic.oup.com/icesjms/article-pdf/82/2/fsaf014/61741528/fsaf014.pdf)
# # - **Inconsistent treatment of data‑limited stocks**: For Category‑2 surplus‑production assessments, failing to use PE diagnostics can make these look more reliable than they are, weakening the precautionary buffer that WKLIFE and ICES seek to maintain for data‑limited advice. [ppl-ai-file-upload.s3.amazonaws](https://ppl-ai-file-upload.s3.amazonaws.com/web/direct-files/attachments/11254397/3b640e6f-75f7-479e-aa1f-f7dba5219cd7/wklife_paper_submitted.pdf)
# # 
# # In short, not using PE diagnostics increases the chance of systematically biased advice, under‑estimated risk, and delayed management response, particularly for data‑limited or structurally complex stocks.
# # 
# # The SP v SSB & Yield v SSB plots ask related but different questions: yield vs SSB focuses on realised catches at a given biomass under the observed exploitation pattern, whereas surplus production vs SSB focuses on the stock’s underlying production dynamics and productivity at that biomass. 
# # 
# # ## Yield vs SSB
# # 
# # - A yield–SSB curve shows equilibrium or realised **yield** as a function of spawning biomass under a given fishing regime or selectivity pattern. 
# # - Interpreting it is about trade‑offs between catch and biomass: reference points like \(F_{\text{MSY}}\), \(B_{\text{MSY}}\), and rebuilding targets are often derived here, so “where are we on the yield curve?” becomes the key inference. 
# # 
# # 
# # ## Surplus production vs SSB
# # 
# # - A surplus‑production–SSB plot shows the net **production** of the stock (growth + recruitment – natural losses – catches) as a function of biomass, typically peaking near \(B_{\text{MSY}}\) and changing sign around equilibrium. [patchwork.data-imaginist](https://patchwork.data-imaginist.com/articles/guides/layout.html)
# # - Inference here is about productivity and process error: whether the production curve is symmetric/asymmetric, whether recent points lie above or below the expected curve, and how environmental or regime shifts affect surplus production at a given biomass. 
# # 
# # 
# # ## Key inference differences
# # 
# # - Yield–SSB is catch‑centric and management oriented (what yield can be taken safely at a given SSB given current selectivity and \(F\)); surplus‑production–SSB is biology‑centric (how much the stock can produce at that SSB regardless of current \(F\)).
# # 
# # - Deviations of points from the yield curve may indicate changes in exploitation or implementation error, whereas deviations from the surplus‑production curve mainly indicate changes in stock productivity, process error, or model misspecification. 


