# SSB for a given recruitment level (returns both lower and upper solutions)

rickerSSBMaxRec<-function(b) {return(1/b)}

#' Generic for SSB at a given recruitment level for stock-recruitment models
setGeneric("rickerRec", function(params, rec, ...) standardGeneric("rickerRec"))
setGeneric("bevertonRec", function(params, rec, ...) standardGeneric("bevertonRec"))
setGeneric("segRegRec", function(params, rec, ...) standardGeneric("segRegRec"))

# Helper for rickerRec
.rickerRec_numeric <- function(params, rec) {
  rickerFun <- function(S) as.numeric(params["a"]) * S * exp(-as.numeric(params["b"]) * S) - rec
  lower <- tryCatch(uniroot(rickerFun, lower=.Machine$double.eps, upper=1/as.numeric(params["b"]))$root, error=function(e) NA)
  upper <- tryCatch(uniroot(rickerFun, lower=1/as.numeric(params["b"]), upper=10*(1/as.numeric(params["b"])))$root, error=function(e) NA)
  c(lower=lower, upper=upper)
}

setMethod("rickerRec", signature(params="FLPar", rec="numeric"),
  function(params, rec, ...) {
    sapply(rec, function(r) .rickerRec_numeric(params, r))
  })
setMethod("rickerRec", signature(params="FLPar", rec="FLPar"),
  function(params, rec, ...) {
    n_iter <- dims(rec)$iter
    res <- sapply(seq_len(n_iter), function(i) .rickerRec_numeric(params, as.numeric(rec[i])))
    dimnames(res) <- list(sol=c("lower","upper"), iter=seq_len(n_iter))
    FLPar(res)
  })
setMethod("rickerRec", signature(params="FLPar", rec="FLQuant"),
  function(params, rec, ...) {
    n_iter <- dim(rec)[6]
    res <- sapply(seq_len(n_iter), function(i) .rickerRec_numeric(params, as.numeric(rec@.Data[1,1,1,1,1,i])))
    dimnames(res) <- list(sol=c("lower","upper"), iter=seq_len(n_iter))
    FLQuant(res)
  })

# Helper for bevertonRec
.bevertonRec_numeric <- function(params, rec) {
  denom = as.numeric(params["a"]) - as.numeric(params["b"]) * rec
  if (denom <= 0) NA else rec/denom
}

setMethod("bevertonRec", signature(params="FLPar", rec="numeric"),
  function(params, rec, ...) {
    sapply(rec, function(r) .bevertonRec_numeric(params, r))
  })
setMethod("bevertonRec", signature(params="FLPar", rec="FLPar"),
  function(params, rec, ...) {
    n_iter <- dims(rec)$iter
    res <- sapply(seq_len(n_iter), function(i) .bevertonRec_numeric(params, as.numeric(rec[i])))
    names(res) <- seq_len(n_iter)
    FLPar(res)
  })
setMethod("bevertonRec", signature(params="FLPar", rec="FLQuant"),
  function(params, rec, ...) {
    n_iter <- dim(rec)[6]
    res <- sapply(seq_len(n_iter), function(i) .bevertonRec_numeric(params, as.numeric(rec@.Data[1,1,1,1,1,i])))
    names(res) <- seq_len(n_iter)
    FLQuant(res)
  })

# Helper for segRegRec
.segRegRec_numeric <- function(params, rec) {
  m1 <- as.numeric(params["m1"]); m2 <- as.numeric(params["m2"]); b1 <- as.numeric(params["b1"]); c <- as.numeric(params["c"])
  b2 <- m1*c+b1-m2*c
  if (rec <= m1*c+b1) (rec-b1)/m1 else (rec-b2)/m2
}

setMethod("segRegRec", signature(params="FLPar", rec="numeric"),
  function(params, rec, ...) {
    sapply(rec, function(r) .segRegRec_numeric(params, r))
  })
setMethod("segRegRec", signature(params="FLPar", rec="FLPar"),
  function(params, rec, ...) {
    n_iter <- dims(rec)$iter
    res <- sapply(seq_len(n_iter), function(i) .segRegRec_numeric(params, as.numeric(rec[i])))
    names(res) <- seq_len(n_iter)
    FLPar(res)
  })
setMethod("segRegRec", signature(params="FLPar", rec="FLQuant"),
  function(params, rec, ...) {
    n_iter <- dim(rec)[6]
    res <- sapply(seq_len(n_iter), function(i) .segRegRec_numeric(params, as.numeric(rec@.Data[1,1,1,1,1,i])))
    names(res) <- seq_len(n_iter)
    FLQuant(res)
  })

#' Generic for expected recruitment given SSB and an SRR object
setGeneric("recHat", function(object, ssb) standardGeneric("recHat"))

setMethod("recHat", signature(object="FLBRP", ssb="numeric"),
  function(object, ssb) {
    model_name <- as.character(model(object))
    pars <- params(object)
    if(model_name == "ricker") {
      return(as.numeric(pars["a"]) * ssb * exp(-as.numeric(pars["b"]) * ssb))
    } else if(model_name == "bevholt") {
      return((as.numeric(pars["a"]) * ssb) / (1 + as.numeric(pars["b"]) * ssb))
    } else if(model_name == "segreg") {
      m1 <- as.numeric(pars["m1"]); m2 <- as.numeric(pars["m2"]); b1 <- as.numeric(pars["b1"]); c <- as.numeric(pars["c"])
      b2 <- m1*c+b1-m2*c
      Rc <- m1*c+b1
      return(ifelse(ssb <= c, m1*ssb+b1, m2*ssb+b2))
    } else {
      stop("Unknown or unsupported SRR model in FLBRP")
    }
  })

setMethod("recHat", signature(object="FLSR", ssb="numeric"),
  function(object, ssb) {
    model_name <- as.character(model(object))
    pars <- params(object)
    if(model_name == "ricker") {
      return(as.numeric(pars["a"]) * ssb * exp(-as.numeric(pars["b"]) * ssb))
    } else if(model_name == "bevholt") {
      return((as.numeric(pars["a"]) * ssb) / (1 + as.numeric(pars["b"]) * ssb))
    } else if(model_name == "segreg") {
      m1 <- as.numeric(pars["m1"]); m2 <- as.numeric(pars["m2"]); b1 <- as.numeric(pars["b1"]); c <- as.numeric(pars["c"])
      b2 <- m1*c+b1-m2*c
      Rc <- m1*c+b1
      return(ifelse(ssb <= c, m1*ssb+b1, m2*ssb+b2))
    } else {
      stop("Unknown or unsupported SRR model in FLSR")
    }
  })

setMethod("recHat", signature(object="FLBRP", ssb="FLPar"),
  function(object, ssb) {
    model_name <- as.character(model(object))
    pars <- params(object)
    n_iter <- dims(ssb)$iter
    res <- rep(NA_real_, n_iter)
    if(model_name == "ricker") {
      for(i in seq_len(n_iter)) res[i] <- as.numeric(pars["a"]) * as.numeric(ssb[i]) * exp(-as.numeric(pars["b"]) * as.numeric(ssb[i]))
    } else if(model_name == "bevholt") {
      for(i in seq_len(n_iter)) res[i] <- (as.numeric(pars["a"]) * as.numeric(ssb[i])) / (1 + as.numeric(pars["b"]) * as.numeric(ssb[i]))
    } else if(model_name == "segreg") {
      m1 <- as.numeric(pars["m1"]); m2 <- as.numeric(pars["m2"]); b1 <- as.numeric(pars["b1"]); c <- as.numeric(pars["c"])
      b2 <- m1*c+b1-m2*c
      for(i in seq_len(n_iter)) res[i] <- if(as.numeric(ssb[i]) <= c) m1*as.numeric(ssb[i])+b1 else m2*as.numeric(ssb[i])+b2
    } else {
      stop("Unknown or unsupported SRR model in FLBRP")
    }
    return(FLPar(res, dimnames=list(iter=seq_len(n_iter))))
  })

setMethod("recHat", signature(object="FLSR", ssb="FLPar"),
  function(object, ssb) {
    model_name <- as.character(model(object))
    pars <- params(object)
    n_iter <- dims(ssb)$iter
    res <- rep(NA_real_, n_iter)
    if(model_name == "ricker") {
      for(i in seq_len(n_iter)) res[i] <- as.numeric(pars["a"]) * as.numeric(ssb[i]) * exp(-as.numeric(pars["b"]) * as.numeric(ssb[i]))
    } else if(model_name == "bevholt") {
      for(i in seq_len(n_iter)) res[i] <- (as.numeric(pars["a"]) * as.numeric(ssb[i])) / (1 + as.numeric(pars["b"]) * as.numeric(ssb[i]))
    } else if(model_name == "segreg") {
      m1 <- as.numeric(pars["m1"]); m2 <- as.numeric(pars["m2"]); b1 <- as.numeric(pars["b1"]); c <- as.numeric(pars["c"])
      b2 <- m1*c+b1-m2*c
      for(i in seq_len(n_iter)) res[i] <- if(as.numeric(ssb[i]) <= c) m1*as.numeric(ssb[i])+b1 else m2*as.numeric(ssb[i])+b2
    } else {
      stop("Unknown or unsupported SRR model in FLSR")
    }
    return(FLPar(res, dimnames=list(iter=seq_len(n_iter))))
  })

