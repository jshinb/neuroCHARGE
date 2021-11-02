##
calc.relimp.lm.js <- function (object, type = "lmg", groups = NULL, groupnames = NULL, 
          always = NULL, ...) 
{
  require(relaimpo)
  lm <- object
  if (missing(lm)) 
    stop("object missing or incorrect")
  if (is.null(lm$terms)) 
    stop("object does not contain a terms component")
  if ("weight" %in% names(list(...)) || "weights" %in% names(list(...))) 
    stop(paste("For lm objects,", "\n", "only weights defined within the lm object are accepted.", 
               sep = ""))
  weights <- lm$weights
  if (is.null(lm$terms)) 
    stop("object does not contain a terms component")
  if (!inherits(lm, "lm")) 
    stop("object is not of class lm")
  if (inherits(lm, "mlm")) 
    stop("relaimpo does not work on multiresponse models")
  if (inherits(lm, "glm")) {
    if (!(lm$family$family == "gaussian" & lm$family$link == 
          "identity")) 
      stop("relaimpo works on linear models only (glms must be gaussian with identity link)")
  }
  terms <- lm$terms
  resp <- attr(terms, "response")
  if (attr(terms, "intercept") != 1) 
    stop("model must contain intercept")
  p <- length(labels(terms))
  if (is.null(always)) 
    numalw <- integer(0)
  if (is.character(always)) {
    if (any(!always %in% labels(terms))) 
      stop(paste("All names in always must refer to terms in the model."))
    else numalw <- which(labels(terms) %in% always) + 1
  }
  if (is.numeric(always)) {
    if (any(!always %in% 2:(p + 1))) 
      stop(paste("Numbers in always must refer to terms in model, \n i.e. here: 2 (for first non-intercept term) to ", 
                 (p + 1), ".", sep = ""))
    else {
      numalw <- always
      always <- labels(terms)[always - 1]
    }
  }
  if (!is.null(groups)) {
    if (any(c("betasq", "pratt") %in% type)) 
      stop("Metrics betasq and pratt do not work with groups.")
    if (!is.list(groups)) {
      if ((is.numeric(groups) || is.integer(groups)) && 
          !all(groups %in% 2:(p + 1))) 
        stop(paste("Numbers in groups must refer to columns 2 to ", 
                   p + 1, " in Y,X1,...,Xp.", sep = ""))
      if (length(groups) <= 1) 
        stop("groups must list groups of more than one regressor.")
      if (!(is.character(groups) || is.numeric(groups))) 
        stop("groups must be numeric or character")
      groups <- list(groups)
    }
    else {
      lapply(groups, function(obj) {
        if (is.numeric(obj) || is.integer(obj)) {
          if (!all(obj %in% 2:(p + 1))) 
            stop(paste("Numbers in elements of groups must refer to columns 2 to ", 
                       p + 1, " in Y,X1,...,Xp.", sep = ""))
        }
        else if (!is.character(obj)) 
          stop("groups must be numeric or character")
      })
    }
    if (!is.null(groupnames)) {
      if (!length(groups) == length(groupnames)) 
        stop(paste("groupnames must have one entry for each group.", 
                   "\n", "There are", length(groups), "groups and", 
                   length(groupnames), "group names."))
    }
  }
  ogroups <- groups
  ngroups <- groups
  WW <- NULL
  if (any(type != "lmg") && max(attr(terms, "order")) > 1) 
    stop("Higher order terms supported for lmg only!")
  if (max(attr(terms, "order")) >= 2) {
    if (!is.null(groups)) 
      stop("Currently, groups and interactions cannot be used at the same time.")
    if (max(attr(terms, "order")) > 2) 
      stop("Currently, only 2-variable interactions are supported.")
    if (sum(attr(lm$terms, "order") > 1) > 1) {
      WW.mat <- attr(lm$terms, "factors")[-1, attr(lm$terms, 
                                                   "order") > 1]
    }
    else {
      WW.mat <- as.matrix(attr(lm$terms, "factors")[-1, 
                                                    attr(lm$terms, "order") > 1])
      colnames(WW.mat) <- colnames(attr(lm$terms, "factors"))[attr(lm$terms, 
                                                                   "order") > 1]
    }
    WW.HK <- sum(attr(lm$terms, "order") == 1)
    WW <- cbind(posWW = col(WW.mat)[WW.mat == 1] + WW.HK, 
                posHK = row(WW.mat)[WW.mat == 1])
    order <- attr(lm$terms, "order")
    WW <- list(WW = WW, Order = order)
    WW$WW <- WW$WW[, c(1, 2, 2)]
    DATA = cbind(lm$model[, 1], model.matrix(lm)[, -1])
    colnames(DATA)[1] <- colnames(lm$model)[1]
    g_facnames <- character(0)
    factor.var <- which(as.vector(attr(lm$terms, "dataClasses") == 
                                    "factor")) - 1
    if (!is.null(groups)) {
      groups <- lapply(groups, function(obj) {
        if (is.character(obj)) 
          which(labels(terms) %in% obj) + 1
        else obj
      })
      ngroups <- groups
    }
    if (any(factor.var %in% WW$WW[, 2]) || length(lm$assign) > 
        length(WW[[2]]) + 1) {
      if (length(lm$assign) > length(WW[[2]]) + 1) {
        hilf <- matrix(0, 0, 3)
        for (we in unique(WW$WW[, 1])) {
          hes <- which(lm$assign %in% WW$WW[WW$WW[, 1] == 
                                              we, 2]) - 1
          hess <- lm$assign[hes + 1]
          hilf <- rbind(hilf, matrix(c(rep(we, length(hes)), 
                                       hes, hess), length(hes), 3, byrow = FALSE))
        }
        WW$WW <- hilf
      }
      y.modelm <- model.matrix(lm)
      g_names <- attr(y.modelm, "dimnames")[[2]]
      g_assign <- attr(y.modelm, "assign")
      g_facpos <- as.data.frame(table(g_assign))
      if (!is.null(groups)) {
        ogroups <- lapply(groups, function(obj) {
          c(which(g_assign %in% (obj - 1)))
        })
        groups <- lapply(ogroups, function(obj) {
          g_names[obj]
        })
      }
      for (k in 2:length(g_facpos[, 1])) {
        fac.name <- labels(terms)[as.numeric(g_facpos[k, 
                                                      1]) - 1]
        fac.varn <- g_names[g_assign == g_facpos[k, 1]]
        if (g_facpos[k, 2] == 1) {
          colnames(DATA)[g_assign == g_facpos[k, 1]] <- fac.name
          if (!is.null(groups)) 
            groups <- lapply(groups, function(obj) {
              if (fac.varn %in% obj) 
                obj[obj == fac.varn] <- fac.name
              obj
            })
          next
        }
        if (any(always %in% fac.name)) {
          always <- c(always[always != fac.name], fac.varn)
          next
        }
        if (!is.null(groups)) {
          schonda <- FALSE
          for (obj in groups) {
            if (all(fac.varn %in% obj)) 
              schonda <- TRUE
          }
          if (schonda) 
            next
        }
        g_facnames <- c(g_facnames, fac.name)
        if (is.null(groups)) {
          groups <- list(fac.varn)
          ngroups <- list(which(labels(terms) %in% fac.name) + 
                            1)
        }
        else {
          groups[[length(groups) + 1]] <- fac.varn
          ngroups[[length(ngroups) + 1]] <- which(labels(terms) %in% 
                                                    fac.name) + 1
        }
      }
      if (!is.null(groups)) 
        groups <- lapply(groups, function(obj) {
          which(colnames(DATA) %in% obj)
        })
      if (is.null(ogroups)) 
        groupnames <- g_facnames
      else {
        if (is.null(groupnames)) 
          groupnames <- c(paste("G", 1:length(ogroups), 
                                sep = ""), g_facnames)
        else groupnames <- c(groupnames, g_facnames)
      }
    }
    hilf <- which(colnames(DATA) %in% always) - 1
    if (any(hilf %in% WW$WW[, 2])) {
      for (we in intersect(hilf, unique(WW$WW[, 2]))) {
        if (!(is.matrix(WW$WW))) 
          WW$WW <- matrix(WW$WW, 1, 3)
        WW$WW <- WW$WW[which(!WW$WW[, 2] == we), ]
      }
      if (!(is.matrix(WW$WW))) 
        WW$WW <- matrix(WW$WW, 1, 3)
      if (nrow(WW$WW) == 0) 
        WW <- NULL
    }
    if (!is.null(WW)) {
      if (any(hilf %in% which(lm$assign %in% WW$WW[, 1]))) 
        stop("Interaction must not be in always, \n \n                 if not all corresponding lower level effects are also in always.")
      WW$WW <- WW$WW[, c(1, 3)]
      if (!(is.matrix(WW$WW))) 
        WW$WW <- matrix(WW$WW, 1, 2)
      for (we in sort(unique(WW$WW[, 1]))) {
        hilf <- unique(WW$WW[which(WW$WW[, 1] == we), 
                             2])
        WW$WW <- rbind(WW$WW[which(WW$WW[, 1] < we), 
        ], cbind(rep(we, length(hilf)), hilf), WW$WW[which(WW$WW[, 
                                                                 1] > we), ])
      }
    }
    if (!is.null(ngroups)) {
      ngroups <- append(ngroups, as.list(setdiff(2:(p + 
                                                      1), c(list2vec(ngroups), numalw))))
      ngroups <- lapply(ngroups, function(obj) {
        if (length(obj) > 1) {
          if (any((obj - 1) %in% WW$WW)) 
            stop("Groups must not contain effects involved in interactions.")
          else obj <- obj[1]
        }
        obj
      })
    }
    if (length(numalw) > 0 && !is.null(ngroups)) 
      ngroups <- lapply(ngroups, function(obj) {
        obj - rowSums(matrix(obj, length(obj), length(numalw), 
                             byrow = F) > matrix(numalw, length(obj), length(numalw), 
                                                 byrow = T))
      })
    if (!is.null(ngroups)) {
      if (length(ngroups) == 0) 
        ngroups <- NULL
    }
    y <- do.call("calc.relimp.default.intern", list(DATA, 
                                                    weights = weights, groups = groups, groupnames = groupnames, 
                                                    WW = WW, always = always, ngroups = ngroups, ...))
  }
  else {
    # DATA = cbind(lm$model[, 1], model.matrix(lm)[, -1])
    # DATA <- as.data.frame(lm$model[, c(resp, which(rowSums(attr(terms,"factors")) > 0))])
    DATA <- as.data.frame(model.matrix(lm)[, c(resp,which(rowSums(attr(terms,"factors")) > 0))])
    DATA[,1] <- lm$model[,resp]
    colnames(DATA)[1] <- names(lm$model)[resp]
    rownames(DATA) = rownames(model.matrix(lm))
    y <- do.call("calc.relimp", list(DATA, weights = weights, 
                                     type = type, groups = groups, groupnames = groupnames, 
                                     always = always, ...))
  }
  y
}
##
environment(calc.relimp.lm.js) <- asNamespace('relaimpo')
assignInNamespace("calc.relimp.lm", calc.relimp.lm.js, ns = "relaimpo")
