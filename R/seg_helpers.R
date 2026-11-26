################################################################################
# Segment dating check for xd_check(): optimal offset path per series x reference
# (dynamic programming on cumulative Pearson contributions)

# Help from Claude, 24 Nov 2026

# Optimal offset path (Viterbi). One offset per series year; maximizes the summed
# cor.contr minus penalty * |jump| at every change of offset.
# S: years x lags matrix of cor.contr (NA -> 0); lags: integer lag of each column.
lag_path <- function(S, lags, penalty = 4) {
  n <- nrow(S); K <- ncol(S)
  V <- matrix(0, n, K); B <- matrix(0L, n, K)
  V[1, ] <- S[1, ]
  sw <- penalty * abs(outer(lags, lags, "-"))   # rows = from, cols = to
  if (n > 1) {
    for (t in 2:n) {
      cand <- V[t - 1, ] - sw
      B[t, ] <- max.col(t(cand), ties.method = "first")
      V[t, ] <- cand[cbind(B[t, ], seq_len(K))] + S[t, ]
    }
  }
  path <- integer(n); path[n] <- which.max(V[n, ])
  if (n > 1) for (t in n:2) path[t - 1] <- B[t, path[t]]
  lags[path]
}

# Offset segments for ONE series vs ONE reference.
# lag.dfs: the per-series list of lagged data.tables built in xd_check
# (cols: year [= chron year], series, std.series, std.chron, cor.contr, lag).
seg_path <- function(lag.dfs, penalty = 4, min.seg = 20, p.thresh = 0.05) {
  if (!length(lag.dfs)) return(NULL)
  al <- data.table::rbindlist(lag.dfs, fill = TRUE)
  al[, lag := as.integer(lag)]
  if (!0L %in% al$lag) return(NULL)            # need as-dated for comparison
  al[, series.year := year - lag]

  S <- data.table::dcast(al, series.year ~ lag, value.var = "cor.contr")
  yrs  <- S$series.year
  lags <- as.integer(colnames(S)[-1])
  S <- as.matrix(S[, -1]); S[is.na(S)] <- 0

  path <- data.table::data.table(series.year = yrs,
                                 offset = lag_path(S, lags, penalty))
  segs <- path[, .(seg.start = min(series.year), seg.end = max(series.year), seg.len = .N),
               by = .(offset, grp = data.table::rleid(offset))]
  segs <- segs[offset != 0L & seg.end - seg.start + 1 >= min.seg]
  if (!nrow(segs)) return(NULL)
  segs[, grp := NULL]
  # is the whole rest of the series off?
  segs[, whole.series := seg.start == min(yrs)]

  # Spearman over each segment: offset vs as dated
  l0 <- al[lag == 0L, .(series.year, chron0 = std.chron)]
  segs[, c("as.dated.corr", "as.dated.p", "offset.corr", "offset.p") := {
    k <- offset; a <- seg.start; b <- seg.end
    s <- l0[al[lag == k & series.year %between% c(a, b)], on = "series.year"]
    ct <- suppressWarnings(cor.test(s$std.series, s$std.chron, method = "spearman"))
    r0 <- tryCatch(suppressWarnings(cor.test(s$std.series, s$chron0, method = "spearman",
                                             use = "complete.obs")),
                   error = \(e) NA_real_)
    list(round(r0$estimate[[1]], 2), round(r0$p.value, 3),
         round(ct$estimate[[1]], 2), round(ct$p.value, 3))
  }, by = .(offset, seg.start)]

  segs[, series := al$series[1]]
  # Subset to just the segs that have a better (and significant) offse.corr than as.dated.corr
  segs[(offset.corr > as.dated.corr | is.na(as.dated.corr)) & offset.p < p.thresh]
}

# Agreement across references: a segment is supported by another reference if
# that reference flags the same offset over >= min.frac of the segment's years.
seg_agree <- function(seg.res, n.refs, min.frac = 0.5) {
  if (!nrow(seg.res)) return(seg.res)
  seg.res[, seg.id := .I]
  yrs <- seg.res[, .(series.year = seg.start:seg.end),
                 by = .(seg.id, series, reference, offset)]
  m <- yrs[yrs, on = .(series, offset, series.year),
           allow.cartesian = TRUE, nomatch = 0L] # i.* = focal segment
  cov <- m[, .(n.yrs = .N), by = .(seg.id = i.seg.id, other = reference)]
  cov[seg.res, on = "seg.id", len := i.seg.end - i.seg.start + 1]
  agree <- cov[n.yrs / len >= min.frac,
               .(n.agree = data.table::uniqueN(other),
                 agree.refs = paste(sort(unique(other)), collapse = ", ")),
               by = seg.id]
  seg.res[agree, on = "seg.id",
          `:=`(n.agree = i.n.agree, agree.refs = i.agree.refs)]
  seg.res[, agreement := paste(n.agree, n.refs, sep = "/")]
  data.table::setorder(seg.res, series, -n.agree, seg.start)
  seg.res[, .(series, reference, offset, seg.start, seg.end, seg.len, whole.series,
              as.dated.corr, as.dated.p, offset.corr, offset.p, agreement, agree.refs)]
}
