# utils.R — the single definition of the null-coalescing helper used across the engine.
#
# `a %||% b` returns `b` when `a` is NULL or a single NA, otherwise `a`. It is sourced after vendor/
# (whose files carry their own NULL-only copies, kept byte-identical to their upstream sources), so this
# is the definition every engine function sees. The length-1 guard matters: config sub-lists are
# length > 1, and `is.na(list)` inside `||` errors on R >= 4.3.
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a
