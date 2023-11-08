#' Applies method $ to different classes
#'
#' Enables slots of objects in this package to be referenced easily.
#' @docType methods
#' @name getter
#' @param x Object.
#' @param name Name of slot.
#'
#' @keywords internal
NULL

#' @rdname getter
setMethod("$",
          "mdIVW",
          function(x, name) slot(x, name))
