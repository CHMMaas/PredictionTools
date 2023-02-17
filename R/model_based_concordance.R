#' @title mb.c
#' @description This function calculates the model-based concordance based on predicted probabilities.
#'
#'
#' @param p.hat Predicted probabilities
#'
#' @return The output of the mb.c function is:
#'
#' mb.c
#'
#' Model-based concordance
#'
#' @export
#'
#' @examples
#' library(PredictionTools)
#' set.seed(1)
#' n <- 100
#' p.hat <- runif(n, 0, 1)
#' PredictionTools::mb.c(p.hat=p.hat)
mb.c <- function(p.hat){
  stopifnot("p.hat must be numeric" = is.numeric(p.hat))
  stopifnot("p.hat must be a vector" = is.vector(p.hat))

  n<-length(p.hat)
  ord<-order(p.hat)
  p.hat<-p.hat[ord]
  q.hat<-1-p.hat
  V1<-(p.hat*(cumsum(q.hat)-q.hat)+q.hat*(sum(p.hat)-cumsum(p.hat)))/(n-1)
  V2<-(p.hat*(sum(q.hat)-q.hat)+q.hat*(sum(p.hat)-p.hat))/(n-1)
  mb.c<-sum(V1)/sum(V2)
  return(mb.c)
}
