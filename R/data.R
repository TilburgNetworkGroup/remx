#' @title networks
#' 
#' @description A set of 15 randomly generated networks. See https://github.com/TilburgNetworkGroup/remulate. 
#' 
#' @format A list containing 15 data frames, each data frame contains:
#' \describe{
#' \item{time}{the time at which each observation occurred}
#' \item{actor1}{the index representing the sender of the interaction}
#' \item{actor2}{the index represending the receiver of the interaction}
#' }
"networks"

#' @title stream
#' 
#' @description A network with 5000 events. See https://github.com/TilburgNetworkGroup/remulate. 
#' 
#' @format A data frame with 3 columns for time, sender and receiver
#' \describe{
#' \item{time}{the time at which each observation occurred}
#' \item{actor1}{the index representing the sender of the interaction}
#' \item{actor2}{the index represending the receiver of the interaction}
#' }
"stream"