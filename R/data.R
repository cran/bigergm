#' A toy network to play `bigergm` with.
#'
#' This network has a clear cluster structure.
#' The number of clusters is four, and which cluster each node belongs to is defined in the variable "block".
#'
#' @format A `statnet`'s network class object. It has three nodal features.
#' \describe{
#'   \item{block}{block membership of each node}
#'   \item{x}{a covariate. It has 10 labels.}
#'   \item{y}{a covariate. It has 10 labels.}
#'   ...
#' }
#' `1` and `2` are not variables with any particular meaning.
#'
"toyNet"

#' A network of friendships between students at Reed College.
#' The data was collected by Facebook and provided as part of Traud et al. (2012)
#' @references 
#' Traud, Mucha, Porter (2012). Social Structure of Facebook Network. 
#' Physica A: Statistical Mechanics and its Applications, 391, 4165-4180
#'
#' @format A `statnet`'s network class object. It has three nodal features.
#' \describe{
#'   \item{doorm}{anonymized dorm in which each node lives.}
#'   \item{gender}{gender of each node.}
#'   \item{high.school}{anonymized highschool to which each node went to.}
#'   \item{year}{year of graduation of each node.}
#'   ...
#' }
#'
"reed"

#' A network of friendships between students at Rice University.
#'
#' The data was collected by Facebook and provided as part of Traud et al. (2012)
#' @references 
#' Traud, Mucha, Porter (2012). Social Structure of Facebook Network. 
#' Physica A: Statistical Mechanics and its Applications, 391, 4165-4180
#' @format A `statnet`'s network class object. It has three nodal features.
#' \describe{
#'   \item{doorm}{anonymized dorm in which each node lives.}
#'   \item{gender}{gender of each node.}
#'   \item{high.school}{anonymized highschool to which each node went to.}
#'   \item{year}{year of graduation of each node.}
#' }
#'
"rice"

#' Bali terrorist network
#'
#' @description
#' The network corresponds to the contacts between the 17 terrorists 
#' who carried out the bombing in Bali, Indonesia in 2002.
#' The network is taken from Koschade (2006).
#'
#' @references 
#' Koschade, S. (2006). A social network analysis of Jemaah Islamiyah: The applications to counter-terrorism and intelligence. 
#' Studies in Conflict and Terrorism, 29, 559--575.
#' 
#' @format A `statnet`'s network class object. 
#'
"bali"

#' Van de Bunt friendship network
#'
#' @description
#' Van de Bunt (1999) and Van de Bunt et al. (1999) 
#' collected data on friendships between 32 freshmen at a European university at 7 time points.
#' Here, the last time point is used.
#' A directed edge from student \code{i} to \code{j} indicates that student \code{i} 
#' considers student \code{j} to be a ``friend" or ``best friend".
#'
#' @references 
#' Van de Bunt, G. G. (1999). Friends by choice. An Actor-Oriented Statistical Network Model for Friendship Networks through Time. Thesis Publishers, Amsterdam.
#' 
#' Van de Bunt, G. G., Van Duijn, M. A. J., and T. A. B. Snijders (1999). Friendship Networks Through Time: An Actor-Oriented Statistical Network Model. Computational and Mathematical Organization Theory, 5, 167--192.
#' 
#' @format A `statnet`'s network class object. 
#'
"bunt"


#' Kapferer collaboration network
#'
#' @description
#' The network corresponds to collaborations between 39 workers in a tailor shop in Africa:
#' an undirected edge between workers \code{i} and \code{j} indicates that the workers collaborated.
#' The network is taken from Kapferer (1972).
#'
#' @references 
#' Kapferer, B. (1972). Strategy and Transaction in an African Factory. Manchester University Press, Manchester, U.K.
#' 
#' @format A `statnet`'s network class object. 
#'
"kapferer"


