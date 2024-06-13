#' Twitter (X) network of U.S. state legislators
#'
#' @description
#' The network includes the Twitter (X) following interactions between U.S. state legislators. 
#' The data was collection by Gopal et al. (2022) and Kim et al. (2022).
#' For this network, we only include the largest connected component of state legislators that were active on Twitter in the six months leading up to and including the insurrection at the United States Capitol on January 6, 2021.
#' All state senate and state representatives for states with a bicameral system are included and all state legislators for state (Nebraska) with a unicameral system are included.
#' 
#' @name state_twitter
#' 
#' @references 
#' Gopal, Kim, Nakka, Boehmke, Harden, Desmarais. 
#' The National Network of U.S. State Legislators on Twitter. 
#' Political Science Research & Methods, Forthcoming.
#' 
#' Kim, Nakka, Gopal, Desmarais,Mancinelli, Harden, Ko, and Boehmke (2022). 
#' Attention to the COVID-19 pandemic on Twitter: Partisan differences among U.S. state legislators. Legislative Studies Quarterly 47, 1023–1041.
#'
#' @format A `statnet`'s network class object. It has the following categorical attributes for each state legislator.
#' \describe{
#'   \item{gender}{factor stating whether the legislator is 'female' or 'male'.}
#'   \item{party}{party affiliation of the legislator, which is 'Democratic', 'Independent' or 'Republican'.}
#'   \item{race}{race with the following levels: 'Asian or Pacific Islander', 'Black', 'Latino', 'MENA(Middle East and North Africa)','Multiracial', 'Native American', and 'White'.}
#'   \item{state}{character of the state that the legislator represents.}
#' }
#' @usage
#' data(state_twitter)
NULL


#' A toy network to play `bigergm` with.
#'
#' @description
#' This network has a clear cluster structure.
#' The number of clusters is four, and which cluster each node belongs to is defined in the variable "block".
#' @name toyNet
#' @format A `statnet`'s network class object. It has three nodal features.
#' \describe{
#'   \item{block}{block membership of each node}
#'   \item{x}{a covariate. It has 10 labels.}
#'   \item{y}{a covariate. It has 10 labels.}
#'   ...
#' }
#' `1` and `2` are not variables with any particular meaning.
#' @usage
#' data(toyNet)
NULL


#' A network of friendships between students at Reed College.
#' 
#' @name reed
#' @description
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
#' data(reed)
NULL

#' A network of friendships between students at Rice University.
#'
#'
#' @name rice
#' @description
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
#' data(rice)
NULL

#' Bali terrorist network
#'
#' @description
#' The network corresponds to the contacts between the 17 terrorists 
#' who carried out the bombing in Bali, Indonesia in 2002.
#' The network is taken from Koschade (2006).
#' @name bali
#'
#' @references 
#' Koschade, S. (2006). A social network analysis of Jemaah Islamiyah: The applications to counter-terrorism and intelligence. 
#' Studies in Conflict and Terrorism, 29, 559--575.
#' 
#' @format A `statnet`'s network class object. 
#' data(bali)
NULL

#' Van de Bunt friendship network
#'
#' @description
#' Van de Bunt (1999) and Van de Bunt et al. (1999) 
#' collected data on friendships between 32 freshmen at a European university at 7 time points.
#' Here, the last time point is used.
#' A directed edge from student \code{i} to \code{j} indicates that student \code{i} 
#' considers student \code{j} to be a ``friend" or ``best friend".
#' @name bunt
#'
#' @references 
#' Van de Bunt, G. G. (1999). Friends by choice. An Actor-Oriented Statistical Network Model for Friendship Networks through Time. Thesis Publishers, Amsterdam.
#' 
#' Van de Bunt, G. G., Van Duijn, M. A. J., and T. A. B. Snijders (1999). Friendship Networks Through Time: An Actor-Oriented Statistical Network Model. Computational and Mathematical Organization Theory, 5, 167--192.
#' 
#' @format A `statnet`'s network class object. 
#' data(bunt)
NULL

#' Kapferer collaboration network
#'
#' @description
#' The network corresponds to collaborations between 39 workers in a tailor shop in Africa:
#' an undirected edge between workers \code{i} and \code{j} indicates that the workers collaborated.
#' The network is taken from Kapferer (1972).
#'
#' @references 
#' Kapferer, B. (1972). Strategy and Transaction in an African Factory. Manchester University Press, Manchester, U.K.
#' @name kapferer
#' @format A `statnet`'s network class object. 
#' data(kapferer)
NULL



