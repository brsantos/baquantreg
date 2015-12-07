#' Expenditures with durable goods in Brazil.
#'
#' A dataset containing the expenditures with durable goods for the state of
#' Maranhão, in Brazil, from the "Consumer Expenditure Survey", which took
#' place between 2008 and 2009.
#'
#' @format A data frame with 2240 rows and 6 variables:
#' \describe{
#'   \item{expenditure}{expenditure, in Brazilian reais}
#'   \item{gender}{1 for male, 2 for female}
#'   \item{age}{Age, in years.}
#'   \item{race}{1 for white, 2 for non-white}
#'   \item{credit_card}{1 if the observation has credit card, 2 otherwise.}
#'   \item{education}{Education, in years.}
#' }
#' @source \url{http://www.ibge.gov.br/english/estatistica/populacao/
#' condicaodevida/pof/2008_2009/default.shtm}
#' @references Santos and Bolfarine (2015) - Bayesian quantile regression
#' analysis for continuous data with a discrete component at zero.
#' \emph{Preprint}. \url{http://arxiv.org/abs/1511.05925}
"BrazilDurableGoods"