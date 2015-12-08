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

#' Proportion of households with access to electricity in Brazil
#'
#' A dataset containing a sample of cities in Brazil, with the proportion of
#' households with access to electricity with data from the Brazilian census
#' of 2000, and some sociodemographic variables.
#'
#' @format A data frame with 500 rows and 6 variables:
#' \describe{
#'   \item{prop_elec}{propotion of households with access to electricity.}
#'   \item{region}{2 for Northeast, 3 for Southeast.}
#'   \item{population}{Population, in years.}
#'   \item{income_percap}{Income per capita, in Brazilian reais.}
#'   \item{hdi}{Human development index.}
#'   \item{pop_density}{Population density.}
#' }
#' @source \url{http://www.ibge.gov.br/english/estatistica/populacao/
#' default_censo_2000.shtm}
#' @references Santos and Bolfarine (2015). Bayesian analysis for zero-or-one
#'  inflated proportion data using quantile regression. \emph{Journal of
#'  Statistical Computation and Simulation}. \url{http://www.tandfonline.com
#'  /doi/abs/10.1080/00949655.2014.986733}
"BrazilElectricity"