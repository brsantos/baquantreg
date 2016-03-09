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

#' Gini index for Brazilian states in three censuses.
#'
#' A dataset containing data about the Gini index of 27 Brazilian states for
#'  three different censuses, 1991, 2000 and 2010, where each row represents
#'  the data about a specific state in a specific census year.
#'
#' @format A data frame with 81 rows and 4 variables:
#' \describe{
#'   \item{UFN}{Name of the state}
#'   \item{GINI}{Gini index}
#'   \item{INCPC}{Average income per capita}
#'   \item{EDUC}{Average years of education}
#'   \item{YEAR}{Year of the census}
#' }
#' @source \url{http://www.ibge.gov.br/english/estatistica/populacao/
#' default_censo_2000.shtm}
#' @references Santos and Bolfarine (2016). On Bayesian quantile regression
#'  and outliers. \url{http://arxiv.org/abs/1601.07344}
"BrazilGini"


#' Results of the Brazilian presidential election in 2014, with percentage of
#'  votes for Dilma Rousseff in each city of the country.
#'
#' A dataset containing data about the presidential election of 2014 in the
#'  total 5501 cities, where a few cities were not included in this dataset,
#'  for which it was not possible to match the information for votes and the
#'  other predictor variables.
#'
#' @format A data frame with 5501 rows and 13 variables:
#' \describe{
#'   \item{city}{Name of the city.}
#'   \item{UF}{Code for the state.}
#'   \item{UF_name}{Code of two letters for the state.}
#'   \item{Region}{Region of the city: 1 = North, 2 = Northeast,
#'   3 = Southeast, 4 = South, 5 = Midwest.}
#'   \item{income_delta}{Relative difference between income per capita in 2010
#'    and 2000.}
#'   \item{Gini_delta}{Relative difference between the Gini index in 2010 and
#'    2000.}
#'   \item{Years_educ}{Years of education.}
#'   \item{HDI}{Human development index.}
#'   \item{avepay_BF}{Average payment received for the participant families in
#'    the federal program for indirect cash transfer, called Bolsa-Familia.}
#'   \item{population}{Population number in 2010.}
#'   \item{latitude}{Latitude.}
#'   \item{longitude}{Longitude.}
#' }
#' @source \url{}
"BrazilElection2014"