#' Calculate pSup
#'
#' Calculate the proportion of each transcript in the Supernatant.
#' @details The input data should be the wide format in terms of fractions.
#' Function get_paras_median from this package is called in this function. Output
#' a data frame that includes a variable pSup calculated for each transcript.
#' @param wide_data A wider data frame
#' @param mixing_params output of get_mixing_params
#' @param chains A number, default to 4
#' @param iter A number, default to 1000
#' @param control A list of parameters, default to list(adapt_delta = 0.85)
#' @seealso [rstan::sampling()], `browseVignettes("rstan")`
#' @keywords pSup proportion supernatant pellet
#' @import dplyr
#' @importFrom dplyr left_join
#' @export
calculate_pSup <- function(wide_data, mixing_params, chains = 4,
                           iter = 1000, control = list(adapt_delta = 0.85))
{
  paras_median <- mixing_params %>%
    filter(Term=="x50") %>%
    select(Condition,Variable,Value) %>%
    pivot_wider(names_from=Variable,values_from=Value)
  wide_data_pSup <- wide_data %>%
    dplyr::left_join(select(paras_median,Condition,scaling_factor_sup,scaling_factor_pellet)) %>%
    mutate(pSup=scaling_factor_sup*Sup/(scaling_factor_sup*Sup+scaling_factor_pellet*Pellet)) %>%
    dplyr::rename("Total.counts"="Tot","Sup.counts"="Sup","Pellet.counts"="Pellet") %>%
    select(ORF,Condition,Total.counts,Sup.counts,Pellet.counts,scaling_factor_sup,scaling_factor_pellet,pSup)
}
