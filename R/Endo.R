#'
#' @docType package
#' @name Endo
#' @format This package includes multiple datasets:
#' \describe{
#'    \item{IPMS_counts}{a dataframe of the processed Orbitrap counts data
#'                      (peptide counts) collected from 4 groups of duplicate
#'                      samples:
#'                      \item{(1) - TMEM with a myc tag (2 samples)}
#'                      \item{(2) - TMEM with a V5 tag (2 samples)}
#'                      \item{(3) - GFP with a myc tag (2 samples)}
#'                      \item{(4) - BirA with a V5 tag (2 samples)}}
#'    \item{CRAPome_results}{a dataframe containing the results of a CRAPome
#'                          query based on the wrangled IPMS counts data}
#'    \item{GO_CC_results_raw}{a dataframe read from a .txt file containing the
#'                            raw results of the cell component geneontology.org
#'                            query}
#'    \item{GO_CC_results}{a dataframe containing the data-cleaned results of
#'                        the cell compnent geneontology.org query}
#'    \item{GO_BP_results_raw}{a dataframe read from a .txt file containing the
#'                            raw results of the biological process
#'                            geneontology.org query}
#'    \item{GO_BP_results}{a dataframe containing the data-cleaned results of
#'                        the biological processes geneontology.org query}
#'    \item{GO_MF_results_raw}{a dataframe read from a .txt file containing the
#'                            raw results of the molecular function
#'                            geneontology.org query}
#'    \item{GO_MF_results}{a dataframe containing the data-cleaned results of
#'                        the molecular function geneontology.org query}
#'    \item{FIREpHly_WT}{a dataframe containing the microscopy data for
#'                      wild-type murine hippocampal neurons, measuring
#'                      ratiometric fluorescence and puncta size metrics of
#'                      endolysosomes}
#'    \item{FIREpHly_Mut}{a dataframe containing the microscopy data for
#'                        Tmem184b-mutant murine hippocampal neurons, measuring
#'                        ratiometric fluorescence and puncta size metrics of
#'                        endolysosomes}
#'    \item{tided_puncta_data}{a dataframe containing the data cleaned pH
#'                            (fluorescence) and puncta metric microscopy data}
#'    \item{PreCrapomeCandidateList}{a list containing the list of
#'                                  heuristic, putative TMEM184B interaction
#'                                  candidates}
#'    \item{PostCrapomeCandidateList}{a list containing the list of
#'                                    heuristic, putative TMEM184B interaction
#'                                    candidates after filtering by common
#'                                    false-positives in the CRAPome}
#' }
