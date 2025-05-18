
devtools::document(pkg = 'C:/Users/Erik/Desktop/Programming/R/Bio/Endo')
devtools::build(pkg = 'C:/Users/Erik/Desktop/Programming/R/Bio/Endo', path = 'C:/Users/Erik/Desktop/Programming/R/Bio/Endo', binary = T, vignettes = F)

devtools::check(pkg = 'C:/Users/Erik/Desktop/Programming/R/Bio/Endo', document = T, vignettes = F, cran = T)
devtools::load_all(export_all = T)

devtools::install('C:/Users/Erik/Desktop/Programming/R/Bio/Endo', reload = T, build = T, dependencies = T)
devtools::install('C:/Users/Erik/Desktop/Programming/R/Bio/Endo', reload = T, build = F, dependencies = T)

# pkgbuild::check_build_tools(debug = TRUE)

usethis::use_news_md()

library(testthat)
usethis::use_github_action("check-standard")
install.packages("covr")
library(covr)

usethis::use_news_md()
usethis::use_vignette("Endo")
testthat::test_that("Endo")

roxygen2::roxygenize(clean = TRUE, package.dir = 'C:/Users/Erik/Desktop/Programming/R/Bio/Endo')
# data load ----
library(tidyverse)
IPMS_counts <- utils::read.csv(
  "C:/Users/Erik/Desktop/Programming/R/Bio/Endo/20220307_4552_4553_BirA_GFP_Samples_View_Report.csv"
)
usethis::use_data(IPMS_counts)

CRAPome_results <- readr::read_table(
  file = 'C:/Users/Erik/Desktop/Programming/R/Bio/Endo/CRAPome_results_272.txt'
)
usethis::use_data(CRAPome_results)

GO_CC_results <- utils::read.table(
  'C:/Users/Erik/Desktop/Programming/R/Bio/Endo/GO_CC_complete_fisher_FDR_all_136.txt',
  sep = '\t',
  skip = 11,
  header = T,
  quote = '')
usethis::use_data(GO_CC_results)

GO_BP_results <- utils::read.table(
  'C:/Users/Erik/Desktop/Programming/R/Bio/Endo/GO_BP_complete_fisher_FDR_all_136.txt',
  sep = '\t',
  skip = 11,
  header = T,
  quote = '')
usethis::use_data(GO_BP_results)

GO_MF_results <- utils::read.table(
  'C:/Users/Erik/Desktop/Programming/R/Bio/Endo/GO_MF_complete_fisher_FDR_all_136.txt',
  sep = '\t',
  skip = 11,
  header = T,
  quote = '')
usethis::use_data(GO_MF_results)

library(readxl)
FIREpHly_WT <- readxl::read_xlsx(
  'C:\\Users\\Erik\\Downloads\\2024-11-18_FIRE-pHly_Nikon Analysis .xlsx',
  col_names = F,
  sheet = 1)
usethis::use_data(FIREpHly_WT)

FIREpHly_Mut <- readxl::read_xlsx(
  'C:\\Users\\Erik\\Downloads\\2024-11-18_FIRE-pHly_Nikon Analysis .xlsx',
  col_names = F,
  sheet = 2)
usethis::use_data(FIREpHly_Mut)

# testing ----
usethis::use_testthat()
usethis::use_test()

results <- usethis::use_cran_comments()

usethis::use_travis()

usethis::use_github_action()

#install rhub
install.packages("rhub")

# check and release ----
#sign in
rhub::validate_email()
# check
devtools::check_rhub()

devtools::check_win_devel()

cran_checks <- rhub::check_for_cran()

#check package viability on all OSs
devtools::check_win_devel()

results$cran_summary()

devtools::test_coverage()

#release to cran
devtools::release()
devtools::spell_check()

##
gitcreds::gitcreds_set(url = 'https://www.github.com/eriklarsen4/Endo')
rhub::rhub_check(gh_url = 'https://www.github.com/eriklarsen4/Endo')
usethis::create_github_token()

usethis::use_spell_check(vignettes = TRUE, lang = 'en-US', error = FALSE)

# misc ----
install.packages("annotationDbi", lib = 'C:/Users/Erik/AppData/Local/R/win-library/4.3')
BiocManager::install("AnnotationDbi")
BiocManager::install("AnnotationDbi", lib = 'C:/Users/Erik/AppData/Local/R/win-library/4.3')

# ----
