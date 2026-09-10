# dsSemiOPBART
dsSemiOPBART: A Federated Linear Regression and Soft Bayesian Addictive Regression Tree for for DataSHIELD serverSh


cd ~/Datashield_package_dev/FedBART/dsSemiOPBART/SoftBart

R CMD INSTALL --preclean --no-multiarch .



# devtools::clean_dll(
#   "~/Datashield_package_dev/FedBART/dsSemiOPBART/SoftBart"
# )

.rs.restartR()

remove.packages("SoftBart")

unlink(
  file.path(
    .libPaths()[1],
    "SoftBart"
  ),
  recursive = TRUE,
  force = TRUE
)

unlink(
  file.path(
    .libPaths()[1],
    "SoftBart.rdb"
  ),
  recursive = TRUE,
  force = TRUE
)

unlink(
  file.path(
    .libPaths()[1],
    "SoftBart.rds"
  ),
  recursive = TRUE,
  force = TRUE
)

devtools::install(
  "~/Datashield_package_dev/FedBART/dsSemiOPBART/SoftBart",
  upgrade = FALSE
)

# devtools::clean_dll(
#   "~/Datashield_package_dev/FedBART/dsSemiOPBART/SoftBart"
# )


# library(SoftBart)
# hypers = Hypers(X = df_imputed[c("CRP","ESR")], Y = df_imputed[["IP.probability_0_4_BM"]])
# 
# opts = Opts()
# 
# smopbart_forest = MakeForest(hypers, opts, FALSE)