load_bascule = function(path="~/GitHub/bascule") {
  devtools::load_all(path)
}

load_simbascule = function(path="~/GitHub/simbascule") {
  devtools::load_all(path)
}

import_py = function(envname="bascule-env", path="~/GitHub/pybascule/") {
  reticulate::use_condaenv(envname)
  py <<- reticulate::import_from_path("pybascule", path)
}

load_deps = function(base_path="~/GitHub/") {
  import_py(path=paste0(base_path, "pybascule"))
  load_bascule(path=paste0(base_path, "bascule"))
  load_simbascule(path=paste0(base_path, "simbascule"))
  library(ggplot2)
}
