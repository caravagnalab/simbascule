load_basilica = function(path="~/GitHub/basilica") {
  devtools::load_all(path)
}

load_simbasilica = function(path="~/GitHub/simbasilica") {
  devtools::load_all(path)
}

import_py = function(envname="basilica-env", path="~/GitHub/pybasilica/") {
  reticulate::use_condaenv(envname)
  py <<- reticulate::import_from_path("pybasilica", path)
}

load_deps = function(base_path="~/GitHub/") {
  import_py(path=paste0(base_path, "pybasilica"))
  load_basilica(path=paste0(base_path, "basilica"))
  load_simbasilica(path=paste0(base_path, "simbasilica"))
  library(ggplot2)
}
