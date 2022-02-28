# Set library path (function from https://milesmcbain.xyz/hacking-r-library-paths/)
"set_lib_paths" <- function(lib_vec) {
    
    lib_vec <- normalizePath(lib_vec, mustWork = TRUE)
    
    shim_fun <- .libPaths
    shim_env <- new.env(parent = environment(shim_fun))
    shim_env$.Library <- character()
    shim_env$.Library.site <- character()
    
    environment(shim_fun) <- shim_env
    shim_fun(lib_vec)
    
}

