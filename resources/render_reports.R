# Batch rendering of .R scripts to markdown reports (GitHub-friendly)
packages_needed <- c("rmarkdown", "tools", "rprojroot", "knitr")
if (any(!packages_needed %in% rownames(installed.packages()))) {
  install.packages(setdiff(packages_needed, rownames(installed.packages())))
}
lapply(packages_needed, library, character.only = TRUE)

# find your project root
proj_root <- rprojroot::find_rstudio_root_file()

# treat proj_root as the base for all figure paths
knitr::opts_knit$set(
  root.dir = proj_root,
  base.dir = proj_root
)

# list all .R scripts in /code
scripts <- list.files(file.path(proj_root, "code"),
                      pattern = "\\.R$", full.names = TRUE)

for (script in scripts) {
  script_base <- tools::file_path_sans_ext(basename(script))
  output_md  <- paste0(script_base, ".md")
  
  # define the *relative* figure folder including the default "figure-gfm" subdir
  fig_rel <- file.path("resources",
                       paste0(script_base, "_files"),
                       "figure-gfm",
                       "")
  # make sure it exists under your repo root
  dir.create(file.path(proj_root, fig_rel),
             recursive = TRUE,
             showWarnings = FALSE)
  
  # set knitr so that knitting happens from proj_root
  knitr::opts_knit$set(root.dir = proj_root)
  # tell all chunks to dump their plots into our new folder
  knitr::opts_chunk$set(fig.path = fig_rel)
  
  # render with knit_root_dir as proj_root, and disable self_contained so images stay external
  # wrap with tryCatch to help debug
  tryCatch({
    rmarkdown::render(
      input          = script,
      # output_format  = rmarkdown::github_document(),
      output_file    = output_md,
      output_dir     = proj_root,
      knit_root_dir  = proj_root,
      envir          = new.env(parent = globalenv())
    )
  }, error = function(e) {
    message("ERROR in ", basename(script), ": ", conditionMessage(e))
    stop(e)
  })
  
  # --- MOVE the entire <script>_files/ folder into resources/ ---
  orig_dir <- file.path(proj_root, paste0(script_base, "_files"))
  dest_dir <- file.path(proj_root, "resources", paste0(script_base, "_files"))
  
  # remove any empty placeholder so rename can succeed
  if (dir.exists(dest_dir))
    unlink(dest_dir, recursive = TRUE)
  
  # now move it wholesale
  if (dir.exists(orig_dir)) {
    file.rename(orig_dir, dest_dir)
  }
  
  # fix up any straggling paths in the .md
  md_path  <- file.path(proj_root, output_md)
  md_lines <- readLines(md_path)
  
  # 1) strip out any absolute proj_root prefix
  md_lines <- gsub(
    pattern     = paste0(proj_root, "/"),
    replacement = "",
    x           = md_lines,
    fixed       = TRUE
  )
  
  # 2) rewrite the remainder into resources/…
  md_lines <- gsub(
    pattern     = paste0(script_base, "_files/figure-gfm/"),
    replacement = fig_rel,
    x           = md_lines,
    fixed       = TRUE
  )
  
  writeLines(md_lines, md_path)
}

# remove leftover HTML files
htmls <- list.files(proj_root, "\\.html$", full.names = TRUE)
if (length(htmls)) file.remove(htmls)
