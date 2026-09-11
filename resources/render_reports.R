# Batch rendering of .R scripts to GitHub-friendly markdown reports

packages_needed <- c(
  "rmarkdown",
  "tools",
  "rprojroot",
  "knitr",
  "xfun"
)

missing_packages <- setdiff(
  packages_needed,
  rownames(installed.packages())
)

if (length(missing_packages)) {
  install.packages(missing_packages)
}

# Find project root
proj_root <- rprojroot::find_rstudio_root_file()

# List all .R scripts in /code
scripts <- list.files(
  file.path(proj_root, "code"),
  pattern = "\\.R$",
  full.names = TRUE
)

# Function executed in a fresh R session for each report
render_report <- function(script, proj_root, output_md) {
  
  knitr::opts_knit$set(
    root.dir = proj_root
  )
  
  rmarkdown::render(
    input         = script,
    output_file   = output_md,
    output_dir    = proj_root,
    knit_root_dir = proj_root,
    envir         = globalenv()
  )
}

for (script in scripts) {
  
  script_base <- tools::file_path_sans_ext(basename(script))
  output_md   <- paste0(script_base, ".md")
  
  tryCatch({
    
    # Render each report in its own clean R process
    xfun::Rscript_call(
      render_report,
      args = list(
        script    = script,
        proj_root = proj_root,
        output_md = output_md
      ),
      options = "--vanilla"
    )
    
  }, error = function(e) {
    
    message(
      "ERROR in ",
      basename(script),
      ": ",
      conditionMessage(e)
    )
    
    stop(e)
  })
  
  # Move generated figure directory from project root to resources/
  orig_dir <- file.path(
    proj_root,
    paste0(script_base, "_files")
  )
  
  dest_dir <- file.path(
    proj_root,
    "resources",
    paste0(script_base, "_files")
  )
  
  if (dir.exists(orig_dir)) {
    
    # Remove previous version so the entire new directory can replace it
    if (dir.exists(dest_dir)) {
      unlink(dest_dir, recursive = TRUE)
    }
    
    # Make sure resources/ itself exists
    dir.create(
      dirname(dest_dir),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # Move generated directory
    moved <- file.rename(orig_dir, dest_dir)
    
    if (!moved) {
      stop(
        "Could not move ",
        orig_dir,
        " to ",
        dest_dir
      )
    }
  }
  
  # Rewrite figure paths in generated markdown
  md_path <- file.path(proj_root, output_md)
  
  if (file.exists(md_path)) {
    
    md_lines <- readLines(md_path)
    
    # Strip any absolute project-root prefix
    md_lines <- gsub(
      pattern     = paste0(proj_root, "/"),
      replacement = "",
      x           = md_lines,
      fixed       = TRUE
    )
    
    # Point generated figure references into resources/
    md_lines <- gsub(
      pattern = paste0(
        script_base,
        "_files/figure-gfm/"
      ),
      replacement = paste0(
        "resources/",
        script_base,
        "_files/figure-gfm/"
      ),
      x = md_lines,
      fixed = TRUE
    )
    
    writeLines(md_lines, md_path)
  }
}

# Remove leftover HTML files
htmls <- list.files(
  proj_root,
  "\\.html$",
  full.names = TRUE
)

if (length(htmls)) {
  file.remove(htmls)
}