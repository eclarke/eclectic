#!/usr/bin/env Rscript
#
# Rebuilds the package docs after a commit.
#
# Move this file after modifications to:
#   .git/hooks/post-commit
#
suppressMessages(stopifnot(require(git2r)))

here <- repository('.')
stopifnot(git2r::head(here)@name == "gh-pages")
staticdocs::build_site(site_path='.', launch = FALSE)
