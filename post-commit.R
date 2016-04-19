#!/usr/bin/env Rscript
#
# Rebuilds the package docs after a commit.
#
# Move this file after modifications to:
#   .git/hooks/post-commit
#
suppressMessages(stopifnot(require(git2r)))

here <- repository('.')
checkout(here, branch='gh-pages')
message("Updating docs...")
git2r::merge(here, "master", strategy.option="theirs")

suppressMessages(staticdocs::build_site(site_path='.', launch = FALSE))
message("Static docs built; pushing...")
push(here)
message("Static docs updated.")
