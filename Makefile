
docs:
	Rscript -e "devtools::document(roclets=c('rd', 'collate', 'namespace', 'vignette'))"

gh-pages:
	git checkout gh-pages
	git merge master -X theirs -m "merge master"

site:docs gh-pages
	Rscript -e "staticdocs::build_site(site_path='.', launch=FALSE)"
	git commit -am 'updated docs'
	git checkout master

update: site
	git push
	git checkout gh-pages
	git push
	git checkout master
