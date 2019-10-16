all: sync

sync:
	@rsync -avu --include="report/***" --include="cache/***" --exclude="*" --exclude="*.Rds" atlas:phd/research/articletwo/ .
	#@rsync -avu --include="report/***" --exclude="*" atlas:phd/research/articletwo/ .

run:
	@ssh atlas grun.py -q highmem.q -n ppca_1008 -c "$HOME/.bin/rmd.sh articletwo.Rmd $HOME/phd/research/articletwo"

clean:
	@rm -Rf ./rds ./cache

$(MAIN).tar.gz: $(FIGURES) $(SOURCES)
	tar -czf $(MAIN).tar.gz $(FIGURES) $(SOURCES)
