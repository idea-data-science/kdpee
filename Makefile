check:
	Rscript -e 'devtools::check()'

clean:
	rm -f src/kdpee.o src/kdpee.so
