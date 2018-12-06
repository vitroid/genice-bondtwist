.DELETE_ON_ERROR:
test: 1h.btwc.test
%.test:
	make $*
	diff $* ref/$*
1h.btwc: genice_bondtwist/formats/bondtwist.py Makefile
	genice 1h -r 2 2 2 -f bondtwist > $@
check:
	./setup.py check
install:
	./setup.py install
pypi: check
	./setup.py sdist bdist_wheel upload
clean:
	-rm $(ALL) *~ */*~ *.btwc
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
