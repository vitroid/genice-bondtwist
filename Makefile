.DELETE_ON_ERROR:
test: T2.btwi
T2.btwi: genice_bondtwist/formats/bondtwist.py Makefile
	genice T2 -f bondtwist > $@
check:
	./setup.py check
install:
	./setup.py install
pypi: check
	./setup.py sdist bdist_wheel upload
clean:
	-rm $(ALL) *~ */*~ *svg
	-rm -rf build dist *.egg-info
	-find . -name __pycache__ | xargs rm -rf
