.DELETE_ON_ERROR:
test: 1h.btwc.test 1c.btc2 1c.btc3 # CRN1.btc2
%.test:
	make $*
	diff $* ref/$*
%.btwc: genice_bondtwist/formats/bondtwist.py Makefile
	genice $* -r 2 2 2 -f bondtwist > $@
%.btwc.yap: genice_bondtwist/formats/bondtwist.py Makefile
	genice $* -r 2 2 2 -f bondtwist[yaplot] > $@
%.btwc.svg: genice_bondtwist/formats/bondtwist.py Makefile
	genice $* -r 2 2 2 -f bondtwist[svg] > $@
%.btc2: genice_bondtwist/formats/bondtwist2.py Makefile
	genice $* -r 2 2 2 -f bondtwist2 > $@
%.btc3: genice_bondtwist/formats/bondtwist3.py Makefile
	genice $* -r 2 2 2 -f bondtwist3 > $@
%.btc2.yap: genice_bondtwist/formats/bondtwist2.py Makefile
	genice $* -r 2 2 2 -f bondtwist2[yaplot] > $@
%.btc3.yap: genice_bondtwist/formats/bondtwist3.py Makefile
	genice $* -r 2 2 2 -f bondtwist3[yaplot] > $@
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
