
all:	ising_T01.root \
	ising_T24.root



ising_T24.root:
	root -b -q -l 'gen_ntuple.C("$@", 2.4, 0.0)'

ising_T01.root:
	root -b -q -l 'gen_ntuple.C("$@", 0.1, 0.0)'
