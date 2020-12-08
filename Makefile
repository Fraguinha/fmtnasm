all: fmtnasm

.PHONY: clean

fmtnasm.ml : fmtnasm.mll
	ocamllex fmtnasm.mll

fmtnasm : fmtnasm.ml
	ocamlopt fmtnasm.ml -o fmtnasm

clean :
	rm -rf fmtnasm *.ml *.cmi *.cmx *.o
