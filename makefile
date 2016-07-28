all: miev0simple.so

miev0simple.so: main.f90 MIEV0.f ErrPack.f bhmie.f
	f2py-2.7 -c -m miev0simple main.f90 MIEV0.f ErrPack.f bhmie.f

#all: miev0simple.so
#	python compile.py
