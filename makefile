all: miev0simple.so
	f2py-2.7 -c -m miev0simple main.f90 MIEV0.f ErrPack.f bhmie.f
