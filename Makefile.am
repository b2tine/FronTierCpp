AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I autoconf-archive/m4

if EX2D_OPT
	MAYBE_EX2D = example2d
endif

if EX3D_OPT
	MAYBE_EX3D = example3d
endif


SUBDIRS = lib tests cFluid $(MAYBE_EX2D) $(MAYBE_EX3D)
