AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

if PETSC_OPT
MAYBE_IFLUID = iFluid
MAYBE_CRYSTAL = crystal
if CGAL_OPT
MAYBE_AIRFOIL = airfoil
endif
endif

if EX2D_OPT
MAYBE_EX2D = example2d
endif

if EX3D_OPT
MAYBE_EX3D = example3d
endif


SUBDIRS = src \
		  lib \
		  tests \
		  cFluid \
		  $(MAYBE_IFLUID) \
		  $(MAYBE_AIRFOIL) \
		  $(MAYBE_CRYSTAL) \
		  $(MAYBE_EX2D) \
		  $(MAYBE_EX3D)


EXTRA_DIST = build
