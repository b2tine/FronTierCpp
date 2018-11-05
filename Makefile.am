AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = build

if PETSC_OPT
MAYBE_IFLUID = iFluid
MAYBE_CRYSTAL = crystal
#MAYBE_SUBSURF = subsurf
if CGAL_OPT
MAYBE_AIRFOIL = airfoil
MAYBE_PARACHUTE = parachute
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
		  $(MAYBE_PARACHUTE) \
		  $(MAYBE_CRYSTAL) \
		  $(MAYBE_SUBSURF) \
		  $(MAYBE_EX2D) \
		  $(MAYBE_EX3D)

clean-local:
	-rm -rf tags
	-rm -rf TAGS
	-rm -rf CTAGS
	-rm -rf Makefile

