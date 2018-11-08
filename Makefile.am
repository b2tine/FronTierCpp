AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = build

if PETSC_OPT
MAYBE_SOLVER = solver
MAYBE_IFLUID = iFluid
MAYBE_CRYSTAL = crystal
MAYBE_MELTING = melting
#MAYBE_SUBSURF = subsurf
if CGAL_OPT
MAYBE_AIRFOIL = airfoil
#MAYBE_PARACHUTE = parachute
endif
endif

if CGAL_OPT
MAYBE_COLLISION = Collision
endif

if EX2D_OPT
MAYBE_EX2D = example2d
endif

if EX3D_OPT
MAYBE_EX3D = example3d
endif


SUBDIRS = lib \
		  cFluid \
		  $(MAYBE_SOLVER) \
		  $(MAYBE_IFLUID) \
		  $(MAYBE_CRYSTAL) \
		  $(MAYBE_MELTING) \
		  $(MAYBE_SUBSURF) \
		  $(MAYBE_COLLISION) \
		  $(MAYBE_AIRFOIL) \
		  $(MAYBE_PARACHUTE) \
		  $(MAYBE_EX2D) \
		  $(MAYBE_EX3D) \
		  tests


TARFILELIST = \
$$ME/LICENSE \
$$ME/README \
$$ME/build \
$$ME/m4/*.m4 \
$$ME/configure.ac \
$$ME/Makefile.am \
$$ME/include/FronTier.h \
$$ME/src/Makefile.am \
$$ME/src/util/Makefile.am \
$$ME/src/util/*.c \
$$ME/src/util/*.cpp \
$$ME/include/util/*.h \
$$ME/src/intfc/Makefile.am \
$$ME/src/intfc/*.c \
$$ME/src/intfc/*.cpp \
$$ME/src/intfc/*.cuh  \
$$ME/src/intfc/*.cu \
$$ME/include/intfc/*.h  \
$$ME/src/front/Makefile.am \
$$ME/src/front/*.c \
$$ME/src/front/*.cpp \
$$ME/src/front/*.cuh \
$$ME/src/front/*.cu \
$$ME/include/front/*.h \
$$ME/include/front/*.hpp \
$$ME/lib/Makefile.am \
$$ME/lib/gmock/Makefile.am \
$$ME/lib/gmock/gmock/gmock.h \
$$ME/lib/gmock/gtest/gtest.h \
$$ME/lib/gmock/src/*.cc \
$$ME/tests/Makefile.am \
$$ME/tests/*.cpp \
$$ME/example2d/Makefile.am \
$$ME/example2d/*.cpp \
$$ME/example2d/README \
$$ME/example3d/Makefile.am \
$$ME/example3d/*.cpp \
$$ME/example3d/input/* \
$$ME/example3d/README \
$$ME/cFluid/Makefile.am \
$$ME/cFluid/*.h \
$$ME/cFluid/*.cpp \
$$ME/cFluid/in-* \
$$ME/cFluid/input-*/in-* \
$$ME/cFluid/mrun* \
$$ME/iFluid/Makefile.am \
$$ME/iFluid/*.h \
$$ME/iFluid/*.cpp \
$$ME/iFluid/input-*/in-* \
$$ME/iFluid/mrun* \
$$ME/solver/*.h \
$$ME/solver/*.cpp \
$$ME/rgbody/*.h \
$$ME/rgbody/*.cpp \
$$ME/rgbody/in-* \
$$ME/rgbody/input-*/in-* \
$$ME/rgbody/mrun* \
$$ME/Collision/Makefile.am \
$$ME/Collision/*.h \
$$ME/Collision/*.cpp \
$$ME/Collision/input/* \
$$ME/airfoil/Makefile.am \
$$ME/airfoil/*.h \
$$ME/airfoil/*.cpp \
$$ME/airfoil/*.cuh \
$$ME/airfoil/*.cu \
$$ME/airfoil/input-*/* \
$$ME/airfoil/mrun* \
$$ME/parachute/Makefile.am \
$$ME/parachute/*.cpp \
$$ME/parachute/input-*/* \
$$ME/parachute/mrun* \
$$ME/parachute2d/*.h \
$$ME/parachute2d/*.cpp \
$$ME/parachute2d/input-*/* \
$$ME/parachute2d/mrun* \
$$ME/crystal/Makefile.am \
$$ME/crystal/*.h \
$$ME/crystal/*.cpp \
$$ME/crystal/input-*/* \
$$ME/crystal/mrun* \
$$ME/crystal/README \
$$ME/subsurf/Makefile.am \
$$ME/subsurf/*.h \
$$ME/subsurf/*.cpp \
$$ME/subsurf/input-2d/* \
$$ME/melting/Makefile.am \
$$ME/melting/*.h \
$$ME/melting/*.cpp \
$$ME/melting/input/* \
$$ME/melting/mrun*

tar:
	cd ./ && ME=`basename $$PWD` && cd .. && tar -cf "`date +FronTier++.%m_%d_%y.tar`" ${TARFILELIST} && \
	gzip "`date +FronTier++.%m_%d_%y.tar`" && mv `date +FronTier++.%m_%d_%y.tar`.gz $$ME/.

