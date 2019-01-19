AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = build \
			 manifest.py

if PETSC_OPT
MAYBE_SOLVER = solver
MAYBE_IFLUID = iFluid
MAYBE_CRYSTAL = crystal
MAYBE_MELTING = melting
MAYBE_SUBSURF = subsurf
if CGAL_ARMA_OPT
MAYBE_COLLISION = Collision
if FOLDING_OPT
MAYBE_FOLDING = folding
endif
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


SUBDIRS = lib \
		  cFluid \
		  $(MAYBE_SOLVER) \
		  $(MAYBE_IFLUID) \
		  $(MAYBE_CRYSTAL) \
		  $(MAYBE_MELTING) \
		  $(MAYBE_SUBSURF) \
		  $(MAYBE_COLLISION) \
		  $(MAYBE_FOLDING) \
		  $(MAYBE_AIRFOIL) \
		  $(MAYBE_PARACHUTE) \
		  $(MAYBE_EX2D) \
		  $(MAYBE_EX3D) \
		  tests




.PHONY: tar

tar: create-tar

create-tar: manifest.py
	@ME=`basename $$PWD` && cp manifest.py ../. && \
		cd .. && ./manifest.py $$ME && rm -f manifest.py && \
		tar -cf "`date +FronTier++.%m_%d_%y.tar.gz`" -T Manifest.txt && \
		mv `date +FronTier++.%m_%d_%y.tar.gz` Manifest.txt $$ME/.

clean-local:
	rm -f Manifest.txt


