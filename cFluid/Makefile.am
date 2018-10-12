AM_CPPFLAGS = -I$(abs_top_srcdir)/include \
			  -I$(abs_top_srcdir)/include/front \
			  -I$(abs_top_srcdir)/include/intfc \
			  -I$(abs_top_srcdir)/include/util \
			  @HDF_CFLAGS@

AM_LDFLAGS = -L$(abs_top_builddir)/lib -lfront -lintfc -lutil \
			 @HDF_LIBS@

noinst_PROGRAMS = cFluid

cFluid_SOURCES = cFluid.h \
				 cFluid.cpp \
				 cFbase.cpp \
				 cFcartsn.cpp \
				 cFeos.cpp \
				 cFinit.cpp \
				 cFriem.cpp \
				 cFsub.cpp \
				 cFtvd.cpp \
				 cFweno.cpp

