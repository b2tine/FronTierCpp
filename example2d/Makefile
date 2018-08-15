CXX = g++ -g

incs = -I../include/ -I../include/front -I../include/intfc -I../include/util
libincs = -L../libs




example00: example00.cpp FT_libs
	$(CXX) $<  $(incs) $(libincs) -lFronTier -o $@



FT_libs: ../libs/libFronTier.a

clean:
	$(RM) example00
	$(RM) *.o
