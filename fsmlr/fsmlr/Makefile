COMPILER = g++
CFLAGS = -O3 -Wno-deprecated -static
RM = del

# object files
OBJS = 	fsmlr-main.o\
		FSMLR1.o
		
EXE = fsmlr-main.exe

########################
test-svector: $(OBJS)
	$(COMPILER) -o fsmlr-main.exe $(OBJS) $(CFLAGS)

.cpp.o:
	$(COMPILER) -c $< $(CFLAGS)
	
clean:
	$(RM) $(OBJS)
	$(RM) $(EXE)
	
run:
	fsmlr-main.exe -con fragment.con -tmp fragment.tmp -set Abraham.set -cfl Abraham1tvp.CFL -nprop 1 -dsc descr.dsc
	
	
