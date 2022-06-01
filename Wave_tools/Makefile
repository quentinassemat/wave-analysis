CFLAGS = -Wall -g -O4
LDLIBS = -lm

OBJS = analysis libwave.a 

CC = gcc

all: $(OBJS)

analysis: main_analysis.o wave_io.o cmdline.o libwave.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

libwave.a: analysis.o mf3.o vf3.o parameters.o keys.o
	ar rc $@ $^
	ranlib $@

clean:
	- /bin/rm *.o

veryclean: clean
	- /bin/rm $(OBJS)

depend:
	- makedepend -Y *.c 2> /dev/null

# DO NOT DELETE

analysis.o:  wave.h cmdline.h wave_io.h
cmdline.o: cmdline.h
mf3.o: wave.h
keys.o: wave.h
parameters.o: wave.h params.h
vf3.o: wave.h
wave_io.o: wave.h
