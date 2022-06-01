CFLAGS = -Wall -g -O4
LDLIBS = -lm

OBJS = gen_word libwave.a

CC = gcc

all: $(OBJS)

gen_word: main_genword.o wave_io.o cmdline.o libwave.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

libwave.a: parameters.o mf3.o vf3.o prng.o keccak.o
	ar rc $@ $^
	ranlib $@


clean:
	- /bin/rm *.o

veryclean: clean
	- /bin/rm $(OBJS)

depend:
	- makedepend -Y *.c 2> /dev/null

# DO NOT DELETE

main_genword.o: wave.h cmdline.h wave_io.h vf3.o
cmdline.o: cmdline.h
mf3.o: wave.h
parameters.o: wave.h params.h
prng.o: wave.h
vf3.o: wave.h
wave_io.o: wave.h
