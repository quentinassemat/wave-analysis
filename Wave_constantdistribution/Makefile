CFLAGS = -Wall -g -O4
LDLIBS = -lm

OBJS = keygen sign verif showpr libwave.a

CC = gcc

all: $(OBJS)

keygen: main_keygen.o wave_io.o cmdline.o libwave.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

sign: main_sign.o wave_io.o cmdline.o libwave.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

verif: main_verif.o wave_io.o cmdline.o libwave.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

showpr: showpr.o cmdline.o libwave.a
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

libwave.a: keygen.o parameters.o sign.o verif.o mf3.o vf3.o prng.o keccak.o
	ar rc $@ $^
	ranlib $@

clean:
	- /bin/rm *.o

veryclean: clean
	- /bin/rm $(OBJS)

depend:
	- makedepend -Y *.c 2> /dev/null

# DO NOT DELETE

cmdline.o: cmdline.h
keygen.o: wave.h
main_keygen.o: cmdline.h wave.h wave_io.h
main_sign.o: cmdline.h wave.h wave_io.h
main_verif.o: cmdline.h wave.h wave_io.h
mf3.o: wave.h
parameters.o: wave.h params.h
prng.o: wave.h
showpr.o: cmdline.h wave.h
sign.o: wave.h
verif.o: wave.h
vf3.o: wave.h
wave_io.o: wave.h
