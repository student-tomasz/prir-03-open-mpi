CC=mpicc
CFLAGS=
LDFLAGS=

all: integrate

run: all
	mpiexec -n 7 integrate 0 5 1133

clean:
	rm -rf integrate

integrate: integrate.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@
