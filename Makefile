CC=mpicc
CFLAGS=
LDFLAGS=

all: integrate

run: all
	mpiexec -n 7 integrate 0 5 1133

run_sin: all
	mpiexec -n 4 integrate 0 3.1416 123

clean:
	rm -rf integrate

integrate: integrate.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@
