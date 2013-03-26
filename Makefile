CC=mpicc
CFLAGS=
LDFLAGS=

all: integrate

run: all
	mpiexec -n 8 integrate 0 10 1133

run_sin: all
	mpiexec -n 4 integrate 0 3.1416 100

clean:
	rm -rf integrate

integrate: integrate.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@
