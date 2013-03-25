#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct _range_t {
    double a,
           b;
    int n;
} range_t;

void parse_args(const int argc, char *argv[], range_t **glob_rng);
void calculate_range(const range_t *glob_rng, const int procs_num, const int proc_id, range_t **proc_rng);
int proc_is_master(const int proc_id);
int proc_is_last(const int proc_id);
double integrate(double (*f)(double), const double a, const double b, const int n);

double linear(double x)
{
    return 2.0*x;
}

int main(int argc, char *argv[])
{
    range_t *glob_rng,
            *proc_rng;
    double rslt,
           part_rslt;

    int procs_num;
    int proc_id;
    int job_tag = 111;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    parse_args(argc, argv, &glob_rng);
    calculate_range(glob_rng, procs_num, proc_id, &proc_rng);
    part_rslt = integrate(&linear, proc_rng->a, proc_rng->b, proc_rng->n);

    /* licz part_rslt */
    printf("Proces %d obliczyl czesciowy wynik %lf\n", proc_id, part_rslt);

    if (proc_is_master(proc_id)) {
        rslt = part_rslt;
        for (proc_id = 1; proc_id < procs_num; proc_id++) {
            // printf("Czekaj na wiadomosc od procesu %d\n", proc_id);
            MPI_Recv(&part_rslt, 1, MPI_DOUBLE, proc_id, job_tag, MPI_COMM_WORLD, &status);
            rslt += part_rslt;
        }
        printf("*** Calkowity wynik %f ***\n", rslt);
    }
    else {
        // printf("Wyslij wiadomosc z procesu %d\n", proc_id);
        MPI_Send(&part_rslt, 1, MPI_DOUBLE, 0, job_tag, MPI_COMM_WORLD);
    }

    free(glob_rng);
    free(proc_rng);

    MPI_Finalize();
    return EXIT_SUCCESS;
}

void parse_args(const int argc, char *argv[], range_t **glob_rng)
{
    if (argc != 4) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    *glob_rng = malloc(sizeof(**glob_rng));
    (*glob_rng)->a = atof(argv[1]);
    (*glob_rng)->b = atof(argv[2]);
    (*glob_rng)->n = atoi(argv[3]);
}

void calculate_range(const range_t *glob_rng, const int procs_num, const int proc_id, range_t **proc_rng)
{
    int points_per_proc = round((double)glob_rng->n / procs_num);
    double h = (glob_rng->b - glob_rng->a) / glob_rng-> n;

    *proc_rng = malloc(sizeof(**proc_rng));

    (*proc_rng)->n = points_per_proc;
    if (proc_is_last(proc_id)) {
        (*proc_rng)->n = glob_rng->n - points_per_proc * (procs_num-1);
    }

    (*proc_rng)->a = glob_rng->a + h * points_per_proc * proc_id;
    (*proc_rng)->b = (*proc_rng)->a + h * (*proc_rng)->n;

    // printf("Proces %d ma liczyc od %lf do %lf z rozdzielczoscia %d\n", proc_id, (*proc_rng)->a, (*proc_rng)->b, (*proc_rng)->n);
}

int proc_is_master(const int proc_id)
{
    return proc_id == 0;
}

int proc_is_last(const int proc_id)
{
    int procs_num;
    MPI_Comm_size(MPI_COMM_WORLD, &procs_num);

    return proc_id == procs_num-1;
}

double integrate(double (*f)(double), const double a, const double b, const int n)
{
    double rslt = 0.0;
    double h = (b - a) / n;

    int i;
    for (i = 0; i < n; i++) {
        rslt += f(a + h*i + h/2) * h;
    }

    return rslt;
}
