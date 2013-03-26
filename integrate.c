#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef struct _range_t {
    double a,
           b;
    int n;
} range_t;

MPI_Datatype range_mpi_t;

void register_range_mpi_t();
void parse_args(const int argc, char *argv[], range_t **glob_rng);
void calculate_ranges(const range_t *glob_rng, const int procs_num, range_t ***proc_rngs);
int proc_is_master(const int proc_id);
double integrate(double (*f)(double), const double a, const double b, const int n);

double linear(double x)
{
    return 2.0*x;
}

int main(int argc, char *argv[])
{
    range_t *glob_rng;
    range_t **proc_rngs;
    range_t proc_rng;
    double rslt,
           part_rslt;

    int procs_num;
    int proc_id;
    int get_rngs_msg_tag = 112;
    int integrate_msg_tag = 111;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);

    register_range_mpi_t();

    if (proc_is_master(proc_id)) {
        parse_args(argc, argv, &glob_rng);
        calculate_ranges(glob_rng, procs_num, &proc_rngs);

        rslt = 0.0;
        for (proc_id = 1; proc_id < procs_num; proc_id++) {
            // printf("Wyslij zakres do procesu %d\n", proc_id);
            MPI_Send(proc_rngs[proc_id], 1, range_mpi_t, proc_id, get_rngs_msg_tag, MPI_COMM_WORLD);
            // printf("Czekaj na wynik od procesu %d\n", proc_id);
            MPI_Recv(&part_rslt, 1, MPI_DOUBLE, proc_id, integrate_msg_tag, MPI_COMM_WORLD, &status);
            // printf("Otrzymano wynik od procesu %d\n", proc_id);
            rslt += part_rslt;
        }
        printf("*** Calkowity wynik %f ***\n", rslt);

        for (proc_id = 1; proc_id < procs_num; proc_id++) {
            free(proc_rngs[proc_id]);
        }
        free(proc_rngs);
        free(glob_rng);
    }
    else {
        // printf("Proces %d czeka na zakres\n", proc_id);
        MPI_Recv(&proc_rng, 1, range_mpi_t, 0, get_rngs_msg_tag, MPI_COMM_WORLD, &status);
        // printf("Proces %d otrzymal zakres\n", proc_id);
        // printf("Proces %d wysyla czesciowy wynik\n", proc_id);
        part_rslt = integrate(&sin, proc_rng.a, proc_rng.b, proc_rng.n);
        MPI_Send(&part_rslt, 1, MPI_DOUBLE, 0, integrate_msg_tag, MPI_COMM_WORLD);
        printf("Proces %d wyslal czesciowy wynik %lf\n", proc_id, part_rslt);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}

void register_range_mpi_t()
{
    const int blocks_count = 3;
    const int blocks_lengths[blocks_count] = {1, 1, 1};
    MPI_Aint offsets[3] = {0, sizeof(double), 2*sizeof(double)};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(blocks_count, blocks_lengths, offsets, types, &range_mpi_t);
    MPI_Type_commit(&range_mpi_t);
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

void calculate_ranges(const range_t *glob_rng, const int procs_num, range_t ***proc_rngs)
{
    int i;
    range_t *proc_rng;

    int jobs_num = procs_num - 1;
    int points_per_job = round((double)glob_rng->n / jobs_num);
    double h = (glob_rng->b - glob_rng->a) / glob_rng-> n;

    *proc_rngs = malloc(sizeof(**proc_rngs) * procs_num);
    for (i = 1; i < procs_num; i++) {
        (*proc_rngs)[i] = malloc(sizeof(***proc_rngs));
    }

    for (i = 1; i < procs_num; i++) {
        proc_rng = (*proc_rngs)[i];
        proc_rng->n = points_per_job;
        if (i == jobs_num) {
            proc_rng->n = glob_rng->n - points_per_job * (jobs_num-1);
        }
    }

    for (i = 1; i < procs_num; i++) {
        proc_rng = (*proc_rngs)[i];
        proc_rng->a = glob_rng->a + h * points_per_job * (i-1);
        proc_rng->b = proc_rng->a + h * proc_rng->n;
    }

    // for (i = 1; i < procs_num; i++) {
    //     proc_rng = (*proc_rngs)[i];
    //     printf("Proces %d ma liczyc od %lf do %lf z rozdzielczoscia %d\n", i, proc_rng->a, proc_rng->b, proc_rng->n);
    // }
}

int proc_is_master(const int proc_id)
{
    return proc_id == 0;
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
