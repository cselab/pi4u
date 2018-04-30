#include <fitfun.h>
#include "spawner.h"
#include "loglike_theta.h"

static pthread_mutex_t fork_mutex = PTHREAD_MUTEX_INITIALIZER;
static int flag[4096];  // MAX_WORKERS

#define REMOVEDIRS  1
#define PREFIX      "."  // "/scratch"
#define FAIL        -1e12
#define BUFLEN      1024



//#include "engine_tmcmc.h"
//extern data_t data;

void loglike_theta_initialize() {
}





void loglike_finalize() {
}





double loglike_theta(double *x, int n, void *output, int *winfo) {
    char workdir[BUFLEN], bindir[BUFLEN];
    double t, loglike;
    int i;

    // make tmp directory
    char cwd[BUFLEN]; getcwd(cwd, BUFLEN);
    sprintf(bindir, "%s/model", cwd);
    sprintf(workdir, PREFIX"/tmpdir.%d.%d.%d.%d",
            winfo[0], winfo[1], winfo[2], winfo[3]);
    mkdir(workdir, S_IRWXU);

    // info
    printf("Spawner(%d): running in %s with params", getpid(), workdir);
    for (i = 0; i < n; i++) printf(" %.6lf", x[i]);
    printf("\n");
    fflush(0);

    // fork simulation
    t = my_gettime();
    while (pthread_mutex_trylock(&fork_mutex) == EBUSY) usleep(500000);

    int rf = fork();
    if (rf < 0) {
        printf("Fork failed\n"); fflush(0);
    } else if (rf == 0) {



		chdir(workdir);

        // copy necessary stuff
        if (copy_from_dir(bindir) != 0) {
            printf("Error in copy from dir %s\n", bindir);
            abort();
        }


        // write parametes to the simulation's input file
        FILE *finp = fopen("params.txt", "w");
        for (i = 0; i < n; i++) fprintf(finp, "%.16lf\n", x[i]);
        fclose(finp);

        // run simulation
        char line[BUFLEN], *largv[64];
        sprintf(line, "sh ./doall.sh");
        parse(line, largv);
        int fd = open("output.txt", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        dup2(fd, 1); dup2(fd, 2);  // make stdout and stderr go to file
        close(fd);

		execvp(*largv, largv);

		printf("This point must not be reached!\n");
		exit(1);
    }

	pthread_mutex_unlock(&fork_mutex);


    // wait for fork
    int status; waitpid(rf, &status, 0); sync();

    // read results
    char llfile[BUFLEN]; sprintf(llfile, "%s/loglike.txt", workdir);
    FILE * pFile = fopen(llfile, "r");
    if (pFile == NULL) {
        loglike = FAIL;
    } else {
        while (!feof(pFile)) fscanf(pFile, "%lf", &loglike);
        fclose (pFile);
    }
    if (isnan(loglike) || isinf(loglike)) loglike = FAIL;

    // cleanup
    if (REMOVEDIRS && flag[torc_worker_id()] == 0) rmrf(workdir);

    // info
    printf("task(%d.%d.%d.%d): ", winfo[0], winfo[1], winfo[2], winfo[3]);
    for (i = 0; i < n; i++) printf("%.6lf ", x[i]);
    printf(" = %.6lf in %lf secs\n", loglike, my_gettime()-t);
    fflush(0);



    return loglike;
}





double loglike_(double *x, int n, void *output, int *winfo){
	return loglike_theta( x, n, output, winfo );
}
