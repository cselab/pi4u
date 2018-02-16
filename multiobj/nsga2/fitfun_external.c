#include <fcntl.h>
#include <ftw.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include "spawner.h"
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */

static pthread_mutex_t fork_mutex = PTHREAD_MUTEX_INITIALIZER;

#define PARAMSFILE	"params.txt"
#define OUTPUTFILE	"fitness.txt"
#define CONSTRFILE	"constr_viol.txt"
#define REMOVEDIRS	0

// REMOVEDIRS=1 removes the directories..., =0 keeps them

void fitfun(double /*const*/ *x, int N, void *output, int *winfo, double *result, int ncon, double *constraints)
{
	// 0.
	int n  = N;
	int me = getpid();
	char line[1024];
	char *largv[64];
	char taskname[256];
	double t0, t1;
	double res1, res2;

	int gen, chain, step, task;
	gen = winfo[0]; chain = winfo[1]; step = winfo[2]; task = winfo[3];

	// 1.
	//sprintf(taskname, "tmpdir.%d.%d", getpid(), torc_worker_id());

	int job_try = 0;

job_restart:

	sprintf(taskname, "tmpdir.%d.%d.%d.%d", gen, chain, step, task);
	mkdir(taskname, S_IRWXU | S_IRWXG | S_IRWXO);	

	{
	int i;
	printf("spanwer(%d): running task %s with params (", me, taskname);

	for (i = 0; i < n-1; i++) printf("%.6lf,", x[i]);
	printf("%.6lf)\n", x[i]);

	fflush(0);
	}


#if 1
	t0 = my_gettime();
	while (pthread_mutex_trylock(&fork_mutex) == EBUSY) usleep(500*1000);
#endif

	// 2.
#if 1
	int rf = fork();
	if (rf < 0)
	{
		printf("spanwer(%d): fork failed!!!!\n", me);
		perror("fork");
		fflush(0);
	}
#else
	int tries = 0;
retry:
	rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!! [try = %d]\n", me, tries); fflush(0);
		tries++;
		if (tries < 10) {
			sleep(1);
			goto retry;
		}
	}
#endif

	if (rf == 0)
	{
		// 2.1 enter tmp folder
		chdir(taskname);

		// 2.2 copy necessary stuff
		int status = copy_from_dir("../bin");
		if (status!=0)
		{
			printf("Error in copy from dir\n");
			abort();
		}

		// 1.3 write input parametes to the simulation's input file
		FILE *finp = fopen(PARAMSFILE, "w");
		int i;
		for (i = 0; i < n; i++) fprintf(finp, "%.16lf\n", x[i]);
		fclose(finp);

		/* 2. run simulation */
		sprintf(line, "./run.sh");
		parse(line, largv);

#if 1
		int fd = open("output", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);

		dup2(fd, 1);// make stdout go to file
		dup2(fd, 2);// make stderr go to file
		close(fd);	// fd no longer needed
#endif

		status = execvp(*largv, largv);
	}

	pthread_mutex_unlock(&fork_mutex);

	// 3.
#if 1
	int status;
	waitpid(rf, &status, 0);
#else
	if (rf > 0) {
		int status, wpid, ctr = 0;
		while ((wpid = waitpid(rf, &status, WNOHANG)) == 0) {
			ctr++;
			sleep(1);
			if (ctr == 120) {
				printf("XXXX: KILLING %d\n", rf); fflush(0);
				kill(rf, SIGKILL);
				break;
			}
		}
	}
#endif


	// 4a.
	double fitness1, fitness2;

	FILE *pFile;
	//double gval = -1, gres;
	char fitname[256];

	sprintf(fitname, "./%s/"OUTPUTFILE, taskname);
	pFile = fopen(fitname, "r");

	if (pFile == NULL)
	{
		printf("spawner(%d): error with fopen for %s\n", me, fitname); fflush(0);
		//fitness = -1e12; // peh: a bad value here, sse = +inf -> fitfun = -sse = -inf
		fitness1 = 1e12;
		fitness2 = 1e12;
	}
	else
	{
		fscanf(pFile, "%lf", &fitness1 );
		fscanf(pFile, "%lf", &fitness2 );

		fclose (pFile);
	}

	// 4b.
	if (ncon) {
		sprintf(fitname, "./%s/"CONSTRFILE, taskname);
		pFile = fopen(fitname, "r");

		if (pFile == NULL)
		{
			printf("spawner(%d): error with fopen for %s\n", me, fitname); fflush(0);
			for (int i = 0; i < ncon; ++i)
				constraints[i] = -1e12;
		}
		else
		{
			for (int i = 0; i < ncon; ++i)
				fscanf(pFile, "%lf", &constraints[i]);

			fclose (pFile);
		}
	}


	res1 = fitness1;
	res2 = fitness2;

	// inc_nfc();  /* increment function call counter*/

	// 5.
	if (REMOVEDIRS) // 0 for testing, 1 for production
		rmrf(taskname);


	t1 = my_gettime();

	{
	int i;
	char msg[1024], buf[64];
	sprintf(msg, "task(%d.%d.%d.%d):", gen, chain, step, task);
	for (i = 0; i < n; i++)
	{
		sprintf(buf, "%.6lf ", x[i]); 
		strcat(msg, buf);
	}
	sprintf(buf, " = %.6lf \t %.6lf", fitness1, fitness2); 
	strcat(msg, buf);

	sprintf(buf, " in %lf secs\n", t1-t0); 
	strcat(msg, buf);

	printf("%s", msg); fflush(0);
	}

	result[0] = res1;
	result[1] = res2;
	
	job_try++;
	if ((res1==1e12 || res2==1e12 || isnan(res1) || isnan(res2) || (t1-t0)<0.0) && (job_try < 5)) // i.e., if fitness.txt doesn't exist, or if the sim ended in less than 0 seconds real time, and retries have been less than 5.
		goto job_restart;


	return;
}

