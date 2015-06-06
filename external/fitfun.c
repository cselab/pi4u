/*
 *  fitfun.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 600
#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <ftw.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <torc.h>
#include "spawner.c"


int dbg_display = 1;
static pthread_mutex_t fork_mutex = PTHREAD_MUTEX_INITIALIZER;

char bindir[256];
char datdir[256];
int input_method;		// 0: argv, 1: text file, 2: binary file
int output_method;		// 1: text file, 2: binary file 

//char basetmpdir[256];		// empty: current 
//int tmpdir1;			// 0: keep, 1: rmrf 
//int tmpdir2;			// 0: same, 1: unique

//char runscript[256];		// what to run

void sighandler(int signo)
{
	printf("spawner: SIGNAL %d was received!\n", signo); fflush(0);
}

double fitfun(double *TP, int   n, void *output, int *info)
{
	double res;
	int i;
	int me = getpid();	/* spanwer_id : worker_id */
	int rf;
	char line[1024];
	char *largv[64];
	char taskname[256];
	double t0, t1;

#if 0
	sprintf(taskname, "tmpdir.%d.%d", getpid(), 0);
#else
	sprintf(taskname, "tmpdir.%d.%d.%d.%d", info[0], info[1], info[2], info[3]);
#endif

	if (dbg_display) {
		printf("spanwer(%d): running task %s with params (", me, taskname);

		for (i = 0; i < n-1; i++)
			printf("%.16lf,", TP[i]);

		printf("%.16lf)\n", TP[i]);

		fflush(0);
	}

	/* create a (temporary) directory */
	mkdir(taskname, S_IRWXU);

	t0 = my_gettime();
#if 0
	rf = fork();
#else
	while (pthread_mutex_trylock(&fork_mutex) == EBUSY)
		usleep(500*1000);

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
		else {
			printf("give up\n");
			rmrf(taskname);
			return -1.0e6;
		}
	}
#endif

	if (rf == 0) {
		chdir(taskname);	/* enter to the new directory */

		/* 1. PREPROCESSING PHASE - APPLICATION SPECIFIC */
		/* copy some application specific input files*/
		copy_from_dir("../bin");
		//copy_from_dir("../dat");

#if 0
		/* pass/propagate the parameters to the running script */
		for (i = 0; i < n; i++) {
			char tmp[64];
			sprintf(tmp, " %.10lf", TP[i]);
			strcat(line, tmp);
                }
#else
		/* write the parameters into a input file */
		FILE *finp = fopen("params.dat", "w");
		for (i = 0; i < n; i++)
			fprintf(finp, "%.16lf\n", TP[i]);
		fclose(finp);
#endif


		/* 2. EXECUTION OF THE FUNCTION EVALUATION (SIMULATION) SOFTWARE */
		strcpy(line, "");
		strcat(line, "./run.sh");
		parse(line, largv);	/* prepare argv */

#if 1
		int fd = open("output", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);

		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);	// fd no longer needed
#endif
		res = execvp(*largv, largv);
		if (res < 0) exit(1);
	}
	pthread_mutex_unlock(&fork_mutex);

	waitpid(rf, NULL, 0);
	sync();
	t1 = my_gettime();
		
	if (dbg_display)
		printf("spanwer(%d): simcode_time=%lf secs\n", me, t1-t0);fflush(0);

	/* 3. POSTPROCESSING PHASE - APPLICATION SPECIFIC */
	/* read the result of the simulation - from "fitness" */
	FILE *fin;
	double gval = -1.0e6;
	char fitname[256];
	sprintf(fitname, "./%s/fitness", taskname);
	fin = fopen(fitname, "r");
	if (fin == NULL) {
		printf("spawner(%d): error with fopen for %s\n", me, fitname); fflush(0);
		res = gval;
	} else {
		fscanf(fin, "%lf", &gval);	// FIT = gval
		fclose(fin);
		if (dbg_display)
			printf("spanwer(%d): Result for task %s: %.16lf\n", me, taskname, gval); fflush(0);

		res = gval;
	}

	/* remove the temporary directory */
	if (1)
		rmrf(taskname);

	return res;
}
