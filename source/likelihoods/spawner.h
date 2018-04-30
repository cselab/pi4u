#ifndef _SPAWNER_H_
#define _SPAWNER_H_

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

#include <pthread.h>
//#include <torc.h>


double my_gettime();
void parse(char *line, char **argv);
int execute_cmd(int me, char *largv[], char *dir);
int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf);
int rmrf(char *path);
int cp(const char *from, const char *to);
int copy_from_dir(char *name);
int copy_file(const char *dirname, const char *filename);

#endif
