
double my_gettime();
void parse(char *line, char **argv);
int execute_cmd(int me, char *largv[], char *dir);
int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf);
int rmrf(char *path);
int cp(const char *from, const char *to);
int copy_from_dir(char *name);
int copy_file(const char *dirname, const char *filename);

