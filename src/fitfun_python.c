#include <stdio.h>
#include <pthread.h>

#include <Python.h>

//#if PY_MAJOR_VERSION >= 3
//#error "3"
//#else
//#error "2"
//#endif

#if PY_MAJOR_VERSION >= 3
#define PyInt_FromLong PyLong_FromLong
#define PyString_FromString PyBytes_FromString
#endif



#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

static PyObject *myModule, *myFunction;
pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;  /* serializes access to the Python interpreter */
static int dbg_display = 0;

#if PY_MAJOR_VERSION >= 3
void *
#else
void 
#endif
fitfun_initialize(char *fitfun_name)
{
	printf("xxxxxxxxxxxxx %s\n", fitfun_name);

	Py_Initialize();
	PyEval_InitThreads();

	PyGILState_STATE state = PyGILState_Ensure();

        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append(\".\")");

#if PY_MAJOR_VERSION >= 3
//	PyObject* myModuleString = PyUnicode_DecodeFSDefault((char *)"fitfun3");
	PyObject* myModuleString = PyUnicode_DecodeFSDefault(fitfun_name);
#else
//	PyObject* myModuleString = PyString_FromString((char*)"fitfun3");
	PyObject* myModuleString = PyString_FromString(fitfun_name);
#endif
	printf("myModuleString = %p\n", myModuleString);

        myModule = PyImport_Import(myModuleString);
        printf("myModule = %p\n", myModule);

	if (myModule == NULL) {
		printf("myModule = NULL\n"); 
		exit(1);
	}

//	myFunction = PyObject_GetAttrString(myModule,(char*)"fitfun0");
	myFunction = PyObject_GetAttrString(myModule,fitfun_name);
	printf("myFunction = %p\n", myFunction);

//	http://realgonegeek.blogspot.ch/2013/08/how-to-pass-c-array-to-python-solution.html
	import_array();

	PyGILState_Release(state);

#if PY_MAJOR_VERSION >= 3
	return 0;
#endif
}


void fitfun_finalize()
{
	return;
//	PyGILState_STATE state = PyGILState_Ensure();
//        Py_DECREF (myModule);
//        Py_DECREF (myFunction);

//	PyGILState_Release(state);
        Py_Finalize();
}

double fitfun(double *TP, long n, void *output, int *info)
{
	double res;
        int i;
        int me = getpid();	/* spanwer_id : worker_id */
//	int rf;
        char taskname[256];

	sprintf(taskname, "tmpdir.%d.%d.%d.%d", info[0], info[1], info[2], info[3]);

        if (dbg_display) {
                printf("spanwer(%d): running task %s with params (", me, taskname);
                for (i = 0; i < n-1; i++)
                        printf("%.16lf,", TP[i]);
                printf("%.16lf)\n", TP[i]);
                fflush(0);
        }
        npy_intp ld[1];
        ld[0] = n;

        pthread_mutex_lock(&m);
	PyGILState_STATE state = PyGILState_Ensure();

        PyArrayObject *py_array = (PyArrayObject *) PyArray_SimpleNewFromData(1, ld, NPY_DOUBLE, TP);

        PyObject* args = PyTuple_New(2);

#if VERBOSE
        printf("myFunction = %p\n", myFunction);
	printf("args = %p\n", args); fflush(0);
#endif

        PyTuple_SetItem(args, 0, (PyObject *)py_array);
        PyTuple_SetItem(args, 1, PyInt_FromLong(n));


        PyObject* myResult = PyObject_CallObject(myFunction, args);

        double result = PyFloat_AsDouble(myResult);

        if (dbg_display) {
                printf("result = %lf\n", result);
        }

        memcpy(&res, &result, 1*sizeof(double));
        Py_DECREF (py_array);
	PyGILState_Release( state );

        pthread_mutex_unlock(&m);
        return res;
}
