/*
 *  torc.c
 *  TORC_Lite
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */
#include <torc_internal.h>
#include <torc.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define f77fun  1

extern MPI_Comm comm_out;

void torc_waitall()
{
    _torc_block();
}

void torc_waitall2()
{
    _torc_block2();
}

void torc_waitall3()
{
    torc_t *self = _torc_self();

    _lock_acquire (&self->lock);
    --self->ndep;
    _lock_release (&self->lock);

    int remdeps;
    while (1)
    {
        _lock_acquire (&self->lock);
        remdeps = self->ndep;
        _lock_release (&self->lock);
        if (remdeps == 0)
            break;
        else {
            thread_sleep(0);
            //usleep(1*1000);
            //sched_yield();
        }
    }
}

int torc_scheduler_loop(int once)
{
    return _torc_scheduler_loop(once);
}

#define _initialize(rte)                              \
{                                                     \
    rte->homenode = rte->sourcenode = torc_node_id(); \
    _lock_init(&rte->lock);                           \
    rte->target_queue = rte->vp_id = -1;              \
    rte->work_id = -1;                                \
}

#ifdef TORC_STATS
static int invisible_flag = 0;
void torc_set_invisible(int flag)
{
    invisible_flag = flag;
}
#endif

void torc_task_detached(int queue, void (*work)(), int narg, ...)
{
    va_list ap;
    int i;
    torc_t * rte;

    rte = _torc_get_reused_desc();

    _initialize(rte);

    _torc_set_work_routine(rte, work);

    rte->narg       = narg;
    rte->rte_type   = 1;    /* external */
    rte->inter_node = 1;
#if 0
    rte->parent = self;
#else
    rte->parent     = NULL;
#endif
    rte->level      = 0;

    if (narg>MAX_TORC_ARGS) Error("narg > MAX_TORC_ARGS");

    va_start(ap, narg);

    for (i=0; i<narg; i++) {
        rte->quantity[i] = va_arg(ap, int);
        rte->dtype[i]    = va_arg(ap, MPI_Datatype);
        rte->btype[i]    = _torc_mpi2b_type(rte->dtype[i]);
        rte->callway[i]  = va_arg(ap, int);

        if ((rte->callway[i] == CALL_BY_COP)&&(rte->quantity[i] > 1)) rte->callway[i] = CALL_BY_COP2;

#if DBG
        printf("ARG %d : Q = %d, T = %d, C = %x O\n", i, rte->quantity[i], rte->dtype[i], rte->callway[i]); fflush(0);
#endif
    }

    for (i=0; i<narg; i++) {
        if (rte->quantity[i] == 0) {    // peh: 02.07.2015
            VIRT_ADDR dummy = va_arg (ap, VIRT_ADDR);
            continue;
        }
        if (rte->callway[i] == CALL_BY_COP) {
            int typesize;
            MPI_Type_size(rte->dtype[i], &typesize);
            switch (typesize) {
                case 4:
                    rte->localarg[i] = *va_arg (ap, INT32 *);
                    break;
                case 8:
                    rte->localarg[i] = *va_arg (ap, INT64 *);
                    break;
                default:
                    Error("typesize not 4 or 8!");
                    break;
            }
        }
        else if (rte->callway[i] == CALL_BY_COP2) {
            VIRT_ADDR addr = va_arg (ap, VIRT_ADDR);
            int typesize;
            void *pmem;
            MPI_Type_size(rte->dtype[i], &typesize);
            pmem = malloc(rte->quantity[i]*typesize);
            memcpy(pmem, (void *)addr, rte->quantity[i]*typesize);
            rte->localarg[i] = (INT64)pmem; //yyyyyy
        }
        else {
            rte->localarg[i] = va_arg (ap, VIRT_ADDR);    /* pointer (C: PTR, VAL) */
        }
    }

    if (queue == -1) {
        torc_to_rq_end(rte);
    }
    else {
        torc_to_lrq_end(queue, rte);
    }
}

void torc_task(int queue, void (*work)(), int narg, ...)
{
    va_list ap;
    int i;
    torc_t * rte;
    torc_t *self = _torc_self();

    /* Check if rte_init has been called */
    _lock_acquire(&self->lock);
    if (self->ndep == 0) self->ndep = 1;
    _lock_release(&self->lock);

    _torc_depadd(self, 1);

    rte = _torc_get_reused_desc();

    _initialize(rte);

    _torc_set_work_routine(rte, work);

    rte->narg       = narg;
    rte->rte_type   = 1;    /* external */
    rte->inter_node = 1;
    rte->parent     = self;
    rte->level      = self->level + 1;

#ifdef TORC_STATS
    if (!invisible_flag) created[self->vp_id]++;

    if (invisible_flag) rte->rte_type = 2;    /* invisible */
#endif

    if (narg>MAX_TORC_ARGS) Error("narg > MAX_TORC_ARGS");

    va_start(ap, narg);

    for (i=0; i<narg; i++) {
        rte->quantity[i] = va_arg (ap, int);
        rte->dtype[i]    = va_arg (ap, MPI_Datatype);
        rte->btype[i]    = _torc_mpi2b_type(rte->dtype[i]);
        rte->callway[i]  = va_arg (ap, int);

        if ((rte->callway[i] == CALL_BY_COP)&&(rte->quantity[i] > 1)) rte->callway[i] = CALL_BY_COP2;

#if DBG
        printf("ARG %d : Q = %d, T = %d, C = %x O\n", i, rte->quantity[i], rte->dtype[i], rte->callway[i]); fflush(0);
#endif
    }

    for (i=0; i<narg; i++) {
        if (rte->quantity[i] == 0) {    // peh: 02.07.2015
            VIRT_ADDR dummy = va_arg (ap, VIRT_ADDR);
            continue;
        }
        if (rte->callway[i] == CALL_BY_COP) {
            int typesize;
            MPI_Type_size(rte->dtype[i], &typesize);
            switch (typesize) {
                case 4:
                    rte->localarg[i] = *va_arg (ap, INT32 *);
                    break;
                case 8:
                    rte->localarg[i] = *va_arg (ap, INT64 *);
                    break;
                default:
                    Error("typesize not 4 or 8!");
                    break;
            }
        }
        else if (rte->callway[i] == CALL_BY_COP2) {
            VIRT_ADDR addr = va_arg (ap, VIRT_ADDR);
            int typesize;
            void *pmem;
            MPI_Type_size(rte->dtype[i], &typesize);
            pmem = malloc(rte->quantity[i]*typesize);
            memcpy(pmem, (void *)addr, rte->quantity[i]*typesize);
            rte->localarg[i] = (INT64)pmem; //yyyyyy
        }
        else {
            rte->localarg[i] = va_arg (ap, VIRT_ADDR);    /* pointer (C: PTR, VAL) */
        }
    }

    if (queue == -1) {
        torc_to_rq_end(rte);
    }
    else {
        torc_to_lrq_end(queue, rte);
    }
}

void torc_task_ex(int queue, int invisible, void (*work)(), int narg, ...)
{
    va_list ap;
    int i;
    torc_t * rte;
    torc_t *self = _torc_self();

    /* Check if rte_init has been called */
    _lock_acquire(&self->lock);
    if (self->ndep == 0) self->ndep = 1;
    _lock_release(&self->lock);

    _torc_depadd(self, 1);

    rte = _torc_get_reused_desc();

    _initialize(rte);

    _torc_set_work_routine(rte, work);

    rte->narg       = narg;
    rte->rte_type   = 1;    /* external */
    rte->inter_node = 1;
    rte->parent     = self;
    rte->level      = self->level + 1;

#ifdef TORC_STATS
    if (!invisible) created[self->vp_id]++;

    if (invisible) rte->rte_type = 2;    /* invisible */
#endif

    if (narg>MAX_TORC_ARGS) Error("narg > MAX_TORC_ARGS");

    va_start(ap, narg);

    for (i=0; i<narg; i++) {
        rte->quantity[i] = va_arg (ap, int);
        rte->dtype[i]    = va_arg (ap, MPI_Datatype);
        rte->btype[i]    = _torc_mpi2b_type(rte->dtype[i]);
        rte->callway[i]  = va_arg (ap, int);

        if ((rte->callway[i] == CALL_BY_COP)&&(rte->quantity[i] > 1)) rte->callway[i] = CALL_BY_COP2;

#if DBG
        printf("ARG %d : Q = %d, T = %d, C = %x O\n", i, rte->quantity[i], rte->dtype[i], rte->callway[i]); fflush(0);
#endif
    }

    for (i=0; i<narg; i++) {
        if (rte->quantity[i] == 0) {    // peh: 02.07.2015
            VIRT_ADDR dummy = va_arg (ap, VIRT_ADDR);
            continue;
        }
        if (rte->callway[i] == CALL_BY_COP) {
            int typesize;
            MPI_Type_size(rte->dtype[i], &typesize);
            switch (typesize) {
                case 4:
                    rte->localarg[i] = *va_arg (ap, INT32 *);
                    break;
                case 8:
                    rte->localarg[i] = *va_arg (ap, INT64 *);
                    break;
                default:
                    Error("typesize not 4 or 8!");
                    break;
            }
        }
        else if (rte->callway[i] == CALL_BY_COP2) {
            VIRT_ADDR addr = va_arg (ap, VIRT_ADDR);
            int typesize;
            void *pmem;
            MPI_Type_size(rte->dtype[i], &typesize);
            pmem = malloc(rte->quantity[i]*typesize);
            memcpy(pmem, (void *)addr, rte->quantity[i]*typesize);
            rte->localarg[i] = (INT64)pmem; //yyyyyy
        }
        else {
            rte->localarg[i] = va_arg (ap, VIRT_ADDR);    /* pointer (C: PTR, VAL) */
        }
    }

    if (queue == -1) {
        torc_to_rq_end(rte);
    }
    else {
        torc_to_lrq_end(queue, rte);
    }
}

void torc_task_direct(int queue, void (*work)(), int narg, ...)
{
    va_list ap;
    int i;
    torc_t * rte;
    torc_t *self = _torc_self();

    /* Check if rte_init has been called */
    _lock_acquire(&self->lock);
    if (self->ndep == 0) self->ndep = 1;
    _lock_release(&self->lock);

    _torc_depadd(self, 1);

    rte = _torc_get_reused_desc();

    _initialize(rte);

    _torc_set_work_routine(rte, work);

    rte->narg       = narg;
    rte->rte_type   = 20;    /* external - direct execution */
    rte->inter_node = 1;
    rte->parent     = self;
    rte->level      = self->level + 1;

    if (narg>MAX_TORC_ARGS) Error("narg > MAX_TORC_ARGS");

    va_start(ap, narg);

    for (i=0; i<narg; i++) {
        rte->quantity[i] = va_arg (ap, int);
        rte->dtype[i]    = va_arg (ap, MPI_Datatype);
        rte->btype[i]    = _torc_mpi2b_type(rte->dtype[i]);
        rte->callway[i]  = va_arg (ap, int);

        if ((rte->callway[i] == CALL_BY_COP)&&(rte->quantity[i] > 1)) rte->callway[i] = CALL_BY_COP2;

#if DBG
        printf("ARG %d : Q = %d, T = %d, C = %x O\n", i, rte->quantity[i], rte->dtype[i], rte->callway[i]); fflush(0);
#endif
    }

    for (i=0; i<narg; i++) {
        if (rte->quantity[i] == 0) {    // peh: 02.07.2015
            VIRT_ADDR dummy = va_arg (ap, VIRT_ADDR);
            continue;
        }
        if (rte->callway[i] == CALL_BY_COP) {
            int typesize;
            MPI_Type_size(rte->dtype[i], &typesize);
            switch (typesize) {
                case 4:
                    rte->localarg[i] = *va_arg (ap, INT32 *);
                    break;
                case 8:
                    rte->localarg[i] = *va_arg (ap, INT64 *);
                    break;
                default:
                    Error("typesize not 4 or 8!");
                    break;
            }
        }
        else if (rte->callway[i] == CALL_BY_COP2) {
            VIRT_ADDR addr = va_arg (ap, VIRT_ADDR);
            int typesize;
            void *pmem;
            MPI_Type_size(rte->dtype[i], &typesize);
            pmem = malloc(rte->quantity[i]*typesize);
            memcpy(pmem, (void *)addr, rte->quantity[i]*typesize);
            rte->localarg[i] = (INT64)pmem; //yyyyyy
        }
        else {
            rte->localarg[i] = va_arg (ap, VIRT_ADDR);    /* pointer */
        }
    }

    if (queue == -1) {
        torc_to_rq_end(rte);
    }
    else {
        torc_to_lrq_end(queue, rte);
    }
}

double torc_gettime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

int torc_getlevel()
{
    torc_t *self = _torc_self();
    return self->level;
}

int torc_node_id()
{
    return mpi_rank;
}

int torc_num_nodes()
{
    return mpi_nodes;
}

int torc_i_num_workers()
{
    return kthreads;
}

int torc_i_worker_id(void)
{
    return _torc_get_vpid();
}

int torc_num_workers()
{
    if (torc_num_nodes() > 1)
        return _torc_total_num_threads();
    else
        return torc_i_num_workers();
}

int torc_worker_id(void)
{
    if (torc_num_nodes() > 1)
        return local_thread_id_to_global_thread_id(_torc_get_vpid());
    else
        return _torc_get_vpid();
}

struct torc_data *torc_data;

void torc_init(int argc, char *argv[], int ms)
{
    static int torc_initialized = 0;

    if (torc_initialized) return;

    torc_initialized = 1;

    torc_data = calloc(1, sizeof(struct torc_data));
    _torc_opt(argc, argv);
    _torc_env_init();
    _torc_worker(0);
}

#if 1
void *torc_getarg_addr(int arg)
{
    torc_t *self = _torc_self();

    if (torc_node_id() == self->homenode) {
        if (self->callway[arg] == CALL_BY_COP)
            return &(self->localarg[arg]);
        else
            return ((void *)self->localarg[arg]);
    }
    else {
        if (self->callway[arg] == CALL_BY_COP)
            return &(self->temparg[arg]);
        else
            return ((void *)self->temparg[arg]);
    }
}

int torc_get_num_args(int arg)
{
    return _torc_self()->narg;
}

int torc_getarg_callway(int arg)
{
    return _torc_self()->callway[arg];
}

int torc_getarg_count(int arg)
{
    return _torc_self()->quantity[arg];
}

int torc_getarg_size(int arg)
{
    int typesize;

    MPI_Type_size(_torc_self()->dtype[arg], &typesize);
    return typesize;
}
#endif


#if 1

void F77_FUNC_(torc_taskinit, TORC_TASKINIT)()
{
    /* nothing to do */
}

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
void F77_FUNC_(torc_waitall, TORC_WAITALL)()
{
    torc_waitall();
}
#endif

void F77_FUNC_(torc_createf, TORC_CREATEF) (int *pqueue, void (* work) (), int *pnarg, ...)
{
    int queue = *pqueue;
    int narg = *pnarg;

    va_list ap;
    int i;
    torc_t * rte;
    torc_t *self = _torc_self();

    /* Check if rte_init has been called */
    _lock_acquire(&self->lock);
    if (self->ndep == 0) self->ndep = 1;
    _lock_release(&self->lock);

    _torc_depadd(self, 1);

    rte = _torc_get_reused_desc();
    _initialize(rte);

    _torc_set_work_routine(rte, work);

    rte->narg       = narg;
    rte->rte_type   = 1;    /* external */
    rte->inter_node = 1;
    rte->parent     = self;
    rte->level      = self->level + 1;

#ifdef TORC_STATS
    if (!invisible_flag) created[self->vp_id]++;

    if (invisible_flag) rte->rte_type = 2;    /* invisible */
#endif

    if (narg>MAX_TORC_ARGS) Error("narg > MAX_TORC_ARGS");

    va_start(ap, pnarg);

    for (i=0; i<narg; i++) {
        rte->quantity[i] = *va_arg (ap, int *);
        MPI_Fint dt      = *va_arg (ap, MPI_Fint *);
        rte->dtype[i]    = MPI_Type_f2c(dt);
        rte->btype[i]    = _torc_mpi2b_type(rte->dtype[i]);
        rte->callway[i]  = *va_arg (ap, int *);

        if ((rte->callway[i] == CALL_BY_COP)&&(rte->quantity[i] > 1)) rte->callway[i] = CALL_BY_COP2;

#if DBG
        printf("ARG %d : Q = %d, T = %d, C = %x O\n", i, rte->quantity[i], rte->dtype[i], rte->callway[i]); fflush(0);
#endif
    }

    for (i=0; i<narg; i++) {
        if (rte->callway[i] == CALL_BY_COP) {
            int typesize;
            MPI_Type_size(rte->dtype[i], &typesize);
            switch (typesize) {
                case 4:
                    rte->localarg[i] = *va_arg (ap, INT32 *);
                    break;
                case 8:
                    rte->localarg[i] = *va_arg (ap, INT64 *);
                    break;
                default:
                    Error("typesize not 4 or 8!");
                    break;
            }
        }
        else if (rte->callway[i] == CALL_BY_COP2) {
            VIRT_ADDR addr = va_arg (ap, VIRT_ADDR);
            int typesize;
            void *pmem;
            MPI_Type_size(rte->dtype[i], &typesize);
            pmem = malloc(rte->quantity[i]*typesize);
            memcpy(pmem, (void *)addr, rte->quantity[i]*typesize);
            rte->localarg[i] = (INT64)pmem; //yyyyyy
        }
        else {
            rte->localarg[i] = va_arg (ap, VIRT_ADDR);    /* pointer (C: PTR, VAL) */
        }
    }

    if (queue == -1) {
        torc_to_rq_end(rte);
    }
    else {
        torc_to_lrq_end(queue, rte);
    }
}

// this is here to support the new pndl version
void F77_FUNC_(torc_taskf, TORC_TASKF) (void (* work) (), int *ptype, int *pnarg, ...)
{
    int queue = torc_worker_id(); //*pqueue;
    int narg = *pnarg;
    int type = *ptype; // ignored

    va_list ap;
    int i;
    torc_t * rte;
    torc_t *self = _torc_self();

    /* Check if rte_init has been called */
    _lock_acquire(&self->lock);
    if (self->ndep == 0) self->ndep = 1;
    _lock_release(&self->lock);

    _torc_depadd(self, 1);

    rte = _torc_get_reused_desc();
    _initialize(rte);

    _torc_set_work_routine(rte, work);

    rte->narg       = narg;
    rte->rte_type   = 1;    /* external */
    rte->inter_node = 1;
    rte->parent     = self;
    rte->level      = self->level + 1;

#ifdef TORC_STATS
    if ((!invisible_flag)&&(!type)) created[self->vp_id]++;

    if ((invisible_flag)||(type)) rte->rte_type = 2;    /* invisible */
#endif

    if (narg>MAX_TORC_ARGS) Error("narg > MAX_TORC_ARGS");

    va_start(ap, pnarg);

    for (i=0; i<narg; i++) {
        rte->quantity[i] = *va_arg (ap, int *);
        MPI_Fint dt      = *va_arg (ap, MPI_Fint *);
        rte->dtype[i]    = MPI_Type_f2c(dt);
        rte->btype[i]    = _torc_mpi2b_type(rte->dtype[i]);
        rte->callway[i]  = *va_arg (ap, int *);

        if ((rte->callway[i] == CALL_BY_COP)&&(rte->quantity[i] > 1)) rte->callway[i] = CALL_BY_COP2;

#if DBG
        printf("ARG %d : Q = %d, T = %d, C = %x O\n", i, rte->quantity[i], rte->dtype[i], rte->callway[i]); fflush(0);
#endif
    }

    for (i=0; i<narg; i++) {
        if (rte->callway[i] == CALL_BY_COP) {
            int typesize;
            MPI_Type_size(rte->dtype[i], &typesize);
            switch (typesize) {
                case 4:
                    rte->localarg[i] = *va_arg (ap, INT32 *);
                    break;
                case 8:
                    rte->localarg[i] = *va_arg (ap, INT64 *);
                    break;
                default:
                    Error("typesize not 4 or 8!");
                    break;
            }
        }
        else if (rte->callway[i] == CALL_BY_COP2) {
            VIRT_ADDR addr = va_arg (ap, VIRT_ADDR);
            int typesize;
            void *pmem;
            MPI_Type_size(rte->dtype[i], &typesize);
            pmem = malloc(rte->quantity[i]*typesize);
            memcpy(pmem, (void *)addr, rte->quantity[i]*typesize);
            rte->localarg[i] = (INT64)pmem; //yyyyyy
        }
        else {
            rte->localarg[i] = va_arg (ap, VIRT_ADDR);    /* pointer (C: PTR, VAL) */
        }
    }

    if (queue == -1) {
        torc_to_rq_end(rte);
    }
    else {
        torc_to_lrq_end(queue, rte);
    }
}

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
int F77_FUNC_(torc_num_workers, TORC_NUM_WORKERS) (void)
{
    return torc_num_workers();
}
#endif

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
int F77_FUNC_(torc_worker_id, TORC_WORKER_ID) (void)
{
    return torc_worker_id();
}
#endif

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
int F77_FUNC_(torc_node_id, TORC_NODE_ID) (void)
{
    return torc_node_id();
}
#endif

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
int F77_FUNC_(torc_num_nodes, TORC_NUM_NODES) (void)
{
    return torc_num_nodes();
}
#endif

void F77_FUNC_(torc_broadcastf, TORC_BROADCASTF) (void *a, long *count, MPI_Fint *datatype)
{

    MPI_Datatype dt;

    dt = MPI_Type_f2c(*datatype);
    torc_broadcast(a, *count, dt);
}

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
void F77_FUNC_ (torc_enable_stealing, TORC_ENABLE_STEALING) ()
{
    torc_enable_stealing ();
}
#endif

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
void F77_FUNC_ (torc_disable_stealing, TORC_DISABLE_STEALING) ()
{
    torc_disable_stealing ();
}
#endif

int torc_sched_nextcpu(int cpu, int stride)
{
    int res;
    int ncpus = torc_num_workers();

    if (cpu == -1)
        cpu = torc_worker_id();
    else
        cpu = (cpu + stride) % ncpus;

    res = cpu;
    return res;
}

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
int F77_FUNC_(torc_sched_nextcpu, TORC_SCHED_NEXTCPU) (int *cpu, int *stride)
{
    return torc_sched_nextcpu(*cpu, *stride);
}
#endif

void F77_FUNC_(torc_initf, TORC_INITF) (int *mode)
{
    torc_init(0, NULL, *mode);
}

#if F77_FUNC_(f77fun, F77FUN) == f77fun
#else
void F77_FUNC_(torc_finalize, TORC_FINALIZE) ()
{
    torc_finalize();
}
#endif

void F77_FUNC_(fff, FFF) ()
{
    fflush(0);
}

#endif

#undef f77fun
