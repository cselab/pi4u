C
C     Predefined datatypes
C

      INCLUDE 'mpif.h'

C
C
C

      INTEGER MODE_MS
      INTEGER MODE_SPMD

      PARAMETER (MODE_MS   = 0)
      PARAMETER (MODE_SPMD = 1)

      INTEGER CALL_BY_COP
      INTEGER CALL_BY_REF
      INTEGER CALL_BY_RES
      INTEGER CALL_BY_PTR
      INTEGER CALL_BY_VAL
      INTEGER CALL_BY_COP2
      INTEGER CALL_BY_VAD

      PARAMETER (CALL_BY_COP  = x'0001')
      PARAMETER (CALL_BY_REF  = x'0002')
      PARAMETER (CALL_BY_RES  = x'0003')
      PARAMETER (CALL_BY_PTR  = x'0004')
      PARAMETER (CALL_BY_VAL  = x'0001')
      PARAMETER (CALL_BY_COP2 = x'0005')
      PARAMETER (CALL_BY_VAD  = x'0006')

C
C     INTERFACE
C

      INTEGER*4 torc_worker_id
      INTEGER*4 torc_num_workers
      INTEGER*4 torc_i_worker_id
      INTEGER*4 torc_i_num_workers
      INTEGER*4 torc_node_id
      INTEGER*4 torc_num_nodes
      INTEGER*4 torc_sched_nextcpu

      DOUBLE PRECISION torc_gettime

      EXTERNAL torc_finalize
      EXTERNAL torc_taskinit
      EXTERNAL torc_waitall
      EXTERNAL torc_register_task
      EXTERNAL torc_enable_stealing
      EXTERNAL torc_disable_stealing
      EXTERNAL torc_sleep

      EXTERNAL torc_broadcastf
      EXTERNAL torc_initf
      EXTERNAL torc_createf
      EXTERNAL torc_taskf


