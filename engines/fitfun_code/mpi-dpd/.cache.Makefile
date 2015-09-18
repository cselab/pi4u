#h5part = 0
#NVCC = nvcc -I$(CRAY_MPICH2_DIR)/include -L$(CRAY_MPICH2_DIR)/lib -I$(HDF5_DIR)/include -I../
#CXX = CC $(CRAY_CUDATOOLKIT_POST_LINK_OPTS) $(CRAY_CUDATOOLKIT_INCLUDE_OPTS) -L$(HDF5_DIR)/lib -I../

h5 = 0
h5part = 0
MY_MPICH2_DIR=/scratch/daint/chatzidp/usr/mpich313
NVCC = nvcc -I$(MY_MPICH2_DIR)/include -L$(MY_MPICH2_DIR)/lib -I$(HDF5_DIR)/include -DNO_H5PART
CXX = mpic++ $(CRAY_CUDATOOLKIT_POST_LINK_OPTS) $(CRAY_CUDATOOLKIT_INCLUDE_OPTS) -L$(HDF5_DIR)/lib
