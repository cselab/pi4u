/*
 *  abc_subsim.cpp
 *  Pi4U
 *
 *  Created by Lina Kulakova on 1/1/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include <gsl/gsl_sort.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include <cstring>
#include <mutex>
#include <cmath>
//#define _PEH_SORTING_
//#define _FIXED_COV_

#if defined(_USE_TORC_)
#include <mpi.h>
extern "C"
{
#include <torc.h>
}
#else
#include <sys/time.h>
static double torc_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

int torc_node_id()
{
	return 0;
}

int torc_num_workers()
{
	return 1;
}
#endif

#include "abc_subsim.h"
#include "problem.h"


/*
 ********************************************************************************************************************************************
 * Classes and Data for database and covariance matrix
 ********************************************************************************************************************************************
 */

/*
 * peh: (II) class subsim_database
 */

// 3-way database of max size = population
class subsim_database
{
private:
	int db_ind;	// peh: db_ind(ex) to add entries
	std::mutex mutex;

public:
	static subsim_database *me;

	std::vector< gsl_vector * > sample_array;
	gsl_vector * discrepancy_array;

	subsim_database()
	{
		subsim_database::me = this;
	}

	void reset()
	{
		db_ind = 0;
	}

	void update(double *sample_v, double discrepancy_sample_v)
	{
		gsl_vector_view sample;
		sample = gsl_vector_view_array(sample_v, PROBDIM);

		int ind;

		mutex.lock();
		{
		ind = db_ind;
		db_ind += 1;
		}
		mutex.unlock();

		gsl_vector_memcpy(sample_array[ind], &sample.vector);
		gsl_vector_set(discrepancy_array, ind, discrepancy_sample_v);
	}

	static void update_task(double *sample_v, double *p_discrepancy_sample_v)
	{
		subsim_database::me->update(sample_v, *p_discrepancy_sample_v);
	}

	void torc_update(double *sample_v, double discrepancy_sample_v)
	{
		if (torc_node_id() == 0) {
			subsim_database::me->update(sample_v, discrepancy_sample_v);
		}
		else {
#if defined(_USE_TORC_)
			torc_create_direct(0, (void (*)())subsim_database::update_task, 2,
						PROBDIM, MPI_DOUBLE, CALL_BY_COP,
						1, MPI_DOUBLE, CALL_BY_COP,
						sample_v, &discrepancy_sample_v);
			torc_waitall3();
#endif
		}
	}
};

subsim_database * subsim_database::me = NULL;

/*
 * peh: (III) global data
 */

subsim_database db;



/*
 * peh: (IV) global data for covariance matrix
 */

gsl_matrix * /*abc_subset::*/ covariance;

#if 0
static void cov_update_task()
{
#if defined(_USE_TORC_)
	MPI_Bcast(covariance->data, PROBDIM*PROBDIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	printf("after MPI_Bcast on node %d\n", torc_node_id()); fflush(0);
#endif
}

void spmd_cov_update()
{
#if defined(_USE_TORC_)
	int i;
	if (torc_num_nodes() == 1) return;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())cov_update_task, 0);
	}
	torc_waitall();
#endif
}
#else
static void cov_update_task(double *cov_data)
{
#if defined(_USE_TORC_)
	memcpy(covariance->data, cov_data, PROBDIM*PROBDIM*sizeof(double));
	if (torc_node_id() == 1) {
		printf("after MPI_Bcast on node %d [%lf %lf %lf %lf]\n", torc_node_id(), 
			covariance->data[0], covariance->data[1], covariance->data[2], covariance->data[3]); fflush(0);
	}
#endif
}

void spmd_cov_update()
{
#if defined(_USE_TORC_)
	int i;
	if (torc_num_nodes() == 1) return;
	torc_disable_stealing();
	for (i = 0; i < torc_num_nodes(); i++) {
		if (i == torc_node_id()) continue;
		torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())cov_update_task, 1,
				PROBDIM*PROBDIM, MPI_DOUBLE, CALL_BY_VAL,
				covariance->data);
	}
	torc_waitall();
	printf("after MPI_Bcast on node %d [%lf %lf %lf %lf]\n", torc_node_id(), 
		covariance->data[0], covariance->data[1], covariance->data[2], covariance->data[3]); fflush(0);
#endif
}
#endif


/*
 ********************************************************************************************************************************************
 * Some helping functions
 ********************************************************************************************************************************************
 */

void gsl_vector_fprintf_as_line(FILE * fp, const gsl_vector * v, const char * format)
{
	for(int i=0; i<v->size; ++i)
	{
		fprintf(fp, format, gsl_vector_get(v, i));
		fprintf(fp, " ");
	}
}

void gsl_matrix_fprintf_as_matrix(FILE * fp, const gsl_matrix * mat, const char * format)
{
	for(int i=0; i<mat->size1; ++i)
	{
		for(int j=0; j<mat->size2; ++j)
		{
			fprintf(fp, format, gsl_matrix_get(mat, i, j));
			fprintf(fp, " ");
		}
		fprintf(fp, "\n");
	}
}

void get_sum(std::vector< gsl_vector * > const& array, gsl_vector * sum)
{
	gsl_vector_set_zero(sum);

	for(int i=0; i<array.size(); ++i)
		gsl_blas_daxpy(1.0, array[i], sum);
}

void get_diag(const gsl_matrix * mat, gsl_vector * diag)
{
	for(int i=0; i<mat->size1; ++i)
		gsl_vector_set(diag, i, gsl_matrix_get(mat, i, i));
}

void make_diag(const gsl_vector * diag, gsl_matrix * mat)
{
	gsl_matrix_set_zero(mat);
	for(int i=0; i<diag->size; ++i)
		gsl_matrix_set(mat, i, i, gsl_vector_get(diag, i));
}

void compute_moments(std::vector< gsl_vector * > const& sample_array, gsl_vector * mean, gsl_matrix * covariance)
{
	// mean
	get_sum(sample_array, mean);
	int n = sample_array.size();
	gsl_vector_scale(mean, 1.0/n);

#if !defined(_FIXED_COV_)
	// covariance
	gsl_vector * diff = gsl_vector_alloc(mean->size);

	gsl_matrix_set_zero(covariance);
//	gsl_matrix_set_all(covariance, 0.1);
	for(int i=0; i<sample_array.size(); ++i)
	{
		gsl_vector_memcpy(diff, sample_array[i]);
		gsl_blas_daxpy(-1.0, mean, diff);
		gsl_matrix_view v = gsl_matrix_view_vector(diff, diff->size, 1);
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &v.matrix, &v.matrix, 1.0, covariance);
	}
	gsl_matrix_scale(covariance, 1.0/(n-1));

	gsl_vector_free(diff);
#else
	gsl_matrix_set_all(covariance, 0.0);
	for(int i=0; i<PROBDIM; ++i)
		gsl_matrix_set(covariance, i, i, 0.1*(upper_bounds[i]-lower_bounds[i]));
#endif
}


int gsl_vector_is_positive(const gsl_vector * v)
{
	for(int i=0; i<v->size; ++i)
		if(gsl_vector_get(v,i)*enforce_positivity[i] < 0)
			return 0;

	return 1;
}


/*
 ****************************************************************************************************************************************
 * Functions of abc_subsim
 ****************************************************************************************************************************************
 */

void abc_subsim::run()
{
	// workspace
	std::vector< gsl_vector * > sample_array_prev;  // peh: local
	for(int i=0; i<population_size; ++i)
	{
		db.sample_array.push_back(gsl_vector_alloc(PROBDIM));
		sample_array_prev.push_back(gsl_vector_alloc(PROBDIM));
	}
	gsl_vector * discrepancy_array_prev = gsl_vector_alloc(population_size);
	db.discrepancy_array = gsl_vector_alloc(population_size);

	gsl_vector * acc_rate_array = gsl_vector_alloc(pivot);
	gsl_vector * mean = gsl_vector_alloc(PROBDIM);

	size_t perm[population_size]; // maybe use gsl_permutation ?

	double qoi_prev;

	// step for level 0
	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! level 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	// restart
	if(restart)
		read_restart_file(db.sample_array, db.discrepancy_array);
	else
		sample_from_prior_in_parallel(db.sample_array, db.discrepancy_array);	// level = 0


#if 0
	return;
#endif

	// compute the qoi
	qoi_prev = 0.0;
	for(int ind=0; ind<population_size; ++ind)
		qoi_prev += gsl_blas_dnrm2(db.sample_array[ind]);
	qoi_prev /= population_size;
	std::cout << "QOI: " << qoi_prev << std::endl; 

	double tolerance_prev = gsl_vector_max(db.discrepancy_array);
	double mean_acc_rate = 0.99;

	// steps for other levels
	int level;
	for(level=1; level<max_level; ++level)
	{
		// get achieved tolerance level
		gsl_sort_index(perm, db.discrepancy_array->data, db.discrepancy_array->stride, db.discrepancy_array->size);
		double tolerance_curr = 0.5*(gsl_vector_get(db.discrepancy_array,perm[pivot]) + gsl_vector_get(db.discrepancy_array,perm[pivot+1]));

		// store previous values
		for(int i=0; i<population_size; ++i)
		{
			gsl_vector_memcpy(sample_array_prev[i], db.sample_array[i]);
		}
		gsl_vector_memcpy(discrepancy_array_prev, db.discrepancy_array);
		tolerance_prev = tolerance_curr;

		write_array_to_file(db.sample_array, db.discrepancy_array, level-1);

		// compute moments
		compute_moments(sample_array_prev, mean, covariance);

		std::cout << "mean: "; gsl_vector_fprintf_as_line(stdout, mean, "%le"); std::cout << std::endl;
		std::cout << "covariance:\n"; gsl_matrix_fprintf_as_matrix(stdout, covariance, "%le");
		std::cout << "weights std: " << gsl_stats_sd(db.discrepancy_array->data, db.discrepancy_array->stride, db.discrepancy_array->size) << std::endl;

		std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! level " << level << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "tolerance " << tolerance_curr << std::endl;


		// propose a covariance - xxx
#if !defined(_FIXED_COV_)
//		gsl_matrix_scale(covariance, 0.1);

//		tune_covariance(sample_array_prev, discrepancy_array_prev, tolerance_curr, db.sample_array,
//				db.discrepancy_array, acc_rate_array, covariance, perm, level);
		tune_covariance_in_parallel(sample_array_prev, discrepancy_array_prev, tolerance_curr, db.sample_array,
				db.discrepancy_array, perm, level);
#endif

		std::cout << "proposal covariance:\n"; gsl_matrix_fprintf_as_matrix(stdout, covariance, "%le");

		spmd_cov_update();	// peh: update covariance on all nodes (switch to SPMD and broadcast)
		db.reset();		// peh: reset the database

#if 0
	return;
#endif

		// run MH
		run_mh_in_parallel(chain_length, perm, sample_array_prev, discrepancy_array_prev, tolerance_curr, acc_rate_array, level);

		double mean_acc_rate = gsl_stats_mean(acc_rate_array->data, acc_rate_array->stride, acc_rate_array->size);
		std::cout << "mean acc_rate: " << mean_acc_rate << std::endl;

		// check stopping criteria
		// TODO: add stagnation criteria (if the tolerance is around the same number for the last several iterations
		// or goes up and down)
		if(tolerance_final/tolerance_curr > 1)
		{
			std::cout << "Final tolerance " << tolerance_final << " achieved, stop." << std::endl;
			break;
		}

		if(gsl_stats_sd(db.discrepancy_array->data, db.discrepancy_array->stride, db.discrepancy_array->size) < weights_std)
		{
			std::cout << "Weights std is below the specified threshold " << weights_std << ", stop." << std::endl;
			break;
		}

		if(mean_acc_rate < acc_rate)
		{
			std::cout << "Acceptance rate is below the specified threshold " << acc_rate << ", stop." << std::endl;
			// return previous values
			for(int i=0; i<population_size; ++i)
				gsl_vector_memcpy(db.sample_array[i], sample_array_prev[i]);
			gsl_vector_memcpy(db.discrepancy_array, discrepancy_array_prev);
			break;
		}

		// compute the qoi
		double qoi_curr = 0.0;
		for(int ind=0; ind<population_size; ++ind)
			qoi_curr += gsl_blas_dnrm2(db.sample_array[ind]);
		qoi_curr /= population_size;
		std::cout << "QOI: " << qoi_curr << std::endl; 

		if(fabs(qoi_prev-qoi_curr)/fabs(qoi_prev) < 0)
		{
			std::cout << "QOI relative change is below 0.01, stop." << std::endl;
			break;
		}

		qoi_prev = qoi_curr;
	}

	// print the solution
	compute_moments(db.sample_array, mean, covariance);

	std::cout << "mean: "; gsl_vector_fprintf_as_line(stdout, mean, "%le"); std::cout << std::endl;
	std::cout << "covariance:\n"; gsl_matrix_fprintf_as_matrix(stdout, covariance, "%le");
	std::cout << "number of evaluations: " << eval_count << std::endl;

	write_array_to_file(db.sample_array, db.discrepancy_array, level-1);

	// free memory
	for(int i=0; i<population_size; ++i)
	{
		gsl_vector_free(db.sample_array[i]);
		gsl_vector_free(sample_array_prev[i]);
	}

	gsl_vector_free(db.discrepancy_array);
	gsl_vector_free(discrepancy_array_prev);

	gsl_vector_free(acc_rate_array);
	gsl_vector_free(mean);
	gsl_matrix_free(covariance);
}


void abc_subsim::read_restart_file(std::vector< gsl_vector * > & sample_array, gsl_vector * discrepancy_array)
{
	std::cout << "Reading restart file " << restart_filename << std::endl;

	std::ifstream ifs;
	ifs.open(restart_filename);
	if(!ifs.good())
	{
		std::cout << "Cannot read file " << restart_filename << ", stop." << std::endl;
		exit(1);
	}

	double t0 = torc_gettime();
	for(int ind=0; ind<population_size; ++ind)
	{
		// read sample
		for(int k=0; k<PROBDIM; ++k)
		{
			double s;
			ifs >> s;
			gsl_vector_set(sample_array[ind], k, s*normalization[k]);
		}

		// read discrepancy
		double d;
		ifs >> d;
		gsl_vector_set(discrepancy_array, ind, d);
	}
	double t1 = torc_gettime();

	ifs.close();
	printf("t1-t0 = %lf seconds\n", t1-t0);
}

#if defined(_PEH_SORTING_)
// peh: add index for sorting as user option
bool mycompare(gsl_vector *a, gsl_vector *b)
{
	return (gsl_vector_get(a, 0) < gsl_vector_get(b, 0));
}

bool mycompare_desc(gsl_vector *a, gsl_vector *b)
{
	return (gsl_vector_get(a, 0) > gsl_vector_get(b, 0));
}

void sort_samples(std::vector< gsl_vector * > & sample_array)
{
	std::sort(sample_array.begin(), sample_array.end(), mycompare);
}

void sort_samples_v2(std::vector< gsl_vector * > & sample_array, int stride)
{
	int n = sample_array.size();

	printf("a: sorting from %d to %d\n", 0, n);
//	std::sort(sample_array.begin(), sample_array.end(), mycompare); // ascending order
	std::sort(sample_array.begin(), sample_array.begin()+n, mycompare);     // ascending order

	int chunks = n / stride;

	int i = stride;
	while (i < n) {
		printf("b: sorting from %d to %d\n", i, i+stride);
		if (i + stride < n)
			std::sort(sample_array.begin()+i, sample_array.begin()+i+stride, mycompare_desc);
		else
			std::sort(sample_array.begin()+i, sample_array.end(), mycompare_desc);

		i += 2*stride;
	}
}

#endif

void abc_subsim::sample_from_prior_in_parallel(std::vector< gsl_vector * > & sample_array, gsl_vector * discrepancy_array)
{	
	double t0 = torc_gettime();
	int info[4];	/* {level, chain, step, task} */ 

#if !defined(_PEH_SORTING_)

	for(int ind=0; ind<population_size; ++ind)
	{
		// get random sample according to prior
		generate_sample_by_prior(sample_array[ind]);
//		std::cout << "Sample # " << ind << ": "; gsl_vector_fprintf_as_line(stdout, sample_array[ind], "%lf"); std::cout << std::endl;
	
		info[0] = 0; info[1] = 0; info[2] = 0; info[3] = ind;
#if defined(_USE_TORC_)
		torc_create(-1, (void (*)())run_problem, 3,
				PROBDIM, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				sample_array[ind]->data, &discrepancy_array->data[ind], info);  // peh: accessible as db.xxx
#else
		run_problem(sample_array[ind]->data, &discrepancy_array->data[ind], info);
#endif
	}

#else	// peh:sort

	// generate random samples, sort them according to their sigma value, submit tasks
	for(int ind=0; ind<population_size; ++ind)
	{
		// get random sample according to prior
		generate_sample_by_prior(sample_array[ind]);
//		std::cout << "Sample # " << ind << ": "; gsl_vector_fprintf_as_line(stdout, sample_array[ind], "%lf"); std::cout << std::endl;
	}

	//sort_samples(sample_array);
	sort_samples_v2(sample_array, torc_num_workers());	// this helps for the sampling phase

	for(int ind=0; ind<population_size; ++ind)
	{
		info[0] = 0; info[1] = 0; info[2] = 0; info[3] = ind;

#if defined(_USE_TORC_)
		torc_create((ind+0)%torc_num_workers(), (void (*)())run_problem, 3,
//		torc_create(-1, (void (*)())run_problem, 3,
				PROBDIM, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				sample_array[ind]->data, &discrepancy_array->data[ind], info);  // peh: accessible as db.xxx
#else
		run_problem(sample_array[ind]->data, &discrepancy_array->data[ind], info);
#endif
	}

#endif


#if defined(_USE_TORC_)
	torc_enable_stealing();
	torc_waitall();
	torc_disable_stealing();
#endif

	double t1 = torc_gettime();
	printf("sample_from_prior_in_parallel: elapsed timet = %lf seconds\n", t1-t0);
}

void abc_subsim::tune_covariance(std::vector< gsl_vector * > const& sample_array_prev,
		const gsl_vector * discrepancy_array_prev, double tolerance_curr, std::vector< gsl_vector * > & sample_array,
		gsl_vector * discrepancy_array_curr, gsl_vector * acc_rate_array,
		gsl_matrix * covariance, const size_t * perm, int level)
{
	double c_opt = (PROBDIM == 1 ? 0.45 : 0.234);
	double scale = (2.38/sqrt(PROBDIM));

	gsl_matrix * test_covariance = gsl_matrix_alloc(PROBDIM, PROBDIM);

	for(int count=0; count<10; ++count)
	{
		printf("trying scale^2 = %f\n", scale*scale);

		gsl_matrix_memcpy(test_covariance, covariance);
		gsl_matrix_scale(test_covariance, scale*scale);

//		ind = perm[random.randint(0,self.pivot-1)]
		int ind = perm[0];


		// peh: this version of run_mh_test is almost identical to the original one; I kept it as it is called sequentially and locally
		abc_subsim::run_mh_test(
					chain_length_test,
					sample_array_prev[ind]->data,
					gsl_vector_get(discrepancy_array_prev, ind),
					test_covariance,
					tolerance_curr,
					&acc_rate_array->data[0],
					0, level);

		double acc_rate = gsl_vector_get(acc_rate_array, 0);
		std::cout << "acc_rate = " << acc_rate << std::endl;

		if(acc_rate <= 0.01)
			acc_rate = 0.01;
		if(acc_rate >= 0.99)
			acc_rate = 0.99;
//		if(acc_rate < 0.15 || acc_rate > 0.4)
		if(acc_rate < 0.3 || acc_rate > 0.8)
			scale *= gsl_cdf_ugaussian_Pinv(0.5*c_opt) / gsl_cdf_ugaussian_Pinv(0.5*acc_rate);
		else
			break;
	}

	printf("final scale^2 = %f\n", scale*scale);
	gsl_matrix_memcpy(covariance, test_covariance);

	// free workspace
	gsl_matrix_free(test_covariance);
}

// Try one scaling parameter for the covariance
void abc_subsim::tune_covariance_task(
		int *pchain_length_test,	// in,1
		const double *sample_v,		// in,PROBDIM
		double *pdiscrepancy_v,		// in,1
		double *ptolerance_curr,	// in,1
		double *acc_rate_array_val,	// out,1
		double *pscale,			// in,1
		int *pnumber,			// in,1
		int *plevel)			// in,1
{
	int chain_length_test = *pchain_length_test;
	double discrepancy_v = *pdiscrepancy_v;
	double tolerance_curr = *ptolerance_curr;
	double scale = *pscale;
	int number = *pnumber;
	int level = *plevel;

	gsl_matrix * test_covariance = gsl_matrix_alloc(PROBDIM, PROBDIM);

	gsl_matrix_memcpy(test_covariance, covariance);
	gsl_matrix_scale(test_covariance, scale*scale);

	abc_subsim::run_mh_test(
				chain_length_test,
				sample_v,
				discrepancy_v,
				test_covariance,
				tolerance_curr,
				acc_rate_array_val,	// out
				number, level);

//	double acc_rate = gsl_vector_get(acc_rate_array, number);
//	std::cout << "acc_rate = " << acc_rate << std::endl;

	// free workspace
	gsl_matrix_free(test_covariance);
}


void abc_subsim::tune_covariance_in_parallel(
		std::vector< gsl_vector * > const& sample_array_prev,
		const gsl_vector * discrepancy_array_prev,
		double tolerance_curr,
		std::vector< gsl_vector * > & sample_array,
		gsl_vector * discrepancy_array_curr,
		const size_t * perm,
		int level)
{
	int attempts = 10;	// user input
	int npoints = 5;	// user input
	int ntries = 5;		// user input
	//double min_scale = 1e-1;
	//double max_scale = 1e+0;

	spmd_cov_update();	// peh: update covariance on all nodes (necessary for parallel tuning)

	std::vector< double > scale(attempts,0);
	double scale_value = 2.0;
	for(int count=attempts-1; count>=0; --count) {
		scale[count] = sqrt(scale_value);
		scale_value *= 0.5;	// user input
	}
	for(int count=0; count<attempts; ++count) {
		printf("scale[%d]^2 = %f\n", count, scale[count]*scale[count]); fflush(0);
	}

	gsl_vector * acc_rate_array = gsl_vector_alloc(attempts);

	double acc_rate_mat[attempts][npoints][ntries];

	chain_length_test = 1;

	int task_id = 0;
	for(int iattempt=0; iattempt<attempts; ++iattempt)
	{
	for(int ipoint=0; ipoint<npoints; ++ipoint)
	{
	for(int itry=0; itry<ntries; ++itry)
	{
		//printf("scale[%d]^2 = %f\n", iattempt, scale[iattempt]*scale[iattempt]);

		int ind = perm[ipoint];
		double discrepancy_array_prev_ind_val = gsl_vector_get(discrepancy_array_prev, ind);

		// ADD A TORC CALL HERE:
#if defined(_USE_TORC_)
		torc_create(-1, (void (*)())tune_covariance_task, 8,
				1, MPI_INT, CALL_BY_COP,
				PROBDIM, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				1, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				&chain_length_test,
				(const double *) (sample_array_prev[ind]->data),
				&discrepancy_array_prev_ind_val,
				&tolerance_curr,
				&acc_rate_mat[iattempt][ipoint][itry],
				//&acc_rate_array->data[iattempt],
				&scale[iattempt],
				//&iattempt,
				&task_id,
				&level);
#else
		tune_covariance_task(	&chain_length_test,
					(const double *) (sample_array_prev[ind]->data),
					&discrepancy_array_prev_ind_val,
					&tolerance_curr,
					//&acc_rate_array->data[iattempt],
					&acc_rate_mat[iattempt][ipoint][itry],
					&scale[iattempt],
					//&iattempt,
					&task_id,
					&level);
#endif
		task_id++;
	}
	}
	}

#if defined(_USE_TORC_)
	torc_waitall();
#endif

	for(int iattempt=0; iattempt<attempts; ++iattempt)
	{
		acc_rate_array->data[iattempt] = 0;
		for(int ipoint=0; ipoint<npoints; ++ipoint)
		for(int itry=0; itry<ntries; ++itry)
			acc_rate_array->data[iattempt] += acc_rate_mat[iattempt][ipoint][itry];

		acc_rate_array->data[iattempt] /= npoints*ntries;
	}

	// Find a good acc rate; if not found use some value (not min, not max)
	double min_acc_rate = 10;
	double max_acc_rate = -10;
	double my_acc_rate = -100;
	int min_index = -1;
	int max_index = -1;
	for(int count=0; count<attempts; ++count)
	{
		double acc_rate = gsl_vector_get(acc_rate_array, count);
		printf("acc_rate[%d] = %.3f\n", count, acc_rate);

		if(acc_rate < min_acc_rate)
		{
			min_index = count;
			min_acc_rate = acc_rate;
		}
		if(acc_rate > max_acc_rate)
		{
			max_index = count;
			max_acc_rate = acc_rate;
		}
	}
	printf("min_index = %d (%f) max_index = %d (%f)\n", min_index, min_acc_rate, max_index, max_acc_rate);

	double low_acc_rate = 0.2;
	double high_acc_rate = 0.5;
	// try to find an acceptance rate in the range [0.3,0.8]	// paper: [0.2,0.4]
	int my_index = -1;
	for(int count=0; count<attempts; ++count)
	{
		double acc_rate = gsl_vector_get(acc_rate_array, count);
		if ((low_acc_rate <= acc_rate) && (acc_rate <= high_acc_rate)) {
			my_index = count;
			break;
		}
	}

	if (my_index == -1) {	// nothing found in the desired range, just pick one expect from the min and max values
		my_index = 0;
		for(int count=0; count<attempts; ++count)
		{
			if(count != min_index && count != max_index)
			{
				my_index = count;
				break;
			}
		}
	}

	printf("final scale[%d]^2 = %f\n", my_index, scale[my_index]*scale[my_index]);

	gsl_matrix_scale(covariance, scale[my_index]*scale[my_index]);

	gsl_vector_free(acc_rate_array);

}

void abc_subsim::tune_covariance_in_parallel_v0(
		std::vector< gsl_vector * > const& sample_array_prev,
		const gsl_vector * discrepancy_array_prev,
		double tolerance_curr,
		std::vector< gsl_vector * > & sample_array,
		gsl_vector * discrepancy_array_curr,
		const size_t * perm,
		int level)
{
	double c_opt = (PROBDIM == 1 ? 0.45 : 0.234);
	int attempts = std::min(std::max(torc_num_workers(),10),10);
	double min_scale = 1e-1;
	double max_scale = 1e+0;

	std::vector< double > scale(attempts,0);
	gsl_vector * acc_rate_array = gsl_vector_alloc(attempts);

	spmd_cov_update();	// peh: update covariance on all nodes (necessary for parallel tuning)

	for(int count=0; count<attempts; ++count)
	{
		scale[count] = (2.38/sqrt(PROBDIM)) * (count*(max_scale-min_scale)/attempts + min_scale);

		printf("scale[%d]^2 = %f\n", count, scale[count]*scale[count]);

		int ind = perm[0];
		double discrepancy_array_prev_ind_val = gsl_vector_get(discrepancy_array_prev, ind);

		// ADD A TORC CALL HERE:
#if defined(_USE_TORC_)
		torc_create(-1, (void (*)())tune_covariance_task, 8,
				1, MPI_INT, CALL_BY_COP,
				PROBDIM, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				1, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				&chain_length_test,
				(const double *) (sample_array_prev[ind]->data),
				&discrepancy_array_prev_ind_val,
				&tolerance_curr,
				&acc_rate_array->data[count],
				&scale[count],
				&count,
				&level);
#else
		tune_covariance_task(	&chain_length_test,
					(const double *) (sample_array_prev[ind]->data),
					&discrepancy_array_prev_ind_val,
					&tolerance_curr,
					&acc_rate_array->data[count],
					&scale[count],
					&count,
					&level);
#endif
	}

#if defined(_USE_TORC_)
	torc_waitall();
#endif

	// Find a good acc rate; if not found use some value (not min, not max)
	double min_acc_rate = 10;
	double max_acc_rate = -10;
	double my_acc_rate = -100;
	int min_index = -1;
	int max_index = -1;
	for(int count=0; count<attempts; ++count)
	{
		double acc_rate = gsl_vector_get(acc_rate_array, count);
		printf("acc_rate[%d] = %.3f\n", count, acc_rate);

		if(acc_rate < min_acc_rate)
		{
			min_index = count;
			min_acc_rate = acc_rate;
		}
		/*else*/
		if(acc_rate > max_acc_rate)
		{
			max_index = count;
			max_acc_rate = acc_rate;
		}
	}
	printf("min_index = %d (%f) max_index = %d (%f)\n", min_index, min_acc_rate, max_index, max_acc_rate);

	double low_acc_rate = 0.3;
	double high_acc_rate = 0.8;
	// try to find an acceptance rate in the range [0.3,0.8]	// paper: [0.2,0.4]
	int my_index = -1;
	for(int count=0; count<attempts; ++count)
	{
		double acc_rate = gsl_vector_get(acc_rate_array, count);
		if ((low_acc_rate <= acc_rate) && (acc_rate <= high_acc_rate)) {
			my_index = count;
			break;
		}
	}

	if (my_index == -1) {	// nothing found in the desired range, just pick one expect from the min and max values
		my_index = 0;
		for(int count=0; count<attempts; ++count)
		{
			if(count != min_index && count != max_index)
			{
				my_index = count;
				break;
			}
		}
	}

	printf("final scale[%d]^2 = %f\n", my_index, scale[my_index]*scale[my_index]);

	gsl_matrix_scale(covariance, scale[my_index]*scale[my_index]);

	gsl_vector_free(acc_rate_array);

}
//#endif

void generate_sample_mh(const gsl_vector * mean, const gsl_matrix * covariance, gsl_vector * sample)
{
	do
	{
		ran_mvn(mean, covariance, sample);
	}
	while(!gsl_vector_is_positive(sample) || !evaluate_prior(sample));
//	while(!gsl_vector_is_positive(sample));
}

//int accept_mh(const gsl_vector * candidate, const gsl_vector * sample, const gsl_vector * data_sim_candidate,
//		const gsl_vector * data_sim_sample, double discrepancy_candidate,  double discrepancy_sample, double tolerance)

int accept_mh(const gsl_vector * candidate, const gsl_vector * sample, double discrepancy_candidate, double discrepancy_sample,
		double tolerance)
{
#ifndef PROBABILISTIC
	double r = evaluate_prior(candidate) / evaluate_prior(sample);
	return (gsl_ran_flat(rng.r[torc_i_worker_id()], 0.0, 1.0) < r)*(discrepancy_candidate/tolerance < 1);
//#else
//	// workspace
//	gsl_vector * distance = gsl_vector_alloc(NDSIM);
//	gsl_vector * err_mean = gsl_vector_calloc(NDSIM); // set to zero
//	gsl_vector * err_var = gsl_vector_alloc(NDSIM);
//	gsl_matrix * err_cov = gsl_matrix_alloc(NDSIM, NDSIM);
//	int isinf;
//
//	// candidate
//	get_error_var(candidate, err_var);
//	make_diag(err_var, err_cov);
//	get_distance(data_sim_candidate, distance);
//
//	double p1 = ran_mvn_logpdf(distance, err_mean, err_cov);
//	double r1 = evaluate_logprior(candidate, &isinf);
//	if(isinf)
//		return 0;
//
//	// sample
//	get_error_var(sample, err_var);
//	make_diag(err_var, err_cov);
//	get_distance(data_sim_sample, distance);
//
//	double p2 = ran_mvn_logpdf(distance, err_mean, err_cov);
//	double r2 = evaluate_logprior(sample, &isinf);
//
//	// free workspace
//	gsl_vector_free(distance);
//	gsl_vector_free(err_mean);
//	gsl_vector_free(err_var);
//	gsl_matrix_free(err_cov);
//
//	double r = gsl_ran_flat(rng.r[torc_i_worker_id()], 0, 1);
//	if(r <= 0 || (p1+r1)-(p2+r2) >= 0)
//		return 1;
//	else
//		return ((p1+r1)-(p2+r2) >= log(r));
#endif
}

void abc_subsim::run_mh_test(
			int chain_length,				// input
			//const gsl_vector * seed,			// input
			const double *seed_v,
			double discrepancy_seed,			// input
			const gsl_matrix * covariance,			// input
			double tolerance,				// input
			double *acc_rate_array_val,			// output to parent
			int seed_ind,					// input: task id
			int level)					// input: level
{

	// parameters
	int ind_begin = seed_ind*chain_length;
	int every = 1;

	// workspace
	gsl_vector * candidate = gsl_vector_alloc(PROBDIM);
	gsl_vector * sample = gsl_vector_alloc(PROBDIM);

	// initial values
	//gsl_vector_memcpy(sample, seed);
	memcpy(sample->data, seed_v, PROBDIM*sizeof(double));
	double discrepancy_sample = discrepancy_seed;

	// Generate chain samples
	int info[4]; /* level, chain, step, task */
	int accept_count = 0;
	int ind = 0;
	for(int count=0; count<every*chain_length; ++count)
	{
		// Step 1: generate a candidate sample
		generate_sample_mh(sample, covariance, candidate);

		info[0] = level; info[1] = seed_ind; info[2] = 0; info[3] = count;	// use info[3] to distinguish them from normal tasks
		double discrepancy_candidate;
		run_problem(candidate->data, &discrepancy_candidate, info);

//		std::cout << "Sample # " << count << ": "; gsl_vector_fprintf_as_line(stdout, candidate, "%lf"); std::cout << std::endl;

		// Step 2: accept/reject the candidate
		if(accept_mh(candidate, sample, discrepancy_candidate, discrepancy_sample, tolerance))
		{
			gsl_vector_memcpy(sample, candidate);
			discrepancy_sample = discrepancy_candidate;
			accept_count += 1;
		}

		//if(count%every == (every-1))
		//{
		//	gsl_vector_memcpy(sample_array[ind_begin+ind], sample);
		//	gsl_vector_set(discrepancy_array, ind_begin+ind, discrepancy_sample);
		//	ind += 1;
		//}
	}

	*acc_rate_array_val = (double)accept_count/(every*chain_length);

	// free workspace
	gsl_vector_free(candidate);
	gsl_vector_free(sample);
}

// the 'run_mh' task
void abc_subsim::run_mh(int *p_chain_length,				// input, 1
			const double * seed_v,				// input, PROBDIM
			double *p_discrepancy_seed,			// input, 1
			double *p_tolerance,				// input, 1
			double *acc_rate_array_val,			// output, 1
			int *p_seed_ind,				// input, 1
			int *p_level)					// input, 1
{
	// task parameters (IN, 1)
	int chain_length = *p_chain_length;
	double discrepancy_seed = *p_discrepancy_seed;
	double tolerance = *p_tolerance;
	int seed_ind = *p_seed_ind;
	int level = *p_level;

#if DBG
	printf("worker %d, seed_v:", torc_worker_id());
	for (int i = 0; i < PROBDIM; i++) printf("%f ", seed_v[i]);
	printf("\n");
#endif

	// parameters
	int ind_begin = seed_ind*chain_length;
	int every = 1;

	// workspace
	gsl_vector * candidate = gsl_vector_alloc(PROBDIM);
	gsl_vector * sample = gsl_vector_alloc(PROBDIM);

	memcpy(sample->data, seed_v, PROBDIM*sizeof(double));

	double discrepancy_sample = discrepancy_seed;

	// Keep the seed
	db.torc_update(sample->data, discrepancy_sample);

	// Generate chain samples
	int accept_count = 0;
	int ind = 0;
	int info[4];	/* level, chain, step, task */
	for(int count=0; count<every*(chain_length-1); ++count)
	{
		// Step 1: generate a candidate sample
		generate_sample_mh(sample, covariance, candidate);

		info[0] = level; info[1] = seed_ind; info[2] = count; info[3] = 0;
		double discrepancy_candidate;
		run_problem(candidate->data, &discrepancy_candidate, info);

//		std::cout << "Sample # " << count << ": "; gsl_vector_fprintf_as_line(stdout, candidate, "%lf"); std::cout << std::endl;

		// Step 2: accept/reject the candidate
		if(accept_mh(candidate, sample, discrepancy_candidate, discrepancy_sample, tolerance))
		{
			//printf("task[%d,%d,%d,%d]: accepted\n", info[0], info[1], info[2], info[3]);
			gsl_vector_memcpy(sample, candidate);
			discrepancy_sample = discrepancy_candidate;
			accept_count += 1;
		}
		else
		{
			//printf("task[%d,%d,%d,%d]: rejected\n", info[0], info[1], info[2], info[3]);
		}


		if(count%every == (every-1))
		{
			db.torc_update(sample->data, discrepancy_sample);
		}
	}

	*acc_rate_array_val = (double)accept_count/(every*chain_length);	// peh: output (OUT,1), stored automatically in the stack of the run_mh_in_parallel routine on node #0

	// free workspace
	gsl_vector_free(candidate);
	gsl_vector_free(sample);

	//printf("Exiting... here (%d,%f,%f,%d)\n", chain_length, discrepancy_seed, tolerance, seed_ind); fflush(0);

}

void abc_subsim::run_mh_in_parallel(int chain_length,
					const size_t * perm,
					//std::vector< gsl_vector * > const& sample_array_prev,
					//const gsl_vector * discrepancy_array_prev,
					std::vector< gsl_vector * > & sample_array_prev,
					gsl_vector * discrepancy_array_prev,
					double tolerance_curr,
					gsl_vector * acc_rate_array,
					int level)
{
//	printf("first = %d, last = %d\n", first, last);
	printf("run_mh_in_parallel: nchains = %d\n", pivot);

	double t0 = torc_gettime();

#if !defined(_PEH_SORTING_)
	for(int number=0; number<pivot; ++number)
	{
		int ind = perm[number];
		double discrepancy_array_prev_ind_val = gsl_vector_get(discrepancy_array_prev, ind);

#if defined(_USE_TORC_)
		torc_create(-1, (void (*)())abc_subsim::run_mh, 7,	// task function with 7 arguments, their description follows:
				1, MPI_INT, CALL_BY_COP,		// chain_length
				PROBDIM, MPI_DOUBLE, CALL_BY_COP,	// sample_array_prev
				1, MPI_DOUBLE, CALL_BY_COP,		// discrepancy_array_prev_ind_val
				1, MPI_DOUBLE, CALL_BY_COP,		// tolerance_curr
				1, MPI_DOUBLE, CALL_BY_RES,		// acc_rate_array
				1, MPI_INT, CALL_BY_COP,		// number
				1, MPI_INT, CALL_BY_COP,		// level
				&chain_length,						// in
				(const double *) (sample_array_prev[ind]->data),	// in: qualified point (perm index)
				&discrepancy_array_prev_ind_val,			// in: qualified discrepancy (perm index)
				&tolerance_curr,					// in
				&acc_rate_array->data[number],				// out: simple index
				&number, &level);
#else
		abc_subsim::run_mh(
				&chain_length,
				(const double *) (sample_array_prev[ind]->data),
				&discrepancy_array_prev_ind_val,
				&tolerance_curr,
				&acc_rate_array->data[number],
				&number, &level);

#endif
	}

#else	// peh:sort  (qualified points according to their sigma value)

	// workspace (temp copy of sample_array_prev + discrepancy_array_prev)
	std::vector< gsl_vector * > sample_array_prev_perm;
	for(int i=0; i<pivot; ++i)
	{
		sample_array_prev_perm.push_back(gsl_vector_alloc(PROBDIM+1));
	}

	// copy selected points into the temporary buffer 
	for(int number=0; number<pivot; ++number)
	{
		int ind = perm[number];
		for (int idx = 0; idx < PROBDIM; idx++) {
			gsl_vector_set(sample_array_prev_perm[number], idx, gsl_vector_get(sample_array_prev[ind], idx));
		}
		gsl_vector_set(sample_array_prev_perm[number], PROBDIM, gsl_vector_get(discrepancy_array_prev, ind));
	}

	// sort them according to sigma
	sort_samples(sample_array_prev_perm);
	//sort_samples_v2(sample_array_prev_perm, torc_num_workers());	// does not help in this case 


	// copy the sorted data back to sample_array_prev + discrepancy_array_prev
	for(int number=0; number<pivot; ++number)
	{
		for (int idx = 0; idx < PROBDIM; idx++) {
			gsl_vector_set(sample_array_prev[number], idx, gsl_vector_get(sample_array_prev_perm[number], idx));
		}
		gsl_vector_set(discrepancy_array_prev, number, gsl_vector_get(sample_array_prev_perm[number], PROBDIM));
		gsl_vector_set(sample_array_prev_perm[number], PROBDIM, gsl_vector_get(discrepancy_array_prev, number));
	}

	// free temporary storage
	for(int i=0; i<pivot; ++i)
	{
		gsl_vector_free(sample_array_prev_perm[i]);
	}

	// submit tasks
	for(int number=0; number<pivot; ++number)
	{
		double discrepancy_array_prev_ind_val = gsl_vector_get(discrepancy_array_prev, number);	// perm -> number

#if defined(_USE_TORC_)
		torc_create((number+0)%torc_num_workers(), (void (*)())abc_subsim::run_mh, 7,	// task function with 6 arguments, their description follows:
//		torc_create(-1, (void (*)())abc_subsim::run_mh, 7,	// task function with 7 arguments, their description follows:
				1, MPI_INT, CALL_BY_COP,		// chain_length
				PROBDIM, MPI_DOUBLE, CALL_BY_COP,	// sample_array_prev
				1, MPI_DOUBLE, CALL_BY_COP,		// discrepancy_array_prev_ind_val
				1, MPI_DOUBLE, CALL_BY_COP,		// tolerance_curr
				1, MPI_DOUBLE, CALL_BY_RES,		// acc_rate_array
				1, MPI_INT, CALL_BY_COP,		// number
				1, MPI_INT, CALL_BY_COP,		// level
				&chain_length,						// in
				(const double *) (sample_array_prev[number]->data),	// in: qualified point (perm index -> number )
				&discrepancy_array_prev_ind_val,			// in: qualified discrepancy (perm index ->  number)
				&tolerance_curr,					// in
				&acc_rate_array->data[number],				// out: simple index
				&number, &level);
#else
		abc_subsim::run_mh(
				&chain_length,
				(const double *) (sample_array_prev[number]->data),
				&discrepancy_array_prev_ind_val,
				&tolerance_curr,
				&acc_rate_array->data[number],
				&number, &level);
#endif
	}

#endif


#if defined(_USE_TORC_)
	torc_enable_stealing();
	torc_waitall();
	torc_disable_stealing();
#endif
	double t1 = torc_gettime();
	printf("run_mh_in_parallel: elapsed time = %lf seconds\n", t1-t0); fflush(0);
}


void abc_subsim::write_array_to_file(std::vector< gsl_vector * > const& sample, const gsl_vector * discrepancy, int level)
{
	printf("%s: level=%d\n", __FUNCTION__, level);
	char filename[256];
	sprintf(filename, "%s_%03d.txt", posterior_filename.c_str(), level);
	FILE * fp = fopen(filename, "w");

	for(int i=0; i<sample.size(); ++i)
	{
		normalize_sample(sample[i]);
		gsl_vector_fprintf_as_line(fp, sample[i], "%.6f");
		fprintf(fp, "%.6f\n", gsl_vector_get(discrepancy, i));
	}

	fclose(fp);
}

int main(int argc, char * argv[])
{
	if(argc < 9)
	{
		std::cout << "Usage: " << argv[0] << " <population_size> <posterior filename> <final tolerance> <max #levels> <acceptance rate> " <<
			"<weights std> <tolerance change> <random seed> [<restart filename>]" << std::endl;
		return -1;
	}

	const char * restart_filename = "";
	int restart = 0;
	if(argc == 10)
	{
		restart = 1;
		restart_filename = argv[9];
	}

	abc_subsim alg(atoi(argv[1]), argv[2], atof(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atoi(argv[8]),
			restart, restart_filename);

	covariance = gsl_matrix_alloc(PROBDIM, PROBDIM);	// peh: declared as global variable (can be static member of abc_subsim)

	fitfun_initialize(NULL);
#if defined(_USE_TORC_)
	torc_register_task((void *)myrandom::set_task);			// peh: task registration for heterogenenous platforms and runtime environments
	torc_register_task((void *)subsim_database::update_task);
	torc_register_task((void *)cov_update_task);
	torc_register_task((void *)run_problem);
	torc_register_task((void *)abc_subsim::run_mh);
	torc_register_task((void *)abc_subsim::tune_covariance_task);

	torc_init(argc, argv, MODE_MW);
#endif

	rng.spmd_set(atoi(argv[8]));	// peh: initialize rng on all nodes

	double t0 = torc_gettime();
	alg.run();
	double t1 = torc_gettime();
	std::cout << "Total elapsed time: " << t1-t0 << " seconds" << std::endl;

	fitfun_finalize();

#if defined(_USE_TORC_)
	torc_finalize();
#endif
//	exit(1);
	return 0;
}
