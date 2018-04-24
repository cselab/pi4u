/*
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _TMCMC_DB_H_
#define _TMCMC_DB_H_



void update_full_db(double point[], double F, double *G, int n, int surrogate);
void init_full_db();

void update_curgen_db(double point[], double F, double prior);
void init_curgen_db();

void update_curres_db(double point[], double F);
void init_curres_db();
void print_full_db();
void print_curgen_db();
void dump_curgen_db(int Gen);
void dump_curres_db(int Gen);
void dump_full_db(int Gen);
void display_curgen_db(int Gen);
int load_curgen_db(int Gen);



#endif
