/**
 * @file    parallel_sort.h
 * @brief   Declares the parallel sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef PARALLEL_SORT_H
#define PARALLEL_SORT_H

#include <mpi.h>

/**
 * @brief   Parallel, distributed sorting over all processors in `comm`. Each
 *          processor has the local input [begin, end).
 *
 * Note that `end` is given one element beyond the input. This corresponds to
 * the API of C++ std::sort! You can get the size of the local input with:
 * int local_size = end - begin;
 *
 * @param begin Pointer to the first element in the input sequence.
 * @param end   Pointer to one element past the input sequence. Don't access this!
 * @param comm  The MPI communicator with the processors participating in the
 *              sorting.
 */
void parallel_sort(int * begin, int* end, MPI_Comm comm);


/*********************************************************************
 *              Declare your own helper functions here               *
 *********************************************************************/
 
/** 
 * @brief  Quick sort routine, implemented iterative using a stack
 *         to local sort the array in the processor
 * @param arr Array to be sorted, 
 * @param l   Starting index, 
 * @param h   Ending index 
 *
 */
void quickSortIterative (int arr[], int l, int h);


/**
 * @brief swap two elements at the two addresses passed, support routine for 
 *           quick sort
 * 
 * @param a address of element to be swapped
 * @param b address of element to be swapped
 */
void swap ( int* a, int* b );

/**
 * @brief  Partion the array in place, to values less than pivot and more than pivot 
 *         return the pivot position, support routine for quick sort iterative
 *
 * @param  arr Array to be partitioned
 * @param l   Starting index, 
 * @param h   Ending index 
 *
 *  returns the pivot index
 *
*/
int partition (int arr[], int l, int h);

#endif // PARALLEL_SORT_H
