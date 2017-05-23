/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>


// implementation of your parallel sorting
void parallel_sort(int * begin, int* end, MPI_Comm comm) {


  // my original rank and comm size
  int orank;
  int osize;
  MPI_Comm_size(comm,&osize);
  MPI_Comm_rank(comm,&orank);

  int dbgrank = -1;

  /* initial seed fixed for random numer genrator */
  srand(1000);


  /* duplicate the communicator so we can use it in the loop and free it*/
  MPI_Comm myComm;
  MPI_Comm_dup(comm, &myComm);
  int myrank  = orank;
  int size = osize;

  /* get the count of the array in comm */
  int rcvcount = end-begin;
 
  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort:begin rank="<<myrank<<",  n/p is = "<<rcvcount<<std::endl;

  int *recvbuf;
  if ( rcvcount )
    recvbuf = (int *) malloc(rcvcount * sizeof(int));
  else
    recvbuf = (int *) malloc(1 * sizeof(int));

  for ( int i = 0;i < rcvcount;i++){
    recvbuf[i] = begin[i];
  }

  int loop = 0;

  while (1) { 

    loop++;
    int m;

    /* if only one processor in a group then that processor exits the group */
    MPI_Comm_size(myComm,&size);
    MPI_Comm_rank(myComm,&myrank);

    if (orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", loop="<<loop<<", with size = "<<size<<std::endl;


    /* total number of elements in the group */
    MPI_Allreduce(&rcvcount,&m,1,MPI_INT,MPI_SUM,myComm);

    /* single element in a group or single processor in a group , exit loop */
    if ( ( size  < 2 ) || ( m < 2 ) ) {

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", BREAK - in communicator of in loop="<<loop<<", performing quickiterativesort"<<std::endl;

      /* Sort the elements in the processor , use iterative quicksort */
      if ( rcvcount > 1 )
	quickSortIterative(recvbuf, 0, rcvcount-1);

      break;
    }

    /* pick a pivot */  
    int pivot;

    if ( m > 1 )
      pivot = rand() % (m-1);
    else
      pivot = 0;

    /* processor with rank 0 in current myComm will send it to all processors in the myComm  - as jinx has different random libs on different nodes */
    MPI_Bcast(&pivot,1, MPI_INT, 0, myComm);

    // floor and ceil values
    int fnp = (int)floor(((double)m)/size);
    int cnp = (int)ceil( ((double)m)/size);
    int rem = m % size;

    /* get the rank of the processor with pivot */
    int pivrank;
    if ( m < size ) {
      fnp = 1;
      cnp = 1;
    }
   
    if ( fnp == cnp ) {
      pivrank = pivot/fnp; 
    } else if ( pivot < rem * cnp ) {
      pivrank = pivot/cnp; 
    } else {
      pivrank = rem + (pivot-rem*cnp)/fnp; 
    }
    

    if ( orank == dbgrank )
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", loop="<<loop<<", pivot rank is = "<<pivrank<<", out of ="<<m<<std::endl;

    /* processor with pivot will broadcast its value to all processors in myComm */
    int pivval;
    int lpivot;

    if ( myrank == pivrank ){

      if ( pivrank < rem ) {
	lpivot = pivot - ( pivrank * cnp );
      } else {
	lpivot = pivot - ( pivrank * fnp ) - rem;
      }
      
      pivval = recvbuf[lpivot];

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", pivot value is = "<<pivval<<", out of ="<<m<<std::endl;
    } 

    // broadcast
    MPI_Bcast(&pivval,1, MPI_INT, pivrank, myComm);

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", pivot value from bcast is = "<<pivval<<", out of ="<<m<<std::endl;

    /* each processor divide the list into numbers more than pivot and those less */ 
    int *sendbuf = (int *)malloc(sizeof(int)*rcvcount);

    int mptr[2];
    mptr[0] = -1;
    mptr[1] = rcvcount;

    for ( int i = 0;i < rcvcount;i++){

      if (recvbuf[i] <= pivval ) {

	mptr[0]++;

	if ( orank == dbgrank)
	  std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", recvbuf["<<i<<"]= "<<recvbuf[i]<<", goes to position ="<<mptr[0]<<std::endl;

	sendbuf[mptr[0]]= recvbuf[i];
      } else  {

	mptr[1]--;

	if ( orank == dbgrank)
	  std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", recvbuf["<<i<<"]= "<<recvbuf[i]<<", goes to position ="<<mptr[1]<<std::endl;

	sendbuf[mptr[1]]= recvbuf[i];
      }

    }

    /* keeps count of the number of elements less than pivot and more then pivot */
    mptr[0]++;
    mptr[1] = rcvcount - mptr[1];

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", count <= pivot="<<mptr[0]<<", count > pivot ="<<mptr[1]<<std::endl;

    /*---------------------------------------------------------------------------------------*/
    /* Now decide what should go to what processor and what should come here from where */
    /*---------------------------------------------------------------------------------------*/

    /* collect the count of values less them pivot and more than pivot on each processor */
    int *amptr = (int*)malloc(sizeof(int)*2*size);
    MPI_Allgather(&mptr[0],2,MPI_INT,amptr,2,MPI_INT,myComm);

    /* count the total min and max values */
    int tlp = 0;
    int rtlp = 0;
    int tgp = 0;
    int rtgp = 0;

    /* from the elements received from all processors get the total min/max and total count */
    for ( int i=0;i<size; i++) {

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", recvd from processor "<<i<<",m<pivot="<<amptr[2*i]<<", m>pivot= "<<amptr[2*i+1]<<std::endl;

      /* count the total min and max values before this processor */
      if ( myrank > i) {
	rtlp = rtlp + amptr[2*i];
	rtgp = rtgp + amptr[2*i+1];
      }

      tlp = tlp + amptr[2*i];
      tgp = tgp + amptr[2*i+1];
    }

    /* The total max+min count of all processors in the comm*/
    m = tlp + tgp;

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", sum of lt/gt recvd from all processors, "<<tlp<<"<=pivot<"<<tgp<<std::endl;

    /* approx max size of elements on each processor if evenly distributed */
    cnp = (int)ceil( ((double)m)/size);
    if ( m < size ) 
      cnp = 1;

    /* calculate Number of procesors for values < pivot for this size */
    int lprocessors = tlp / cnp;
    int gprocessors;

    /* make sure we have at least one processor for values in each set */
    if ( lprocessors == size ){
      lprocessors = size-1;
      gprocessors = 1;
    } else if (  lprocessors == 0 ) {
      lprocessors = 1;
      gprocessors = size-1;
    } else {
      gprocessors = size-lprocessors;
    }
   
    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", lprocessors="<<lprocessors<<", gprocessors="<<gprocessors<<std::endl;

    /* determine how many elements to send where and populate the arrays used in alltoall*/
    int *sdispls = (int *)malloc(size*sizeof(int));
    int *sendcounts = (int *)malloc(size*sizeof(int));
    int *rdispls = (int *)malloc(size*sizeof(int));
    int *recvcounts = (int *)malloc(size*sizeof(int));

    int maxCountBeforeCurrentProcessor;
    rcvcount = 0;

    /* the max count before this processor and in this processor allowed */
    if ( myrank < lprocessors ) {
      maxCountBeforeCurrentProcessor = myrank * (int)floor(((double)tlp)/lprocessors) +  (((int)myrank < tlp % lprocessors) ? myrank : tlp % lprocessors );
      rcvcount =  (int)floor(((double)tlp)/lprocessors) +  (((int)myrank < tlp % lprocessors) ? 1 : 0 );
    } else {
      maxCountBeforeCurrentProcessor = (myrank-lprocessors) * (int)floor(((double)tgp)/gprocessors) +  (((int)(myrank-lprocessors) < tgp % gprocessors) ?(myrank-lprocessors)  : tgp % gprocessors );
      rcvcount =  (int)floor(((double)tgp)/gprocessors) +  (((int)(myrank-lprocessors) < tgp % gprocessors) ? 1 : 0 );
    }

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", maxCountBeforeCurrentProcessor="<<maxCountBeforeCurrentProcessor<<", rcvcount="<<rcvcount<<std::endl;

    /* set up the count and disp array for values > pivot val */
    int distributedFromCurrentProcessor;
    int cumulativeDistributedToltgtProcessors[2] = {0,0};
    int totalElementsinCurrentProcessor;
    int totalElementsBeforeCurrentProcesssor;
    int totalElements;
    int nprocessors;
    int pcount;
    int elementsSentBeforei = 0;

    /* for each processor determine the elements to be sent and the elements to receive */
    for (int i=0;i<size;i++) {

      /* initialize the values depending on if processor in current iteration has to be sent values < pivor ot values > pivot */
      if ( i < lprocessors ) {
	totalElementsinCurrentProcessor = mptr[0];
	distributedFromCurrentProcessor = cumulativeDistributedToltgtProcessors[0];
	totalElementsBeforeCurrentProcesssor = rtlp;
	totalElements = tlp;
	nprocessors  = lprocessors;
	pcount = i;
      } else {
	totalElementsinCurrentProcessor = mptr[1];
	distributedFromCurrentProcessor = cumulativeDistributedToltgtProcessors[1];
	totalElementsBeforeCurrentProcesssor = rtgp;
	totalElements = tgp;
	nprocessors = gprocessors;
	pcount = i-lprocessors;
      }

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", totalElementsBeforeCurrentProcesssor="<<totalElementsBeforeCurrentProcesssor<<", totalElementsinCurrentProcessor="<<totalElementsinCurrentProcessor<<",distributedFromCurrentProcessor="<<distributedFromCurrentProcessor<<std::endl;

      /* To send */
      /* does current processor , still have elements to send ? */
      if  ( (totalElementsinCurrentProcessor-distributedFromCurrentProcessor) > 0 ) {

	int maxBeforeiProcessor = pcount * (int)floor(((double)totalElements)/nprocessors) +  ((pcount < totalElements % nprocessors) ? pcount : totalElements % nprocessors );
	int maxOniProcessor =  (int)floor(((double)totalElements)/nprocessors) +  ((pcount < totalElements % nprocessors) ? 1 : 0 );

	if ( orank == dbgrank)
	  std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<" sending to processor, i="<<i<<",maxBeforeiProcessor ="<<maxBeforeiProcessor<<",maxOniProcessor="<<maxOniProcessor<<std::endl;

	/* If i processor will be filled from previous processor and no elements from current processor go there */  
	if ( totalElementsBeforeCurrentProcesssor >= maxBeforeiProcessor + maxOniProcessor ) {
	  sendcounts[i] = 0;
	} else if (totalElementsBeforeCurrentProcesssor + totalElementsinCurrentProcessor >= maxBeforeiProcessor + maxOniProcessor ) { 
	  /* this i processor will be completely filled with elements from current processor */
	  sendcounts[i] = (totalElementsBeforeCurrentProcesssor > maxBeforeiProcessor )? (maxBeforeiProcessor + maxOniProcessor - totalElementsBeforeCurrentProcesssor) :maxOniProcessor;
	} else if ( totalElementsBeforeCurrentProcesssor + totalElementsinCurrentProcessor < maxBeforeiProcessor + maxOniProcessor ) {
	  /* this i processor will only be filled partially by current processor */
	  sendcounts[i] = ( totalElementsBeforeCurrentProcesssor > maxBeforeiProcessor )? totalElementsinCurrentProcessor:(totalElementsinCurrentProcessor-distributedFromCurrentProcessor);
	} else {
	  std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", ERROR during distribution in iteration "<<i<<std::endl;
	}
  
	distributedFromCurrentProcessor = distributedFromCurrentProcessor + sendcounts[i];

      } else {
	sendcounts[i] = 0;
      }

      /* set up the displacement count */
      if ( i > 0 )
	sdispls[i]=sdispls[i-1]+sendcounts[i-1];
      else
	sdispls[i]=0;


      /* initialize the values depending on if processor in current iteration has to be sent values < pivor ot values > pivot */
      if ( i < lprocessors ) {
	cumulativeDistributedToltgtProcessors[0] = distributedFromCurrentProcessor;
      } else {
	cumulativeDistributedToltgtProcessors[1] = distributedFromCurrentProcessor;
      }

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", sendcounts["<<i<<"]="<<sendcounts[i]<<", sdispls["<<i<<"]="<<sdispls[i]<<std::endl;

      /* To recv */
      /* Processors before me should have how many */
      /* Anything more than that till I reach capacity is to be sent to me */
      /* Other processors will have to send me 0 */

      /* elements in current i processor, which this processor needs to receive depends on if this is a <pivot or >pivot value processor */
      int countIniProcessor;
      if ( myrank < lprocessors ) {
	countIniProcessor =  amptr[2*i];
      } else {
	countIniProcessor =  amptr[2*i+1];
      }

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<" receiving from processor, i="<<i<<",countIniProcessor ="<<countIniProcessor<<", elementsSentBeforei="<<elementsSentBeforei<<",rcvcount="<<rcvcount<<std::endl;

      /* processor before i have filled this current processor */
      if ( elementsSentBeforei > maxCountBeforeCurrentProcessor + rcvcount ) {
	recvcounts[i] = 0;
      } else {

	/* all elements up till current i processor will not fill into this processor */
	if ( elementsSentBeforei + countIniProcessor <= maxCountBeforeCurrentProcessor ) {
	  recvcounts[i] = 0;
	} else {
	  /* elements from current i processor can be sent to current processor */
	  /* if elements from current i processor go to the processor before current processor, */
	  /* then determine how many to current processor */
	  if ( elementsSentBeforei <=  maxCountBeforeCurrentProcessor ) {
	    recvcounts[i] = ((elementsSentBeforei + countIniProcessor - maxCountBeforeCurrentProcessor) <= rcvcount)? (elementsSentBeforei + countIniProcessor - maxCountBeforeCurrentProcessor) : rcvcount;
	  } else {
	    recvcounts[i] = ((rcvcount - elementsSentBeforei + maxCountBeforeCurrentProcessor) >= countIniProcessor )? countIniProcessor : (rcvcount - elementsSentBeforei + maxCountBeforeCurrentProcessor);
	  }
	}

      }

      elementsSentBeforei = elementsSentBeforei + countIniProcessor;

      /* set up the displacement count */
      if ( i > 0 )
	rdispls[i]=rdispls[i-1]+recvcounts[i-1];
      else
	rdispls[i]=0;

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", recvcounts["<<i<<"]="<<recvcounts[i]<<", rdispls["<<i<<"]="<<rdispls[i]<<std::endl;


    }


    /* free the array of count of elements in each processor */
    free (amptr);
 
    /* count the total size required for the recv array in the recv count to recv new elements from all processors in communicator*/
    rcvcount= 0;
    for (int i=0;i<size;i++) {
      rcvcount = rcvcount + recvcounts[i];
    }

    /* free the recv buffer once we have already swapped the elements as min and max than pivot into sendbuf , so we can allocate new mem to this buffer */
    free(recvbuf);

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  added up rcvcount is = "<<rcvcount<<", the size from end-begin is ="<<end-begin<<std::endl;

    /* allocate space in the present processor to recieve elements 
       if the processors should not recv any elements then it is a single processor in the group 
       allocate a token memory so the free functions later on do not need a special check
    */
    if ( rcvcount )
      recvbuf = (int *)malloc(rcvcount*sizeof(int));
    else
      recvbuf = (int *)malloc(1*sizeof(int));

    /*DEBUG*/
    if ( orank == dbgrank ){
      int scount = 0;

      for (int i=0;i<size;i++) {
	scount = scount + sendcounts[i];
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  send in alltoallv sendcounts["<<i<<"]="<<sendcounts[i]<<"-"<<sdispls[i]<<std::endl;
      }

      for (int i=0;i<scount;i++) {
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  send buf in alltoallv sendbuf["<<i<<"]="<<sendbuf[i]<<std::endl;
      }

      for (int i=0;i<size;i++) {
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  recvd in alltoallv recvcounts["<<i<<"]="<<recvcounts[i]<<"-"<<rdispls[i]<<std::endl;
      }

    }

    /* send elements from each processors to every other processor */
    MPI_Alltoallv(sendbuf,sendcounts,sdispls,MPI_INT,recvbuf,
		  recvcounts, rdispls, MPI_INT ,
		  myComm);

    if (orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  completed alltoallv in loop="<<loop<<std::endl;

    /*DEBUG*/
    if ( orank == dbgrank )
      for (int i=0;i<rcvcount;i++) {
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  recvd from alltoallv recvbuf["<<i<<"]="<<recvbuf[i]<<std::endl;
      }


    /* free sendbuf and other arrays used in MPI_Alltoallv, as no more needed */
    free(sendbuf);
    free(sendcounts);
    free(sdispls);
    free(recvcounts);
    free(rdispls);

    /* ---------------------------------------------------------------------------*/
    /* wait for all processors to come in */
    /* ---------------------------------------------------------------------------*/
    MPI_Barrier(myComm);

    /* split the communicator into two groups of those <= pivot value and those > pivot val */
    MPI_Comm tempComm;
    MPI_Comm_split(myComm, (myrank<lprocessors)?1:0,myrank, &tempComm);
    MPI_Comm_free(&myComm);
    MPI_Comm_dup(tempComm, &myComm);
    MPI_Comm_free(&tempComm);

  }


  /* ---------------------------------------------------------------------------*/
  /* finished sorting loop , wait for all processors to come in */
  /* ---------------------------------------------------------------------------*/
  MPI_Barrier(comm);

  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<"Completed Loop for All processors"<<std::endl;

  // restore the original rank and size of MPI_WORL_COMM
  myrank = orank;
  size = osize;


  /* ---------------------------------------------------------------------------*/
  /* rebalancing */
  /* ---------------------------------------------------------------------------*/

  /* Now rebalance the elements in each processor to the original count , 
     copy them back to the input memory location */

  /* collect the count of current values and count of expected values on each processor */
  int *amptr = (int*)malloc(sizeof(int)*2*size);
  int mptr[2];
  mptr[0] = rcvcount; // currently held 
  mptr[1] = (end-begin); // expected by results
 
  /* gather this count from each processor in each processor */
  MPI_Allgather(&mptr[0],2,MPI_INT,amptr,2,MPI_INT,comm);

  /* keep count of total elements, actual elements before current processor and expected elements by results before
     current processor */
  int m = 0;
  int mr = 0;
  int mc = 0;

  for ( int i=0;i<size; i++) {

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", recvd from processor "<<i<<", current elements in processor ="<<amptr[2*i]<<" , elements expected by results ="<<amptr[2*i+1]<<std::endl;

    /* count the total values actually held before this processor */
    if ( myrank > i) 
      mr = mr + amptr[2*i];

    m = m + amptr[2*i];
    mc = mc + amptr[2*i+1];

  }


  /* determine how many elements to send where and populate the arrays used in alltoall*/
  int *sdispls = (int *)malloc(size*sizeof(int));
  int *sendcounts = (int *)malloc(size*sizeof(int));
  int *rdispls = (int *)malloc(size*sizeof(int));
  int *recvcounts = (int *)malloc(size*sizeof(int));

  /* the max elements before the current processor, max elements in current proessor */
  int maxCountBeforeCurrentProcessor;
  rcvcount = 0;
    
  maxCountBeforeCurrentProcessor = myrank * (int)floor(((double)m)/size) +  (((int)myrank < m % size) ? myrank : m % size );
  rcvcount =  (int)floor(((double)m)/size) +  (((int)myrank < m % size) ? 1 : 0 );

  if ( orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", rebalancing : maxCountBeforeCurrentProcessor="<<maxCountBeforeCurrentProcessor<<", rcvcount="<<rcvcount<<std::endl;


  /* set up the count and disp array for values > pivot val */
  int distributedFromCurrentProcessor = 0;
  int totalElementsinCurrentProcessor = mptr[0];
  int totalElementsBeforeCurrentProcesssor = mr;
  int maxBeforeiProcessor = 0;
  int maxOniProcessor = 0;
  int elementsSentBeforei = 0;

  /* for each processor determine the elements to be sent and the elements to receive */
  for (int i=0;i<size;i++) {

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort:rank="<<orank<<","<<myrank<<", rebalancing : totalElementsBeforeCurrentProcesssor="<<totalElementsBeforeCurrentProcesssor<<", totalElementsinCurrentProcessor="<<totalElementsinCurrentProcessor<<",distributedFromCurrentProcessor="<<distributedFromCurrentProcessor<<std::endl;

    /* To send */
    /* does current processor , still have elements to send ? */
    if  ( (totalElementsinCurrentProcessor-distributedFromCurrentProcessor) > 0 ) {

      /* get the values from the values sent by the processor which are the actual begin-end */
      maxBeforeiProcessor = ( i > 0 )? (maxBeforeiProcessor+amptr[2*(i-1)+1]):0;
      maxOniProcessor = amptr[2*i+1];

      if ( orank == dbgrank)
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<" rebalancing : sending to processor, i="<<i<<",maxBeforeiProcessor ="<<maxBeforeiProcessor<<",maxOniProcessor="<<maxOniProcessor<<std::endl;

      /* If i processor will be filled from previous processor and no elements from current processor go there */  
      if ( totalElementsBeforeCurrentProcesssor >= maxBeforeiProcessor + maxOniProcessor ) {
	sendcounts[i] = 0;
      } else if (totalElementsBeforeCurrentProcesssor + totalElementsinCurrentProcessor >= maxBeforeiProcessor + maxOniProcessor ) { 
	/* this i processor will be completely filled with elements from current processor */
	sendcounts[i] = (totalElementsBeforeCurrentProcesssor > maxBeforeiProcessor )? (maxBeforeiProcessor + maxOniProcessor - totalElementsBeforeCurrentProcesssor) :maxOniProcessor;
      } else if ( totalElementsBeforeCurrentProcesssor + totalElementsinCurrentProcessor < maxBeforeiProcessor + maxOniProcessor ) {
	/* this i processor will only be filled partially by current processor */
	sendcounts[i] = ( totalElementsBeforeCurrentProcesssor > maxBeforeiProcessor )? totalElementsinCurrentProcessor:(totalElementsinCurrentProcessor-distributedFromCurrentProcessor);
      } else {
	std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", ERROR during rebalancing in iteration "<<i<<std::endl;
      }
  
      distributedFromCurrentProcessor = distributedFromCurrentProcessor + sendcounts[i];

    } else {
      sendcounts[i] = 0;
    }

    /* set up the displacement count */
    if ( i > 0 )
      sdispls[i]=sdispls[i-1]+sendcounts[i-1];
    else
      sdispls[i]=0;


    /* To recv */
    /* Processors before me should have how many */
    /* Anything more than that till I reach capacity is to be sent to me */
    /* Other processors will have to send me 0 */

    /* elements in current i processor, got from the elements sent by the processor in allgather */
    int countIniProcessor = amptr[2*i];

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<" rebalancing: receiving from processor, i="<<i<<",countIniProcessor ="<<countIniProcessor<<", elementsSentBeforei="<<elementsSentBeforei<<",rcvcount="<<rcvcount<<std::endl;

    /* processor before i have filled this current processor */
    if ( elementsSentBeforei > maxCountBeforeCurrentProcessor + rcvcount ) {
      recvcounts[i] = 0;
    } else {

      /* all elements up till current i processor will not fill into this processor */
      if ( elementsSentBeforei + countIniProcessor <= maxCountBeforeCurrentProcessor ) {
	recvcounts[i] = 0;
      } else {
	/* elements from current i processor can be sent to current processor */
	/* if elements from current i processor go to the processor before current processor, then determine how many to current processor */
	if ( elementsSentBeforei <=  maxCountBeforeCurrentProcessor ) {
	  recvcounts[i] = ((elementsSentBeforei + countIniProcessor - maxCountBeforeCurrentProcessor) <= rcvcount)? (elementsSentBeforei + countIniProcessor - maxCountBeforeCurrentProcessor) : rcvcount;
	} else {
	  recvcounts[i] = ((rcvcount - elementsSentBeforei + maxCountBeforeCurrentProcessor) >= countIniProcessor )? countIniProcessor : (rcvcount - elementsSentBeforei + maxCountBeforeCurrentProcessor);
	}
      }

    }

    /* update elements before i processor */
    elementsSentBeforei = elementsSentBeforei + countIniProcessor;

    /* set up the displacement count */
    if ( i > 0 )
      rdispls[i]=rdispls[i-1]+recvcounts[i-1];
    else
      rdispls[i]=0;

    if ( orank == dbgrank)
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", rebalancing : "<<", sendcounts["<<i<<"]="<<sendcounts[i]<<", sdispls["<<i<<"]="<<sdispls[i]<<", recvcounts["<<i<<"]="<<recvcounts[i]<<", rdispls["<<i<<"]="<<rdispls[i]<<std::endl;


  }
  
  /* free the array of count of elements in each processor */
  free (amptr);

  /* count the total size required for the recv array in the recv count */
  rcvcount= 0;
  for (int i=0;i<size;i++) {
    rcvcount = rcvcount + recvcounts[i];
  }
  
  if ( orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", rebalancing : added up rcvcount is = "<<rcvcount<<", the size from begin-end is ="<<end-begin<<std::endl;

  /* allocate space in the present processor to recieve elements */
  int *nrecvbuf;
  if ( rcvcount )
    nrecvbuf = (int *)malloc(rcvcount*sizeof(int));
  else 
    nrecvbuf = (int *)malloc(1*sizeof(int));

  //DEBUG
  if ( orank == dbgrank ){
    int scount = 0;

    for (int i=0;i<size;i++) {
      scount = scount + sendcounts[i];
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", rebalancing : send in alltoallv sendcounts["<<i<<"]="<<sendcounts[i]<<"-"<<sdispls[i]<<std::endl;
    }

    for (int i=0;i<scount;i++) {
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<",  rebalancing : send buf in alltoallv sendbuf["<<i<<"]="<<recvbuf[i]<<std::endl;
    }

    for (int i=0;i<size;i++) {
      std::cout<<"DEBUG:parallel_sort::parallel_sort: rank="<<orank<<","<<myrank<<", rebalancing : recvd in alltoallv recvcounts["<<i<<"]="<<recvcounts[i]<<"-"<<rdispls[i]<<std::endl;
    }

  }


  /* send elements from each processors to every other processor for balancing */
  MPI_Alltoallv(recvbuf, sendcounts,sdispls, MPI_INT , nrecvbuf,
		recvcounts, rdispls, MPI_INT ,
		comm);

  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort:rank="<<orank<<", completed rebalancing alltoallv"<<std::endl;

  /* ---------------------------------------------------------------------------*/
  /* wait for all processors to come in */
  /* ---------------------------------------------------------------------------*/
  MPI_Barrier(comm);

  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort:rank="<<orank<<", freeing memory "<<std::endl;

  /* free sendbuf and other arrays used in MPI_Alltoallv, as no more needed */
  free(recvbuf);
  free(sendcounts);
  free(sdispls);
  free(recvcounts);
  free(rdispls);

  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort:rank="<<orank<<", before swapping elements"<<std::endl;

  /* swap the elements back to the begin pointer */
  int np = end-begin;
  for (int i=0;i<np;i++) {
    begin[i] = nrecvbuf[i];
  }

  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort:rank="<<orank<<", completed swapping elements="<<np<<std::endl;

  /* free the recv buffer once we have already swapped the elements as min and max than pivot into sendbuf */
  free (nrecvbuf);
  MPI_Comm_free(&myComm);

  if (orank == dbgrank)
    std::cout<<"DEBUG:parallel_sort::parallel_sort:rank="<<orank<<", completed parallel_sort routine"<<std::endl;

}















/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/
/*
  I have used the code in these iterative quick sort routines with some modifications 
  from http://www.geeksforgeeks.org
  as the programming instructions permitted us to use it for local sorting only
*/

/* swap two elements at the two addresses passed
 */
void swap ( int* a, int* b )
{
  int t = *a;
  *a = *b;
  *b = t;
}
 
/* Partion the array in place, to values less than pivot and more than pivot 
   return the pivot position 
*/
int partition (int arr[], int l, int h)
{
  int x = arr[h];
  int i = (l - 1);
 
  for (int j = l; j <= h- 1; j++)
    {
      if (arr[j] <= x)
        {
	  i++;
	  swap (&arr[i], &arr[j]);
        }
    }
  swap (&arr[i + 1], &arr[h]);
  return (i + 1);
}
 
/* 
   Quick sort routine, implemented iterative using a stack
   arr --> Array to be sorted, 
   l   --> Starting index, 
   h  --> Ending index */
void quickSortIterative (int arr[], int l, int h)
{
  // Create an auxiliary stack
  int stack[ h - l + 1 ];
 
  // initialize top of stack
  int top = -1;
 
  // push initial values of l and h to stack
  stack[ ++top ] = l;
  stack[ ++top ] = h;
 
  // Keep popping from stack while is not empty
  while ( top >= 0 )
    {
      // Pop h and l
      h = stack[ top-- ];
      l = stack[ top-- ];
 
      // Set pivot element at its correct position in sorted array
      int p = partition( arr, l, h );
 
      // If there are elements on left side of pivot, then push left
      // side to stack
      if ( p-1 > l )
        {
	  stack[ ++top ] = l;
	  stack[ ++top ] = p - 1;
        }
 
      // If there are elements on right side of pivot, then push right
      // side to stack
      if ( p+1 < h )
        {
	  stack[ ++top ] = p + 1;
	  stack[ ++top ] = h;
        }
    }
}


// ...
