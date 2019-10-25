/*

Ubuntu:
-------
gcc -I/home/heiland/dev/libRR_1.3/include/rr/C -fopenmp  omp_rr.c -L/home/heiland/dev/libRR_1.3/lib -lroadrunner_c_api -o omp_rr

OSX:
-----
clang -Xpreprocessor -fopenmp -m64  -I/Users/heiland/dev/roadrunner-osx-10.9-cp36m/include/rr/C -L/Users/heiland/dev/roadrunner-osx-10.9-cp36m/lib -lroadrunner_c_api  -L/usr/local/opt/libomp/lib -lomp  omp_rr.c -o omp_rr
export DYLD_LIBRARY_PATH=/Users/heiland/dev/roadrunner-osx-10.9-cp36m/lib
http://sys-bio.github.io/roadrunner/c_api_docs/html/group__helper_routines.html
http://sys-bio.github.io/roadrunner/c_api_docs/html/rrc__api_8h.html
http://sys-bio.github.io/roadrunner/c_api_docs/html/group__initial_conditions.html


clang -I/Users/heiland/dev/roadrunner-osx-10.9-cp36m/include/rr/C -L/Users/heiland/dev/roadrunner-osx-10.9-cp36m/lib -lroadrunner_c_api  -L/usr/local/opt/libomp/lib -lomp  omp_rr.c -o omp_rr

*/
#undef __cplusplus
#define STATIC_RRC
#include <stdio.h>
#include <stdlib.h>
#include "rrc_api.h"
#include "rrc_types.h"
#include "rrc_utilities.h"

#include <omp.h>

int main (int argc, char *argv[]) {
   /* RRHandle rrHandle;
    RRCDataPtr result;
   int index;
   int col;
   int row;
   */

   if(argc < 4) {
      printf("Provide args: <name of sbml file> <num cells> <num threads>\n");
      exit(1);
   }
  int ncells = atoi(argv[2]);
  printf("------  # of cells = %d\n\n", ncells);
  int nthreads = atoi(argv[3]);
  printf("------  # of threads = %d\n\n", nthreads);

  omp_set_num_threads(nthreads);

  RRHandle rrHandle;
  RRHandle rrHandleArray[42];
  RRVectorPtr vptr;
  RRCDataPtr result;  // start time, end time, and number of points
  int idx_oxygen;

//   #pragma omp parallel for
  #pragma omp parallel for private(rrHandle, vptr, result, idx_oxygen)
  for (int icell=1; icell <= ncells; icell++)
  {
   idx_oxygen = 3;
   printf ("------------  cell %d  ---------------\n", icell);

   printf ("Starting Test Program %s\n", argv[0]);
   // RRHandle rrHandle = createRRInstance();
   rrHandle = createRRInstance();
   if (!loadSBML (rrHandle, "feedback.xml")) {
      printf ("Error while loading SBML file\n");
      printf ("Error message: %s\n", getLastError());
      getchar ();
      exit (0);
   }
   // printf("API = %s\n", getAPIVersion ());
   // printf("getNumberOfIndependentSpecies= %d\n", getNumberOfIndependentSpecies(rrHandle));
   // printf("getNumInstantiatedIntegrators= %d\n", getNumInstantiatedIntegrators(rrHandle));

   // RRVector* v= getFloatingSpeciesInitialConcentrations(rrHandle);
   // struct RRVector vec;
   // setFloatingSpeciesConcentrations(rrHandle, &vec);
   // RRStringArrayPtr = getFloatingSpeciesInitialConditionIds(rrHandle);
   // double myvec[9];
   // char* ids;
   // ids = getFloatingSpeciesInitialConditionIds(rrHandle);
   // printf("getFloatingSpeciesInitialConditionIds = %s\n",ids);

   int r = getNumberOfReactions(rrHandle);
   int m = getNumberOfFloatingSpecies(rrHandle);
   int b = getNumberOfBoundarySpecies(rrHandle);
   int p = getNumberOfGlobalParameters(rrHandle);
   int c = getNumberOfCompartments(rrHandle);

   printf ("Number of reactions = %d\n", r);
   printf ("Number of floating species = %d\n", m);  // 4
   printf ("Number of boundary species = %d\n", b);  // 0
   printf ("Number of compartments = %d\n", c);  // 1

   printf ("Floating species names:\n");
   printf ("-----------------------\n");
   // cout<<stringArrayToString(getFloatingSpeciesIds())<<endl<<endl;
   printf("%s\n\n",stringArrayToString(getFloatingSpeciesIds(rrHandle)));

   // What is this? How does it differ from the one above??
   // printf ("Initial Floating species names:\n");
   // printf ("-------------------------------\n");
   // // cout<<stringArrayToString(getFloatingSpeciesInitialConditionIds())<<endl;
   // printf("%s\n\n",stringArrayToString(getFloatingSpeciesInitialConditionIds(rrHandle)));

/*   rrc_api.h:C_DECL_SPEC RRVectorPtr rrcCallConv getFloatingSpeciesConcentrations(RRHandle handle); 
rrc_api.h:C_DECL_SPEC bool rrcCallConv setFloatingSpeciesInitialConcentrations (RRHandle handle, const RRVectorPtr vec);
*/
   printf ("Floating species conc:\n");
   printf ("-------------------------------\n");
   // RRVectorPtr vptr = getFloatingSpeciesConcentrations(rrHandle);
   vptr = getFloatingSpeciesConcentrations(rrHandle);

/*
rrc_types.h:
typedef struct RRVector
{
    int             Count;
    double*         Data; 
} *RRVectorPtr; 
*/
   printf("%d\n",vptr->Count);
   for (int idx=0; idx<vptr->Count; idx++)
      printf("%d %f\n",idx, vptr->Data[idx]);

   // idx_oxygen = 3;
   vptr->Data[idx_oxygen] += 0.1;
   setFloatingSpeciesConcentrations(rrHandle, vptr);

   vptr = getFloatingSpeciesConcentrations(rrHandle);
   printf("%d\n",vptr->Count);
   for (int idx=0; idx<vptr->Count; idx++)
      printf("%d %f\n",idx, vptr->Data[idx]);

   result = simulateEx (rrHandle, 0, 10, 10);  // start time, end time, and number of points
   // int index = 0;
   // // Print out column headers... typically time and species.
   // for (int col = 0; col < result->CSize; col++)
   // {
   //    printf ("%10s", result->ColumnHeaders[index++]);
   //    if (col < result->CSize - 1)
   //    {
   //       printf ("\t");
   //    }
   // }
   // printf ("\n");

   // index = 0;
   // // Print out the data
   // for (int row = 0; row < result->RSize; row++)
   // {
   //    for (int col = 0; col < result->CSize; col++)
   //    {
   //       printf ("%10f", result->Data[index++]);
   //       if (col < result->CSize -1)
   //       {
   //          printf ("\t");
   //       }
   //    }
   // printf ("\n");
   // }

   //Cleanup
   freeRRCData (result);
   freeRRInstance (rrHandle);
   // getchar ();

  }  // end for
//  }  // end pragma

   exit (0);
}
