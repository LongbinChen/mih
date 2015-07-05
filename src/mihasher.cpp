#include <algorithm>
#include "mihasher.h"
#include "myhdf5.h"
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "io.h"
#include <sys/stat.h> 
#include <fcntl.h>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/*
 * Inputs: query, numq, dim1queries
 */

/*
 * Outputs: results, numres, stats
 *
 *   results: an array of indices (1-based) of the K-nearest neighbors
 *   for each query. So the array includes K*numq uint32 integers.
 *
 *   numres: includes the number of database entries that fall at any
 *   specific Hamming distance from the query until the K nearest
 *   neighbors are reached. So from this array you can figure out the
 *   Hamming distances of the K-nearest neighbors.
 */

void mihasher::batchquery(UINT32 *results, UINT32 *numres, qstat *stats, UINT8 *queries, UINT32 numq, int dim1queries)
{
    counter = new bitarray;
    counter->init(N);

    UINT32 *res  = new UINT32[K*(D+1)];
    UINT64 *chunks = new UINT64[m];

    UINT32 *presults = results;
    UINT32 *pnumres = numres;
    qstat *pstats = stats;
    UINT8 *pq = queries;

    for (int i=0; i<numq; i++) {
	query(presults, pnumres, pstats, pq, chunks, res);

	presults += K;
	pnumres += B+1;
	pstats ++;
	pq += dim1queries;
    }
	
    delete [] res;
    delete [] chunks;

    delete counter;
}


// Temp variables: chunks, res -- I did not want to malloc inside
// query, so these arrays are passed from outside

void mihasher::query(UINT32 *results, UINT32* numres, qstat *stats, UINT8 *query, UINT64 *chunks, UINT32 *res)
{
    UINT32 maxres = K ? K : N;			// if K == 0 that means we want everything to be processed.
						// So maxres = N in that case. Otherwise K limits the results processed.

    UINT32 n = 0; 				// number of results so far obtained (up to a distance of s per chunk)
    UINT32 nc = 0;				// number of candidates tested with full codes (not counting duplicates)
    UINT32 nd = 0;                      	// counting everything retrieved (duplicates are counted multiple times)
    UINT32 nl = 0;				// number of lookups (and xors)
    UINT32 *arr;
    int size = 0;
    UINT32 index;
    int hammd;
    clock_t start, end;

    start = clock();

    counter->erase();
    memset(numres, 0, (B+1)*sizeof(*numres));

    split(chunks, query, m, mplus, b);
    
    int s;			// the growing search radius per substring

    int curb = b;		// current b: for the first mplus substrings it is b, for the rest it is (b-1)

    for (s = 0; s <= d && n < maxres; s++) {
	for (int k=0; k<m; k++) {
	    if (k < mplus)
		curb = b;
	    else
		curb = b-1;
	    UINT64 chunksk = chunks[k];
	    nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s

	    UINT64 bitstr = 0; 			// the bit-string with s number of 1s
	    for (int i=0; i<s; i++)
	    	power[i] = i;			// power[i] stores the location of the i'th 1
	    power[s] = curb+1;			// used for stopping criterion (location of (s+1)th 1)

	    int bit = s-1;			// bit determines the 1 that should be moving to the left
	    // we start from the left-most 1, and move it to the left until it touches another one

	    while (true) {			// the loop for changing bitstr
	    	if (bit != -1) {
	    	    bitstr ^= (power[bit] == bit) ? (UINT64)1 << power[bit] : (UINT64)3 << (power[bit]-1);
	    	    power[bit]++;
	    	    bit--;
	    	} else { // bit == -1
	    	    /* the binary code bitstr is available for processing */
	    	    arr = H[k].query(chunksk ^ bitstr, &size); // lookup
	    	    if (size) {			// the corresponding bucket is not empty
	    		nd += size;
	    		for (int c = 0; c < size; c++) {
	    		    index = arr[c];
	    		    if (!counter->get(index)) { // if it is not a duplicate
	    			counter->set(index);
	    			hammd = match(codes + (UINT64)index*(B_over_8), query, B_over_8);
	    			nc++;
	    			if (hammd <= D && numres[hammd] < maxres) {
	    			    res[hammd * K + numres[hammd]] = index + 1;
	    			}
				numres[hammd]++;
	    		    }
	    		}
	    	    }
	    	    /* end of processing */

	    	    while (++bit < s && power[bit] == power[bit+1]-1) {
	    		bitstr ^= (UINT64)1 << (power[bit]-1);
	    		power[bit] = bit;
	    	    }
	    	    if (bit == s)
	    		break;
	    	}
	    }

	    n = n + numres[s*m+k]; // This line is very tricky ;)
	    // The k'th substring (0 based) is the last chance of an
	    // item at a Hamming distance of s*m+k to be
	    // found. Because if until the k'th substring, an item
	    // with distance of s*m+k is not found, then it means that
	    // all of the substrings so far have a distance of (s+1)
	    // or more, and the remaining substrings have a distance
	    // of s or more (total > s*m+k).
	    
	    if (n >= maxres)
		break;
	}
    }
	    
    end = clock();

    stats->ticks = end-start;
    stats->numcand = nc;
    stats->numdups = nd;
    stats->numlookups = nl;

    n = 0;
    for (s = 0; s <= D && n < K; s++ ) {
	for (int c = 0; c < numres[s] && n < K; c++)
	    results[n++] = res[s*K + c];
    }

    UINT32 total = 0;
    stats->maxrho = -1;
    for (int i=0; i<=B; i++) {
	total += numres[i];
	if (total >= K && stats->maxrho == -1)
	    stats->maxrho = i;
    }
    stats->numres = n;
}

mihasher::mihasher(int _B, int _m)
{
    B = _B;
    B_over_8 = B/8;
    m = _m;
    b = ceil((double)B/m);
 
    D = ceil(B/2.0);		// assuming that B/2 is large enough radius to include all of the k nearest neighbors
    d = ceil((double)D/m);
   
    mplus = B - m * (b-1);
    // mplus     is the number of chunks with b bits
    // (m-mplus) is the number of chunks with (b-1) bits

    xornum = new UINT32 [d+2];
    xornum[0] = 0;
    for (int i=0; i<=d; i++)
	xornum[i+1] = xornum[i] + choose(b, i);
    
    H = new SparseHashtable[m];
    // H[i].init might fail
    for (int i=0; i<mplus; i++)
	H[i].init(b);
    for (int i=mplus; i<m; i++)
	H[i].init(b-1);
}

void mihasher::setK(int _K)
{
    K = _K;
}

mihasher::~mihasher()
{
    delete[] xornum;
    delete[] H;
}

void mihasher::populate(UINT8 *_codes, UINT32 _N, int dim1codes)
{
    N = _N;
    codes = _codes;

    int k = 0;
//#pragma omp parallel shared(k)
    {
	UINT64 * chunks = new UINT64[m];
//#pragma omp for
	for (k=0; k<m; k++) {
	    UINT8 * pcodes = codes;

	    for (UINT64 i=0; i<N; i++) {
		split(chunks, pcodes, m, mplus, b);
		
		H[k].count_insert(chunks[k], i);

		if (i % (int)ceil(N/1000) == 0) {
		    printf("%.2f%%\r", (double)i/N * 100);
		    fflush(stdout);
		}
		pcodes += dim1codes;
	    }

	    // for (int k=0; k<m; k++)
	    // 	H[k].allocate_mem_based_on_counts();
	    
	    pcodes = codes;
	    for (UINT64 i=0; i<N; i++) {
		split(chunks, pcodes, m, mplus, b);

		H[k].data_insert(chunks[k], i);

		if (i % (int)ceil(N/1000) == 0) {
		    printf("%.2f%%\r", (double)i/N * 100);
		    fflush(stdout);
		}
		pcodes += dim1codes;
	    }
	}
	 
	delete [] chunks;
    }

    // N = _N;
    // codes = _codes;
    // UINT64 * chunks = new UINT64[m];

    // UINT8 * pcodes = codes;
    // for (UINT64 i=0; i<N; i++, pcodes += dim1codes) {
    // 	split(chunks, pcodes, m, mplus, b);

    // 	for (int k=0; k<m; k++)
    // 	    H[k].lazy_insert(chunks[k], i);
    // 	if (i % (int)ceil(N/1000) == 0) {
    // 	    printf("%.2f%%\r", (double)i/N * 100);
    // 	    fflush(stdout);
    // 	}
    // }

    // for (int k=0; k<m; k++)
    // 	H[k].cleanup_insert(codes, m, k, mplus, b, dim1codes);
    // delete [] chunks;
}

void mihasher::build_index(const char* infile, const char* outfile) {

    UINT32 NQ = 0, Q0 = 0, Q1 = 0;


	/* Loading the codes and queries from the input file */	
    UINT8 *codes_db;
    int dim1codes;
    UINT8 *codes_query;
    
    printf("Loading codes... ");
    fflush(stdout);

    int fd = open(infile, O_RDONLY, (mode_t)0400);
    if (fd == -1)
      return ;
    off_t data_size = lseek(fd, 0, SEEK_END);

    N = data_size / B;
    B_over_8 = B / 8;

    codes_db = (UINT8*)malloc((size_t)N * (B/8) * sizeof(UINT8));
    load_bin_codes(infile, "B");
    dim1codes = B/8;

    printf("done.\n");
   
    fflush(stdout);

        printf("done.\n");
    /* Done with the inputs */

    printf("N = %.0e |", (double)N);
    printf(" NQ = %d, range [%d %d) |", NQ, Q0, Q1);
    printf(" B = %d |", B);
    printf(" K = %4d |", K);
    printf(" m = %2d |", m);
    printf("\n");
  /* Run multi-index hashing for K-nearest neighbor search and store the required stats */
    mihasher *MIH = NULL;
    clock_t start0, end0;
    time_t start1, end1;
    qstat *stats = (qstat*) new qstat[NQ];

    result_t result;
    result.n = N;
    result.nq = NQ;
    result.k = K;
    result.b = B;
    result.m = m;
    result.q0 = Q0;
    result.q1 = Q1;
    result.wt = -1;
    result.cput = -1;
    result.vm = -1;
    result.rss = -1;
    result.res = NULL;
    result.nres = NULL;
    result.stats = NULL;

    result.res = (UINT32 **) malloc(sizeof(UINT32*)*NQ);
    result.res[0] = (UINT32 *) malloc(sizeof(UINT32)*K*NQ);
    for (size_t i=1; i<NQ; i++)
	result.res[i] = result.res[i-1] + K;

    result.nres = (UINT32 **) malloc(sizeof(UINT32*)*NQ);
    result.nres[0] = (UINT32 *) malloc(sizeof(UINT32)*(B+1)*NQ);
    for (size_t i=1; i<NQ; i++)
	result.nres[i] = result.nres[i-1] + (B+1);

    result.stats = (double **) malloc(sizeof(double*)*NQ);
    result.stats[0] = (double *) malloc(sizeof(double)*STAT_DIM*NQ);
    for (size_t i=1; i<NQ; i++)
	result.stats[i] = result.stats[i-1] + STAT_DIM;
    
   
    printf("Populating %d hashtables with %.0e entries...\n", m, (double)N);
    fflush (stdout);
    start1 = time(NULL);
    start0 = clock();
    
    populate(codes_db, N, dim1codes);
	    
    end0 = clock();
    end1 = time(NULL);
    
    double ct = (double)(end0-start0) / (CLOCKS_PER_SEC);
    double wt = (double)(end1-start1);
    
    printf("done. | cpu %.0fm%.0fs | wall %.0fm%.0fs\n", ct/60, ct-60*int(ct/60), wt/60, wt-60*int(wt/60));

}

void mihasher::query() {
	UINT32 NQ = 0;
 	clock_t start0, end0;
    time_t start1, end1;
   	result_t result;
    result.n = N;
    result.nq = NQ;
    result.k = K;
    result.b = B;
    result.m = m;
    result.wt = -1;
    result.cput = -1;
    result.vm = -1;
    result.rss = -1;
    result.res = NULL;
    result.nres = NULL;
    result.stats = NULL;



    printf("query... ");
    fflush (stdout);
    
    start1 = time(NULL);
    start0 = clock();
    
//    MIH->batchquery(result.res[0], result.nres[0], stats, codes_query, NQ, dim1queries);
    
    end0 = clock();
    end1 = time(NULL);
    
    result.cput = (double)(end0-start0) / (CLOCKS_PER_SEC) / NQ;
    result.wt = (double)(end1-start1) / NQ;
    process_mem_usage(&result.vm, &result.rss);
    result.vm  /= double(1024*1024);
    result.rss /= double(1024*1024);
    printf("done | cpu %.3fs | wall %.3fs | VM %.1fgb | RSS %.1fgb     \n", result.cput, result.wt, result.vm, result.rss);

    double *pstats_d = result.stats[0];
 
    return ;

}

void mihasher::load_bin_codes(const char *filename, const char *varStr) {
    hid_t       file, space, dset;          /* Handles */
    herr_t      status;
    hsize_t     dims[10];

    /* Open file and dataset using the default properties */
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, varStr, H5P_DEFAULT);
    space = H5Dget_space(dset);

    /* Get info of the size of the array */
    int ndims = H5Sget_simple_extent_dims(space, dims, NULL);

    printf("dims is %d by %d\n", dims[0], dims[1]);
     fflush (stdout);
    assert(ndims == 2);
    
    hsize_t start[2] = {0, 0};
    hsize_t count[2];
    N = dims[0];
    count[0] = N;
	B = dims[1];
    count[1] = B;
    status = H5Sselect_hyperslab(space, H5S_SELECT_SET, start, NULL, count, NULL);
    printf("status %d\n", status);

    int total_size = (size_t)N * (B/8) * sizeof(UINT8);

    printf("total size is %d\n", total_size);
    fflush(stdout);
    codes = (UINT8*)malloc((size_t)N * (B/8) * sizeof(UINT8));
 
    /*
     * Read the data using the default properties.
     */
    printf("start0: %d, B: %d, N: %d\n", 0, B, N);
    // codes = codes-sizeof(*codes)*start0*(*B); // change the pointer so that it points to the beginning of the selected set.
    status = H5Dread(dset, H5T_NATIVE_UCHAR, H5S_ALL, space, H5P_DEFAULT, codes);
    printf("status %d\n", status);

    /*
     * Close and release resources.
     */
    status = H5Sclose(space);
    status = H5Dclose(dset);
    status = H5Fclose(file);

    return;
}

void mihasher::load_double_matrix(const char *filename, const char *varStr, UINT8 *matrix, int *nrow, int *ncol) {
    hid_t       file, space, dset;          /* Handles */
    herr_t      status;
    hsize_t     dims[10];

    /* Open file and dataset using the default properties */
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen(file, varStr, H5P_DEFAULT);
    space = H5Dget_space(dset);

    /* Get info of the size of the array */
    int ndims = H5Sget_simple_extent_dims(space, dims, NULL);
    assert(ndims == 2);
    
    hsize_t count[2];
    if (*nrow == 0)
	*nrow = dims[0];
    else
	assert(*nrow <= dims[0]);
    if (*ncol == 0)
	*ncol = dims[1];
    else
	assert(*ncol == dims[1]);

    /*
     * Read the data using the default properties.
     */
    // matrix = matrix-sizeof(*matrix)*start0*B;
    status = H5Dread(dset, H5T_NATIVE_UCHAR, H5S_ALL, space, H5P_DEFAULT, matrix);

    /*
     * Close and release resources.
     */
    status = H5Sclose(space);
    status = H5Dclose(dset);
    status = H5Fclose(file);
}

void mihasher::process_mem_usage(double *vm_usage, double *resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   *vm_usage     = 0.0;
   *resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   *vm_usage     = vsize / 1024.0;
   *resident_set = rss * page_size_kb;
}



void mihasher::saveVarRef(hobj_ref_t *ref, hid_t file, int i, int dim1, int dim2, const char *varStrMain, const char *varStr, void *var, const char *type) {
    hsize_t dims[2] = {dim1, dim2};
    hid_t  space = H5Screate_simple(2, dims, NULL);
    herr_t status;
    char str[80];

    sprintf(str, "/refs/%s%d.%s", varStrMain, i, varStr);

    printf("Writing %s with size [%d, %d]\n", str, dim1, dim2);

    hid_t dset = 0;
    if (!strcmp(type, "uint32")) {
	dset = H5Dcreate(file, str, H5T_NATIVE_UINT32, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, (UINT32*)var);
    } else if (!strcmp(type, "double")) {
	dset = H5Dcreate(file, str, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (double*)var);
    }
    status = H5Sclose(space);
    status = H5Dclose(dset);
    
    status = H5Rcreate(ref, file, str, H5R_OBJECT, -1);
}

void mihasher::saveRes(const char *filename, const char *varStr, const result_t *result, int n, int overwrite /* = 1 */) {  
    hid_t       file, space, dset, group;          /* Handles */
    herr_t      status;
    hsize_t     dims1[1], maxdims1[1], chunkdims1[1];
    char str[80];

    if (access(filename, F_OK) == -1) {
	overwrite = 1;
    }

    /* Open file and dataset using the default properties */
    if (overwrite == 1) {
	printf("Creating the file %s\n", filename);
    	file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Creating a group for references */
	group = H5Gcreate(file, "/refs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(group);
    } else if (overwrite == 2) {
    	file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }

    /* Creating the result hdf5 compound type */
    hid_t restype = H5Tcreate(H5T_COMPOUND, sizeof(resulth5_t));
    status = H5Tinsert (restype, "wt",    HOFFSET(resulth5_t, wt),      H5T_NATIVE_DOUBLE);
    status = H5Tinsert (restype, "cput",  HOFFSET(resulth5_t, cput),    H5T_NATIVE_DOUBLE);
    status = H5Tinsert (restype, "vm",    HOFFSET(resulth5_t, vm),      H5T_NATIVE_DOUBLE);
    status = H5Tinsert (restype, "rss",   HOFFSET(resulth5_t, rss),     H5T_NATIVE_DOUBLE);
    status = H5Tinsert (restype, "n",     HOFFSET(resulth5_t, n),       H5T_NATIVE_INT);
    status = H5Tinsert (restype, "m",     HOFFSET(resulth5_t, m),       H5T_NATIVE_INT);
    status = H5Tinsert (restype, "b",     HOFFSET(resulth5_t, b),       H5T_NATIVE_INT);
    status = H5Tinsert (restype, "nq",    HOFFSET(resulth5_t, nq),      H5T_NATIVE_INT);
    status = H5Tinsert (restype, "q0",    HOFFSET(resulth5_t, q0),      H5T_NATIVE_INT);
    status = H5Tinsert (restype, "q1",    HOFFSET(resulth5_t, q1),      H5T_NATIVE_INT);
    status = H5Tinsert (restype, "k",     HOFFSET(resulth5_t, k),       H5T_NATIVE_INT);
    status = H5Tinsert (restype, "res",   HOFFSET(resulth5_t, res),     H5T_STD_REF_OBJ);
    status = H5Tinsert (restype, "nres",  HOFFSET(resulth5_t, nres),    H5T_STD_REF_OBJ);
    status = H5Tinsert (restype, "stats", HOFFSET(resulth5_t, stats),   H5T_STD_REF_OBJ);

    if (overwrite == 1) {
	/* Copying results into hdf5 result type, and saving the references */
	resulth5_t *result5t = new resulth5_t[n];
	for (int i=0; i<n; i++) {
	    memcpy(&result5t[i], &result[i], HOFFSET(resulth5_t, res));
	    saveVarRef(&result5t[i].res,  file, i, result[i].nq, result[i].k, varStr, "res", result[i].res[0], "uint32");
	    saveVarRef(&result5t[i].nres, file, i, result[i].nq, result[i].b+1, varStr, "nres", result[i].nres[0], "uint32");
	    if (result[i].stats != NULL)
		saveVarRef(&result5t[i].stats, file, i, result[i].nq, STAT_DIM, varStr, "stats", result[i].stats[0], "double");
	    else
		result5t[i].stats = 0;
	}
	dims1[0] = n;
	maxdims1[0] = H5S_UNLIMITED;
	chunkdims1[0] = 1;

	/* Modify dataset creation properties, i.e. enable chunking  */
	hid_t prop = H5Pcreate(H5P_DATASET_CREATE);
	status = H5Pset_chunk(prop, 1, chunkdims1);
	
	space = H5Screate_simple (1, dims1, maxdims1);
	dset = H5Dcreate(file, varStr, restype, space, H5P_DEFAULT, prop, H5P_DEFAULT);
	status = H5Dwrite(dset, restype, H5S_ALL, H5S_ALL, H5P_DEFAULT, result5t);

	free(result5t);

	status = H5Pclose(prop);
	// status = H5Sclose(space);
	status = H5Dclose(dset);
    } else if (overwrite == 2) {
	dset = H5Dopen(file, varStr, H5P_DEFAULT);

	space = H5Dget_space(dset);
	int ndims = H5Sget_simple_extent_dims(space, dims1, NULL);
	assert(ndims == 1);
	int nold = (int)dims1[0];

	hsize_t dims1extend[1] = {n};
	hsize_t dims1new[1] = {(dims1[0]+dims1extend[0])};

	status = H5Dextend (dset, dims1new);

	/* Select a hyperslab in extened portion of dataset  */
	space = H5Dget_space (dset);
	status = H5Sselect_hyperslab(space, H5S_SELECT_SET, dims1, NULL, dims1extend, NULL);

	/* Copying results into hdf5 result type, and saving the references */
	resulth5_t *result5t = new resulth5_t[n];
	result5t = result5t - nold;
	for (int i=0; i<n; i++) {
	    memcpy(&result5t[nold+i], &result[i], HOFFSET(resulth5_t, res));
	    saveVarRef(&result5t[nold+i].res,  file, i+dims1[0], result[i].nq, result[i].k, varStr, "res",  result[i].res[0], "uint32");
	    saveVarRef(&result5t[nold+i].nres, file, i+dims1[0], result[i].nq, result[i].b+1, varStr, "nres", result[i].nres[0], "uint32");
	    if (result[i].stats != NULL)
		saveVarRef(&result5t[nold+i].stats, file, i+dims1[0], result[i].nq, STAT_DIM, varStr, "stats", result[i].stats[0], "double");
	    else
		result5t[nold+i].stats = 0;
	}
	status = H5Dwrite(dset, restype, space, space, H5P_DEFAULT, result5t);

	free(result5t+nold);

	// status = H5Sclose(space);
	status = H5Dclose(dset);
    }

    status = H5Sclose(space);
    status = H5Tclose(restype);
    status = H5Fclose(file);
}


