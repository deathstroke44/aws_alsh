#include "util.h"
#include<bits/stdc++.h>
using namespace std;
using namespace std::chrono;




int readFVecsFromExternal(
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data) {
  FILE *infile = fopen(fname, "rb");
  if (infile == NULL) {
    std::cout << "File not found" << std::endl;
    return 1;
  }
  
  int rowCt = 0;
  int dimen;
  while (true) {
    if (fread(&dimen, sizeof(int), 1, infile) == 0) {
      break;
    }
    if (dimen != d) {
      std::cout << "N and actual dimension mismatch" << std::endl;
      return 1;
    }
    std::vector<float> v(d);
    if(fread(v.data(), sizeof(float), dimen, infile) == 0) {
      std::cout << "Error when reading" << std::endl;
    };
    
    for (int i=0; i<d; i++) {
      data[rowCt][i] = v[i];
    }

    rowCt++;
    
    if (n != -1 && rowCt >= n) {
      break;
    }
  }
  // std::cout<<"Row count test: "<<rowCt<<std::endl;

  if (fclose(infile)) {
    std::cout << "Could not close data file" << std::endl;
	return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
int read_data(						// read data/query set from disk
	int   n,							// number of data/query objects
	int   d,			 				// dimensionality
	const char *fname,					// address of data/query set
	float **data)						// data/query objects (return)
{
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int i   = 0;
	int tmp = -1;
	while (!feof(fp) && i < n) {
		fscanf(fp, "%d", &tmp);
		for (int j = 0; j < d; ++j) {
			fscanf(fp, " %f", &data[i][j]);
		}
//		fscanf(fp, "\n");

		++i;
	}
//	assert(feof(fp) && i == n);
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R)							// ground truth results (return)
{
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Could not open %s\n", fname);
		return 1;
	}

	int tmp1 = -1;
	int tmp2 = -1;
	fscanf(fp, "%d %d\n", &tmp1, &tmp2);
//	assert(tmp1 == qn && tmp2 == MAXK);
	assert(tmp2 == MAXK);

	for (int i = 0; i < qn; ++i) {
		for (int j = 0; j < MAXK; ++j) {
			fscanf(fp, "%d %f ", &R[i][j].id_, &R[i][j].key_);
		}
		fscanf(fp, "\n");
	}
	fclose(fp);

	return 0;
}

// -----------------------------------------------------------------------------
int read_ground_truthV2(				// read ground truth results from disk
	int qn,
	int d,								// number of query objects
	const char *fname,					// address of truth set
	Result **R,
	float **data,
	float **query)							// ground truth results (return)
{
	FILE *infile = fopen(fname, "rb");
	if (infile == NULL) {
		std::cout << "File not found" << std::endl;
		return 1;
	}
	
	int rowCt = 0;
	int dimen;
	while (true) {
		if (fread(&dimen, sizeof(int), 1, infile) == 0) {
			break;
		}
		if (dimen != MAXK) {
			std::cout << "N and actual dimension mismatch" << std::endl;
			return 1;
		}
		std::vector<int> v(MAXK);
		if(fread(v.data(), sizeof(int), dimen, infile) == 0) {
			std::cout << "Error when reading" << std::endl;
		};
		
		for (int i=0; i<MAXK; i++) {
			R[rowCt][i].id_ = v[i]+1;
			R[rowCt][i].key_ = calc_l2_dist(d, query[rowCt], data[v[i]]);
		}

		rowCt++;
		
		if (qn != -1 && rowCt >= qn) {
			break;
		}
	}
	// std::cout<<"Row count test: "<<rowCt<<std::endl;

	if (fclose(infile)) {
		std::cout << "Could not close data file" << std::endl;
		return 1;
	}
	return 0;
}

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	double ret = 0.0f;
	for (int i = 0; i < dim; ++i) {
		ret += (p1[i] * p2[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_angle(				// calc angle
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	return acos(calc_cosangle(dim, p1, p2));
}

// -----------------------------------------------------------------------------
float calc_cosangle(				// calc cos(angle)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	double ret = 0.0f;
	double norm0 = 0., norm1 = 0.;
	for (int i = 0; i < dim; ++i) {
		ret += (p1[i] * p2[i]);
		norm0 += p1[i] * p1[i];
		norm1 += p2[i] * p2[i];
	}
	if(norm0==0 || norm1==0){
		return 0;
	}
	return ret/sqrt(norm0*norm1);
}

// -----------------------------------------------------------------------------
float calc_l2_sqr(					// calc L2 square distance
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float diff = 0.0f;
	float ret  = 0.0f;
	for (int i = 0; i < dim; ++i) {
		diff = p1[i] - p2[i];
		ret += diff * diff;
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L2 distance
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	return sqrt(calc_l2_sqr(dim, p1, p2));
}

// -----------------------------------------------------------------------------
float calc_l1_dist(					// calc L1 distance
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	float ret = 0.0f;
	for (int i = 0; i < dim; ++i) {
		ret += fabs(p1[i] - p2[i]);
	}
	return ret;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list)					// results returned by algorithms
{
	int i = list->size()-1;
	int last = k - 1;
	//loop until list->ithkey >= R[last].key
	while (i >= 0 && R[last].key_ - list->ith_key(i) > EPS) {
		i--;
	}
	return (i + 1) * 100.0f / k;
}

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results
	MinK_List *list)					// results returned by algorithms
{
	int cnt=0;
	for (int p=0; p<list->size();p++){
		for (int r=0;r<k;r++) {
			if(list->ith_id(p)+1==R[r].id_) {
				cnt++;
				break;
			}
		}
	}
	return cnt * 1.0f / k;
}

// -----------------------------------------------------------------------------
float calc_map(					// calc map (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results
	MinK_List *list)					// results returned by algorithms
{
	int cnt=0;
	int _k =k;
	k=min(list->size(),k);
	int ap=0;
	for (int p=1; p<=k;p++){
		bool isR_kExact = false;
		for (int r=0;r<k;r++) {
			if(list->ith_id(p)+1==R[r].id_) {
				isR_kExact=true;
			}
		}
		if (isR_kExact) {
			int ct = 0;
			for (int j=0; j<p; j++) {
				for (int jj=0; jj<p; jj++) {
					if (list->ith_id(j)+1==R[jj].id_) {
						ct++;
						break;
					}
				}
			}
			ap += (double)ct/p;
		}
	}
	return ap * 1.0f / _k;
}

// -----------------------------------------------------------------------------
int get_hits(						// get the number of hits between two ID list
	int   k,							// top-k value
	int   t,							// top-t value
	const Result *R,					// ground truth results 
	MaxK_List *list)					// results returned by algorithms
{
	int i = k - 1;
	int last = t - 1;
	while (i >= 0 && R[last].key_ - list->ith_key(i) > EPS) {
		i--;
	}
	return min(t, i + 1);
}

float calc_weighted_dist2(			// calc inner product
	int   dim,							// dimension
	const float *w,
	const float *p1,					// 1st point
	const float *p2)					// 2nd point
{
	double ret = 0.0f;
	for (int i = 0; i < dim; ++i) {
		ret += w[i]*(p1[i]-p2[i])*(p1[i]-p2[i]);
	}
	return ret;
}


int calc_hamming_dist(			// calc inner product
	int   dim,		
	const uint8_t *p1,					// 1st point
	const uint8_t *p2)					// 2nd point
{
	int tail = dim%8;
	int ret = 0;
	for(int i=0;i<tail;i++){
		ret += get_num_bits8(p1[i]^p2[i]);
	}
	ret += calc_hamming_dist(dim/8, (const uint64_t*)(p1+tail), (const uint64_t*)(p2+tail));
	return ret;
}

int calc_hamming_dist(			// calc inner product
	int   dim,		
	const uint64_t *p1,					// 1st point
	const uint64_t *p2)					// 2nd point
{
	int ret = 0;
	for(int i=0;i<dim;i++){
		ret += get_num_bits64(p1[i]^p2[i]);
	}
	return ret;
}

float calc_ratio(
	int k, 
	const Result *Rs, 
	MinK_List *list)
{
	if(list->size()<k){
		return sqrt(1e9);
	}
	double ret  = (list->ith_key(k-1)+1e-9) / (Rs[k-1].key_+1e-9);
	if(ret<0){
		ret = 1e9;
	} else if(ret<1){
		ret = 1/ret;
	}
	return sqrt(ret);
}
float calc_ratio(
	int k, 
	const Result *Rs, 
	MaxK_List *list)
{
	if(list->size()<k){
		return sqrt(1e9);
	}
	double ret  = (list->ith_key(k-1)+1e-9) / (Rs[k-1].key_+1e-9);
	if(ret<0){
		ret = 1e9;
	} else if(ret<1){
		ret = 1/ret;
	}
	return sqrt(ret);
}


void calc_min_max(
	int n, 
	int qn, 
	int d, 
	const float** data, 
	const float** query, 
	float& maxx, 
	float& minx)
{
	maxx=0;
	minx=1e10;

	for(int i=0;i<n;i++){
		for(int j=0;j<d;j++){
			maxx = std::max(data[i][j], maxx);
			minx = std::min(data[i][j], minx);
		}
	}
	for(int i=0;i<qn;i++){
		for(int j=0;j<d;j++){
			maxx = std::max(query[i][j], maxx);
			minx = std::min(query[i][j], minx);
		}
	}
}