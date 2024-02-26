#include "e2lsh.h"


E2LSH::E2LSH(int d, int K, int L, Scalar r):
    dim(d), K(K), L(L), r(r)
{
    assert(K>0 && L>0);

    //generate random projection vectors    
	std::default_random_engine rng;
	std::normal_distribution<double> normal(0.,1.0);
	std::uniform_real_distribution<double> uniform(0.,r);
    
    a.resize(K*L);
    b.resize(K*L);

	for(int i=0;i<K*L;i++){
		a[i].resize(d);
		for(int j=0;j<d;j++){
			a[i][j] = normal(rng);
		}
		b[i] = uniform(rng);
	}
}

E2LSH::~E2LSH()
{

}

std::vector<uint64_t> E2LSH::hash_data(const Scalar* data)
{
    // FILE* fp;
    // bool flag=true;
    // Creates a file "demo_file"
    // with file access as write-plus mode
    // std::cout<<"Reach hrer "<<flag<<std::endl;
    // fp = fopen("projection.txt", "a+");
    std::vector<uint64_t> ret;
    // fprintf(fp, "--%f, %f, %f--\n", data[0][0], data[0][1], data[0][2]);
    // fprintf(fp, "Data: ");
    // for(int i=0;i<dim;i++) {
    //     if(flag) fprintf(fp, "%f ", data[i]);
    // }
    // if(flag) fprintf(fp, "\n"); 
    for(int l=0;l<L;l++){
        uint64_t sigl = 0;
        // if(flag) fprintf(fp, "%d\n", l);
        for(int k=0;k<K;k++){
            double projection = 0.;
            for(int i=0;i<dim;i++){
                double x = data[i];
                projection += x*a[l*K+k][i];
                // if(flag) fprintf(fp, "Projection temp --%f, %f, %f, %f --\n", x, a[l*K+k][i], x*a[l*K+k][i], projection);
            }
            projection += b[l];
            int projection_bucket = projection/r;

            sigl = hash_combine(sigl, uint64_t(projection_bucket));
            // if(flag) fprintf(fp, "projection final --%f, %d, %d --\n", projection, projection_bucket, sigl);
        }
        ret.push_back(sigl);
    }
    // fclose(fp);
    return ret;
}