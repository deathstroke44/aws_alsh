#pragma once

#include <vector>
#include <cstdlib>
#include "../../util.h"
#include "../../register.h"
#include "../../def.h"

#include "../e2lsh.h"
#include "../s2alsh.h"
#include "../slalsh.h"
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include<bits/stdc++.h>
using namespace std;
using namespace std::chrono;


inline bool read_sigs(const std::string& filename, int N, int K, int L, std::vector<std::vector<uint64_t>>& sigs)
{
	FILE* fp = fopen(filename.c_str(), "rb");
	if(fp==NULL){
		return false;
	}

	int L_, N_;

	// fread(&K_, sizeof(int), 1, fp);
	fread(&L_, sizeof(int), 1, fp);
	fread(&N_, sizeof(int), 1, fp);

	// printf("K=%d, L=%d, N=%d\n", K_, L_, N);
	
	assert(L_>=L*K && N_==N);

	sigs.resize(N);
    for(int j=0;j<N;j++){
        sigs[j].resize(L);
    }

    //column major
    std::vector<uint64_t> tmp(N*K);
	for(int i=0;i<L;i++){
        fread(&tmp[0], sizeof(uint64_t), N*K, fp);

        for(int j=0;j<N;j++){
            uint64_t code = 0;
            for(int k=0;k<K;k++){
                int idx = k*N +j;

                code = hash_combine(code, tmp[idx]);
            }
            sigs[j][i] = code;
        }
	}

	fclose(fp);
	return true;
}

//hp: const float *-> [signatures] with the length max(Ls)*max(Ks)
template<class PF> 
double compute_signatures(
    int n, 
    int d, 
    const float** data, 
    const PF& hp,
    std::vector<int> Ks, 
    std::vector<int> Ls,
    const std::string& hashCacheName
)
{
    //compute hash table for all combination of K and L
    
    int maxK = Ks.back();
    int maxL = Ls.back();
    int nHashers = maxK*maxL;

    // int nsigs = maxK*maxL;

    std::vector<std::vector<uint64_t> > datasigs;
    datasigs.reserve(n); 
    using UNIQUE_FP = std::unique_ptr<FILE, decltype(&fclose)>;

    std::vector<UNIQUE_FP> fps;
    UNIQUE_FP fp(fopen(hashCacheName.c_str(), "wb"), &fclose);
    clock_t t; 
    t = clock();
    for(int i=0;i<n;i++){
        const float* datai = data[i];
        std::vector<uint64_t> datasig = hp(datai);
        // datasigs should have the length maxK*maxL

        if(i%10000==0){
            // printf("%d-th signatures\n", i);
        }

        assert(datasig.size()==nHashers);

        datasigs.emplace_back(datasig);
    }
    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 


    fwrite(&nHashers, sizeof(int), 1, fp.get());
    fwrite(&n, sizeof(int), 1, fp.get());
    std::vector<uint64_t> tmp(n);
    for(int i=0;i<nHashers;i++){
        for(int j=0;j<n;j++){
            tmp[j] = datasigs[j][i];
        }

        fwrite(&tmp[0], sizeof(uint64_t), n, fp.get());
    }
    return time_taken;
}



template<class QF> 
double compute_query_signatures(
    int qn, 
    int d, 
    const float** query,
    const float** weight,  
    const QF& hq,
    std::vector<int> Ks, 
    std::vector<int> Ls,
    const std::string& hashCacheName
)
{
    int maxK = Ks.back();
    int maxL = Ls.back();
    int nHashers = maxK*maxL;

    // int nsigs = maxK*maxL;

    std::vector<std::vector<uint64_t> > querysigs;
    querysigs.reserve(qn); 
    using UNIQUE_FP = std::unique_ptr<FILE, decltype(&fclose)>;

    std::vector<UNIQUE_FP> fps;
    UNIQUE_FP fp(fopen(hashCacheName.c_str(), "wb"), &fclose);
    auto start = high_resolution_clock::now();
    clock_t t; 
    t = clock();

    for(int i=0;i<qn;i++){
        const float* queryi = query[i];
        const float* weighti = weight[i];

        auto querysig = hq(queryi, weighti);
        // datasigs should have the length maxK*maxL

        if(i%10000==0){
            printf("%d-th signatures\n", i);
        }

        assert(querysigs.size()==nHashers);

        querysigs.emplace_back(querysig);   
    }
    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 


    fwrite(&nHashers, sizeof(int), 1, fp.get());
    fwrite(&qn, sizeof(int), 1, fp.get());
    std::vector<uint64_t> tmp(qn);
    for(int i=0;i<nHashers;i++){
        for(int j=0;j<qn;j++){
            tmp[j] = querysigs[j][i];
        }

        fwrite(&tmp[0], sizeof(uint64_t), qn, fp.get());
    }
    return time_taken;
}


//input:
//n: #points
//qn: #queries
//d: dimension
//data: pointer to data
//query: pointer to query
//weight: pointer to weight
//results: pointer to results
//hp: hash function for data    PF:: const Scalar* -> vector<uint64_t>
//hq: hash function for query   QF:: const Scalar* query -> const Scalar* weight -> vector<uint64_t> 
//fp, FILE*, 
//dataHashName: cache to store the signatures of data
//queryHashName: cache to store the signatures of query
template<class PF, class QF>
void fraction_recall( 
    int n,
    int qn, 
    int d, 
    std::vector<int>& Ks, 
    std::vector<int>& Ls,
    const float** data, 
    const float** query,
    const float** weight, 
    const Result** results,
    const PF& hp, 
    const QF& hq, 
    FILE* fp, 
    const std::string& dataHashName,
    const std::string& queryHashName
    )
{
    //try read sigs, if not success, compute from data

    std::vector<std::vector<uint64_t> > dataSigs;
    std::vector<std::vector<uint64_t> > querySigs;
    double indexing_time = 0;
    

    // fprintf(fp, "--%d, %d, %d--\n", n, qn, d);
    // fprintf(fp, "--%f, %f, %f--\n", data[0][0], data[0][1], data[0][2]);
    // fprintf(fp, "--%f, %f, %f--\n", data[1][0], data[1][1], data[1][2]);

    for(int K:Ks){
        for(int L:Ls){
            // printf("test K=%d, L=%d\n", K, L);
            // fprintf(fp, "K=%d, L=%d\n", K, L);

            if(!read_sigs(dataHashName, n, K, L, dataSigs)){
                indexing_time += compute_signatures(n, d, data, hp, Ks, Ls, dataHashName);
                bool isReadSuccessful = read_sigs(dataHashName, n, K, L, dataSigs);
                assert(isReadSuccessful);
            }

            //build an inverted list
            std::vector<std::unordered_map<uint64_t, std::vector<int> > > invLists(L);
            for(int i=0;i<n;i++){
                for(int j=0;j<L;j++){
                    uint64_t sigij = dataSigs[i][j];

                    invLists[j][sigij].push_back(i);
                    // fprintf(fp, "data=%d, L=%d sigij=%d\n", i, j, sigij);
                }
            }



            if(!read_sigs(queryHashName, qn, K, L, querySigs)){
                indexing_time += compute_query_signatures(qn, d, query, weight, hq, Ks, Ls, queryHashName);
                bool isReadSuccessful = read_sigs(queryHashName, qn, K, L, querySigs);
                assert(isReadSuccessful);
            }

            double step = 0.005;
            int stepn = step*n;
            // printf("test1 K=%d, L=%d\n", K, L);
            float recall = 0.0f;
            float map = 0.0f;
            clock_t t; 
            t = clock();
            for(int i=0;i<qn;i++){
                //compute totCnts
                auto& qsigi = querySigs[i];
                // printf("test2 K=%d, L=%d, qn=%d\n", K, L, i);
                std::unordered_set<int> checked;
                MinK_List que(MAXK);
                // fprintf(fp, "00, 00\n");

                for(int j=0;j<L;j++){
                    // printf("%d, %d\n", j, checked.size());
                    uint64_t qsigij  = qsigi[j];
                    // fprintf(fp, "Q=%d, Signature=%d\n", i, qsigi); 
                    //for each hashTable, check
                    auto it = invLists[j].find(qsigij);
                    if(it!=invLists[j].end()){
                        auto& postingList = it->second;
                        for(int dataid:postingList){
                            if(checked.find(dataid)==checked.end()){
                                //not checked!
                                float dw = calc_weighted_dist2(d, weight[i], query[i], data[dataid]);
                                checked.insert(dataid);
                                que.insert(dw, dataid);

                                // if(checked.size()%stepn==1){
                                //     float recall = calc_recall(10, results[i], &que);
                                //     printf("test2 recall=%f\n", recall);
                                //     // fprintf(fp, "Recall %f, %f\n", checked.size()*1./n, recall);
                                // }
                            }
                        }
                    }
                }
                // fprintf(fp, "K=%d, L=%d, Q=%d, Candidate=%d\n", K, L, i, checked.size());
                // int ln = que.size()-1;
                // fprintf(fp, "%dth  query: ", i);
                // while (ln >= 0) {
                //     fprintf(fp, "%d %f || ", que.ith_id(ln), que.ith_key(ln));
                //     ln--;
                // }
                // fprintf(fp, "\n");
                // ln = que.size()-1;
                // while (ln >= 0) {
                //     fprintf(fp, "%d %f || ", results[i][ln].id_, results[i][ln].key_);
                //     ln--;
                // }
                // fprintf(fp, "\n");
                recall += calc_recall(MAXK, results[i], &que);
                map += calc_map(MAXK, results[i], &que);
                
            }
            t = clock() - t; 
            double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
            fprintf(fp, "L=%d K=%d MAXK:%d recall: %f map: %f indexing time: %f search time: %f\n", L, K, MAXK, recall/qn, map/qn, indexing_time, time_taken);
        }
    }
}
