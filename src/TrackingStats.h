//
// Created by Ales Varabyou on 3/4/20.
//

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <iomanip>

#ifndef GTEX_STATS_TRACKINGSTATS_H
#define GTEX_STATS_TRACKINGSTATS_H

#define FIXED_FLOAT(x) std::fixed <<std::setprecision(2)<<(x) // https://stackoverflow.com/questions/5907031/printing-the-correct-number-of-decimal-points-with-cout

class TrackingStats{
public:
    TrackingStats() = default;
    ~TrackingStats() = default;

    // save_all_tx_stats
    std::vector<int> num_tx_per_all_loc,num_tx_per_all_loc_real,num_tx_per_all_loc_intronic,num_tx_per_all_loc_splicing,num_tx_per_all_loc_intergenic,num_tx_per_all_loc_left;
    std::vector<float> sum_tx_per_all_loc,sum_tx_per_all_loc_real,sum_tx_per_all_loc_intronic,sum_tx_per_all_loc_splicing,sum_tx_per_all_loc_intergenic,sum_tx_per_all_loc_left;
    // save_tissue_tx_stats
    std::vector<int> num_tx_per_tissue_loc,num_tx_per_tissue_loc_real,num_tx_per_tissue_loc_intronic,num_tx_per_tissue_loc_splicing,num_tx_per_tissue_loc_intergenic,num_tx_per_tissue_loc_left;
    std::vector<float> sum_tx_per_tissue_loc,sum_tx_per_tissue_loc_real,sum_tx_per_tissue_loc_intronic,sum_tx_per_tissue_loc_splicing,sum_tx_per_tissue_loc_intergenic,sum_tx_per_tissue_loc_left;
    // save_sample_tx_stats
    std::vector<int> num_tx_per_sample_loc,num_tx_per_sample_loc_real,num_tx_per_sample_loc_intronic,num_tx_per_sample_loc_splicing,num_tx_per_sample_loc_intergenic,num_tx_per_sample_loc_left;
    std::vector<float> sum_tx_per_sample_loc,sum_tx_per_sample_loc_real,sum_tx_per_sample_loc_intronic,sum_tx_per_sample_loc_splicing,sum_tx_per_sample_loc_intergenic,sum_tx_per_sample_loc_left;
    // save_tissue_tx_stats2
    std::vector<int> num_tx_per_tissue_loc2,num_tx_per_tissue_loc_real2,num_tx_per_tissue_loc_intronic2,num_tx_per_tissue_loc_splicing2,num_tx_per_tissue_loc_intergenic2,num_tx_per_tissue_loc_left2;
    // save_sample_tx_stats2
    std::vector<int> num_tx_per_sample_loc2,num_tx_per_sample_loc_real2,num_tx_per_sample_loc_intronic2,num_tx_per_sample_loc_splicing2,num_tx_per_sample_loc_intergenic2,num_tx_per_sample_loc_left2;
    // save_cov_fracs
    std::vector<std::tuple<int,float,int,float,int,float>> cov_fracs; // number of real transcripts, sum of tpms of real, number of splicing transcripts, sum of tpms of splicing, num of intronic transcirpts, sum of intronic
    // save_cov_joined
    std::vector<std::tuple<std::vector<float>,std::vector<float>,std::vector<float>>> covs_joined;
    // save_tpm_fracs
    std::vector<std::tuple<int,float,int,float,int,float>> tpm_fracs; // number of real transcripts, sum of tpms of real, number of splicing transcripts, sum of tpms of splicing, num intronic, sum intronic
    // save_tpm_joined
    std::vector<std::tuple<std::vector<float>,std::vector<float>,std::vector<float>>> tpms_joined;
    // save_sample_cov
    std::vector<float> cov_sample_real,cov_sample_intronic,cov_sample_splicing,cov_sample_intergenic,cov_sample_left;
    std::vector<float> fpkm_sample_real,fpkm_sample_intronic,fpkm_sample_splicing,fpkm_sample_intergenic,fpkm_sample_left;
    std::vector<float> tpm_sample_real,tpm_sample_intronic,tpm_sample_splicing,tpm_sample_intergenic,tpm_sample_left;
    // save_sample_loc
    std::vector<std::array<int,9>> num_locs_sample; // real,splice,int,all_real,pol
    // save_tissue_loc
    std::vector<std::tuple<std::string,int,int>> num_locs_tissue;
    // save_sample_txs
    std::map<std::string,std::tuple<int,int,int,int>> sample_txs;
    // save_tissue_tx_stats3
    std::vector<std::tuple<int,int,int,int>> num_tx_per_tissue_loc3;
    // save_sample_tx_stats3
    std::vector<std::tuple<int,int,int,int>> num_tx_per_sample_loc3;
    // save_sample_tx_stats4
    std::map<std::string,std::tuple<int,int,int,int,std::string>> num_tx_per_sample_loc4;
    std::map<std::string,std::array<int,4>> sample2tissue_loc_txs;
    std::pair<std::map<std::string,std::array<int,4>>::iterator,bool> sample2tissue_loc_txs_it;
    // save_tissue_tx_stats5
    std::map<std::string,std::tuple<int,int,int,int>> num_tx_per_tissue_loc5;
    // save_sample_tpm_gauss
    std::vector<std::vector<std::pair<float,float>>> gauss_sample_tpm_per_tissue_loc;
    // save_sample_tx_gauss
    std::vector<std::tuple<float,float,int>> gauss_sample_tx_per_tissue_loc;
    // save_sample_gauss
    std::vector<std::pair<std::tuple<float,float,int,int>,std::vector<std::tuple<float,float,std::vector<float>>>>> gauss_sample_per_tissue_loc;
    // save_sample_tx_stats6
    std::map<std::pair<std::string,std::string>, // tissue and locus
            std::tuple<std::array<int,4>, // 4 integers contain total number of transcripts in a locus
                    std::array<std::set<std::string>,4>, // 4 integers contain number of transcripts in a locus of this tissue specifically
                    std::map<std::string, // sample name
                            std::array<std::vector<float>,4>// tpms of each type
                    >
            >
    > num_tx_per_sample_loc6;
    // save_sample_tx_stats5
    // The following container structure is:
    //  key - sample_loc_id
    //  value tuple:
    //     <0> - int - number of real transcripts
    //     <1> - int - number of splicing transcripts
    //     <2> - int - number of intronic transcripts
    //     <3> - int - number of intergenic transcripts
    //     <4> - float - tpms observed in real transcripts
    //     <5> - float - tpms observed in splicing transcripts
    //     <6> - float - tpms observed in intronic transcripts
    //     <7> - float - tpms observed in intergenic transcripts
    std::map<std::pair<std::string,std::string>,std::tuple<std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>>> num_tx_per_sample_loc5;
    // save_loc_stats
    int num_avg_real_locs_per_tissue=0,num_avg_noise_locs_per_tissue=0,num_avg_undef_locs_per_tissue=0;
    // save_num_tx_tissue
    std::map<std::string,std::tuple<int,int,int,int>> num_tx_per_tissue;

    std::vector<float> tissue_tx_freq_real,tissue_tx_freq_intronic,tissue_tx_freq_splicing,tissue_tx_freq_intergenic;

    void save_stats(std::string base_out_fname){
        std::cout<<"saving stats"<<std::endl;
//        save_all_tx_stats(base_out_fname);
//        save_tissue_tx_stats(base_out_fname);
//        save_sample_tx_stats(base_out_fname);
//        save_sample_cov(base_out_fname);
//        save_loc_stats(base_out_fname);
//        save_tissue_tx_stats2(base_out_fname);
//        save_sample_tx_stats2(base_out_fname);
//        save_cov_fracs(base_out_fname);
//        save_cov_joined(base_out_fname);
//        save_tpm_fracs(base_out_fname);
//        save_tpm_joined(base_out_fname);
//        save_tissue_tx_stats3(base_out_fname);
//        save_sample_tx_stats3(base_out_fname);
//        save_tissue_tx_stats5(base_out_fname);
//        save_num_tx_tissue(base_out_fname);
//        save_sample_tx_gauss(base_out_fname);
//        save_sample_tpm_gauss(base_out_fname);
//        save_sample_gauss(base_out_fname);
        save_tissue_loc(base_out_fname);
        save_sample_loc(base_out_fname);
        save_sample_txs(base_out_fname);
//        save_sample_tx_stats4(base_out_fname);
        save_sample_tx_stats5(base_out_fname);
        save_sample_tx_stats6(base_out_fname);
        save_tx_freq(base_out_fname);
        std::cout<<"done saving stats"<<std::endl;
    }

private:
    void save_all_tx_stats(std::string base_out_fname);

    void save_tissue_tx_stats(std::string base_out_fname);

    void save_sample_tx_stats(std::string base_out_fname);

    void save_tissue_tx_stats2(std::string base_out_fname);

    void save_sample_tx_stats2(std::string base_out_fname);

    // for real loci only - used to compute fraction of real to splicing transcripts per locus
    void save_cov_fracs(std::string base_out_fname);

    void save_cov_joined(std::string base_out_fname);

    void save_tpm_fracs(std::string base_out_fname);

    void save_tpm_joined(std::string base_out_fname);

    void save_sample_cov(std::string base_out_fname);

    void save_sample_loc(std::string base_out_fname);

    // we need to have distribution of loci per tissue and per sample as well to choose from
    void save_tissue_loc(std::string base_out_fname);

    void save_sample_txs(std::string base_out_fname);

    void save_tissue_tx_stats3(std::string base_out_fname);

    void save_sample_tx_stats3(std::string base_out_fname);

    void save_sample_tx_stats4(std::string base_out_fname);

    void save_tissue_tx_stats5(std::string base_out_fname);

    void save_sample_tpm_gauss(std::string base_out_fname);

    void save_sample_tx_gauss(std::string base_out_fname);

    void save_sample_gauss(std::string base_out_fname);

    void save_sample_tx_stats6(std::string base_out_fname);

    void save_sample_tx_stats5(std::string base_out_fname);

    void save_loc_stats(std::string base_out_fname);

    void save_num_tx_tissue(std::string base_out_fname);

    void save_tx_freq(std::string base_out_fname);
};


#endif //GTEX_STATS_TRACKINGSTATS_H
