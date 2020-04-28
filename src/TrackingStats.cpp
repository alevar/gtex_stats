//
// Created by sparrow on 3/4/20.
//

#include "TrackingStats.h"

void TrackingStats::save_all_tx_stats(std::string base_out_fname){
    if(num_tx_per_all_loc.size()!=num_tx_per_all_loc_real.size()) {
        std::cerr<<"all loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_all_loc: "<<num_tx_per_all_loc.size()<<std::endl;
        std::cerr<<"num_tx_per_all_loc_real: "<<num_tx_per_all_loc_real.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_all_loc_real.size()!=num_tx_per_all_loc_intronic.size()) {
        std::cerr<<"all loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_all_loc_real: "<<num_tx_per_all_loc_real.size()<<std::endl;
        std::cerr<<"num_tx_per_all_loc_intronic: "<<num_tx_per_all_loc_intronic.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_all_loc_intronic.size() != num_tx_per_all_loc_splicing.size()) {
        std::cerr<<"all loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_all_loc_intronic: "<<num_tx_per_all_loc_intronic.size()<<std::endl;
        std::cerr << "num_tx_per_all_loc_splicing: " << num_tx_per_all_loc_splicing.size() << std::endl;
        exit(-1);
    }
    else if(num_tx_per_all_loc_splicing.size() != num_tx_per_all_loc_intergenic.size()) {
        std::cerr<<"all loc vectors not equal size"<<std::endl;
        std::cerr << "num_tx_per_all_loc_splicing: " << num_tx_per_all_loc_splicing.size() << std::endl;
        std::cerr<<"num_tx_per_all_loc_intergenic: "<<num_tx_per_all_loc_intergenic.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_all_loc_intergenic.size()!=num_tx_per_all_loc_left.size()) {
        std::cerr<<"all loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_all_loc_intergenic: "<<num_tx_per_all_loc_intergenic.size()<<std::endl;
        std::cerr<<"num_tx_per_all_loc_left: "<<num_tx_per_all_loc_left.size()<<std::endl;
        exit(-1);
    }
    else{
        std::string num_tx_per_all_loc_fname = base_out_fname+".num_tx_per_all_loc";
        std::ofstream num_tx_per_all_loc_ss(num_tx_per_all_loc_fname.c_str());
        num_tx_per_all_loc_ss<<"all,real,splicing,intronic,intergenic,left"<<std::endl;

        for(int i=0;i<this->num_tx_per_all_loc.size();i++){
            num_tx_per_all_loc_ss << num_tx_per_all_loc[i] << "," << num_tx_per_all_loc_real[i] << "," << num_tx_per_all_loc_intronic[i] << "," << num_tx_per_all_loc_splicing[i] << "," << num_tx_per_all_loc_intergenic[i] << "," << num_tx_per_all_loc_left[i] << std::endl;
        }

        num_tx_per_all_loc_ss.close();

        std::string sum_tx_per_all_loc_fname = base_out_fname+".sum_tx_per_all_loc";
        std::ofstream sum_tx_per_all_loc_ss(sum_tx_per_all_loc_fname.c_str());
        sum_tx_per_all_loc_ss<<"all,real,splicing,intronic,intergenic,left"<<std::endl;

        for(int i=0;i<this->sum_tx_per_all_loc.size();i++){
            sum_tx_per_all_loc_ss << sum_tx_per_all_loc[i] << "," << sum_tx_per_all_loc_real[i] << "," << sum_tx_per_all_loc_intronic[i] << "," << sum_tx_per_all_loc_splicing[i] << "," << sum_tx_per_all_loc_intergenic[i] << "," << sum_tx_per_all_loc_left[i] << std::endl;
        }

        sum_tx_per_all_loc_ss.close();
    }
}

void TrackingStats::save_tissue_tx_stats(std::string base_out_fname){
    if(num_tx_per_tissue_loc.size()!=num_tx_per_tissue_loc_real.size()) {
        std::cerr<<"tissue loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc: "<<num_tx_per_tissue_loc.size()<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc_real: "<<num_tx_per_tissue_loc_real.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_tissue_loc_real.size()!=num_tx_per_tissue_loc_intronic.size()) {
        std::cerr<<"tissue loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc_real: "<<num_tx_per_tissue_loc_real.size()<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc_intronic: "<<num_tx_per_tissue_loc_intronic.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_tissue_loc_intronic.size() != num_tx_per_tissue_loc_splicing.size()) {
        std::cerr<<"tissue loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc_intronic: "<<num_tx_per_tissue_loc_intronic.size()<<std::endl;
        std::cerr << "num_tx_per_tissue_loc_splicing: " << num_tx_per_tissue_loc_splicing.size() << std::endl;
        exit(-1);
    }
    else if(num_tx_per_tissue_loc_splicing.size() != num_tx_per_tissue_loc_intergenic.size()) {
        std::cerr<<"tissue loc vectors not equal size"<<std::endl;
        std::cerr << "num_tx_per_tissue_loc_splicing: " << num_tx_per_tissue_loc_splicing.size() << std::endl;
        std::cerr<<"num_tx_per_tissue_loc_intergenic: "<<num_tx_per_tissue_loc_intergenic.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_tissue_loc_intergenic.size()!=num_tx_per_tissue_loc_left.size()) {
        std::cerr<<"tissue loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc_intergenic: "<<num_tx_per_tissue_loc_intergenic.size()<<std::endl;
        std::cerr<<"num_tx_per_tissue_loc_left: "<<num_tx_per_tissue_loc_left.size()<<std::endl;
        exit(-1);
    }
    else{
        std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc";
        std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
        num_tx_per_tissue_loc_ss<<"all,real,splicing,intronic,intergenic,left"<<std::endl;

        for(int i=0;i<this->num_tx_per_tissue_loc.size();i++){
            num_tx_per_tissue_loc_ss << num_tx_per_tissue_loc[i] << "," << num_tx_per_tissue_loc_real[i] << "," << num_tx_per_tissue_loc_intronic[i] << "," << num_tx_per_tissue_loc_splicing[i] << "," << num_tx_per_tissue_loc_intergenic[i] << "," << num_tx_per_tissue_loc_left[i] << std::endl;
        }

        num_tx_per_tissue_loc_ss.close();

        std::string sum_tx_per_tissue_loc_fname = base_out_fname+".sum_tx_per_tissue_loc";
        std::ofstream sum_tx_per_tissue_loc_ss(sum_tx_per_tissue_loc_fname.c_str());
        sum_tx_per_tissue_loc_ss<<"all,real,splicing,intronic,intergenic,left"<<std::endl;

        for(int i=0;i<this->sum_tx_per_tissue_loc.size();i++){
            sum_tx_per_tissue_loc_ss << sum_tx_per_tissue_loc[i] << "," << sum_tx_per_tissue_loc_real[i] << "," << sum_tx_per_tissue_loc_intronic[i] << "," << sum_tx_per_tissue_loc_splicing[i] << "," << sum_tx_per_tissue_loc_intergenic[i] << "," << sum_tx_per_tissue_loc_left[i] << std::endl;
        }

        sum_tx_per_tissue_loc_ss.close();
    }
}

void TrackingStats::save_sample_tx_stats(std::string base_out_fname){
    if(num_tx_per_sample_loc.size()!=num_tx_per_sample_loc_real.size()) {
        std::cerr<<"sample loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_sample_loc: "<<num_tx_per_sample_loc.size()<<std::endl;
        std::cerr<<"num_tx_per_sample_loc_real: "<<num_tx_per_sample_loc_real.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_sample_loc_real.size()!=num_tx_per_sample_loc_intronic.size()) {
        std::cerr<<"sample loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_sample_loc_real: "<<num_tx_per_sample_loc_real.size()<<std::endl;
        std::cerr<<"num_tx_per_sample_loc_intronic: "<<num_tx_per_sample_loc_intronic.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_sample_loc_intronic.size() != num_tx_per_sample_loc_splicing.size()) {
        std::cerr<<"sample loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_sample_loc_intronic: "<<num_tx_per_sample_loc_intronic.size()<<std::endl;
        std::cerr << "num_tx_per_sample_loc_splicing: " << num_tx_per_sample_loc_splicing.size() << std::endl;
        exit(-1);
    }
    else if(num_tx_per_sample_loc_splicing.size() != num_tx_per_sample_loc_intergenic.size()) {
        std::cerr<<"sample loc vectors not equal size"<<std::endl;
        std::cerr << "num_tx_per_sample_loc_splicing: " << num_tx_per_sample_loc_splicing.size() << std::endl;
        std::cerr<<"num_tx_per_sample_loc_intergenic: "<<num_tx_per_sample_loc_intergenic.size()<<std::endl;
        exit(-1);
    }
    else if(num_tx_per_sample_loc_intergenic.size()!=num_tx_per_sample_loc_left.size()) {
        std::cerr<<"sample loc vectors not equal size"<<std::endl;
        std::cerr<<"num_tx_per_sample_loc_intergenic: "<<num_tx_per_sample_loc_intergenic.size()<<std::endl;
        std::cerr<<"num_tx_per_sample_loc_left: "<<num_tx_per_sample_loc_left.size()<<std::endl;
        exit(-1);
    }
    else{
        std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc";
        std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
        num_tx_per_sample_loc_ss<<"all,real,splicing,intronic,intergenic,left"<<std::endl;

        for(int i=0;i<this->num_tx_per_sample_loc.size();i++){
            num_tx_per_sample_loc_ss << num_tx_per_sample_loc[i] << "," << num_tx_per_sample_loc_real[i] << "," << num_tx_per_sample_loc_intronic[i] << "," << num_tx_per_sample_loc_splicing[i] << "," << num_tx_per_sample_loc_intergenic[i] << "," << num_tx_per_sample_loc_left[i] << std::endl;
        }

        num_tx_per_sample_loc_ss.close();

        std::string sum_tx_per_sample_loc_fname = base_out_fname+".sum_tx_per_sample_loc";
        std::ofstream sum_tx_per_sample_loc_ss(sum_tx_per_sample_loc_fname.c_str());
        sum_tx_per_sample_loc_ss<<"all,real,splicing,intronic,intergenic,left"<<std::endl;

        for(int i=0;i<this->sum_tx_per_sample_loc.size();i++){
            sum_tx_per_sample_loc_ss << sum_tx_per_sample_loc[i] << "," << sum_tx_per_sample_loc_real[i] << "," << sum_tx_per_sample_loc_intronic[i] << "," << sum_tx_per_sample_loc_splicing[i] << "," << sum_tx_per_sample_loc_intergenic[i] << "," << sum_tx_per_sample_loc_left[i] << std::endl;
        }

        sum_tx_per_sample_loc_ss.close();
    }
}

void TrackingStats::save_tissue_tx_stats2(std::string base_out_fname){
    std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc2";
    std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
    num_tx_per_tissue_loc_ss<<"total,real,splicing,intronic,intergenic,left"<<std::endl;

    for(int i=0;i<this->num_tx_per_tissue_loc2.size();i++){
        num_tx_per_tissue_loc_ss<<num_tx_per_tissue_loc2[i]<<","<<num_tx_per_tissue_loc_real2[i]<<","<<num_tx_per_tissue_loc_intronic2[i]<<","<<num_tx_per_tissue_loc_splicing2[i]<<","<<num_tx_per_tissue_loc_intergenic2[i]<<","<<num_tx_per_tissue_loc_left2[i]<<std::endl;
    }

    num_tx_per_tissue_loc_ss.close();
}

void TrackingStats::save_sample_tx_stats2(std::string base_out_fname){
    std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc2";
    std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
    num_tx_per_sample_loc_ss<<"total,real,splicing,intronic,intergenic,left"<<std::endl;

    for(int i=0;i<this->num_tx_per_sample_loc2.size();i++){
        num_tx_per_sample_loc_ss<<num_tx_per_sample_loc2[i]<<","<<num_tx_per_sample_loc_real2[i]<<","<<num_tx_per_sample_loc_intronic2[i]<<","<<num_tx_per_sample_loc_splicing2[i]<<","<<num_tx_per_sample_loc_intergenic2[i]<<","<<num_tx_per_sample_loc_left2[i]<<std::endl;
    }

    num_tx_per_sample_loc_ss.close();
}

// for real loci only - used to compute fraction of real to splicing transcripts per locus
void TrackingStats::save_cov_fracs(std::string base_out_fname){
    std::string frac_fname = base_out_fname+".cov_fracs";
    std::ofstream frac_ss(frac_fname.c_str());
    frac_ss << "nt_real,sum_real,nt_splicing,sum_splicing,nt_intronic,sum_intronic" << std::endl;
    for(auto& t : cov_fracs){
        frac_ss << std::get<0>(t) << "," << std::get<1>(t) << "," << std::get<2>(t) << "," << std::get<3>(t) << "," << std::get<4>(t) << "," << std::get<5>(t) << std::endl;
    }
    frac_ss.close();
}

void TrackingStats::save_cov_joined(std::string base_out_fname){
    std::string joined_fname = base_out_fname+".covs_joined";
    std::ofstream joined_ss(joined_fname.c_str());
    joined_ss << "nt_real,nt_splicing,nt_intronic,real,splicing,intronic" << std::endl;
    for(auto& t : covs_joined){
        joined_ss << std::get<0>(t).size() << "," << std::get<1>(t).size() << "," << std::get<2>(t).size() << ",";
        if(std::get<0>(t).size() == 0){
            joined_ss << "0;";
        }
        for(auto& s : std::get<0>(t)){
            joined_ss << s <<";";
        }
        joined_ss.seekp(-1, std::ios_base::end);
        joined_ss << ",";
        if(std::get<1>(t).size() == 0){
            joined_ss << "0;";
        }
        for(auto& s : std::get<1>(t)){
            joined_ss << s <<";";
        }
        joined_ss.seekp(-1, std::ios_base::end);
        joined_ss << ",";
        if(std::get<2>(t).size() == 0){
            joined_ss << "0;";
        }
        for(auto& s : std::get<2>(t)){
            joined_ss << s <<";";
        }
        joined_ss.seekp(-1, std::ios_base::end);
        joined_ss << std::endl;
    }
    joined_ss.close();
}

void TrackingStats::save_tpm_fracs(std::string base_out_fname){
    std::string frac_fname = base_out_fname+".tpm_fracs";
    std::ofstream frac_ss(frac_fname.c_str());
    frac_ss << "nt_real,sum_real,nt_splicing,sum_splicing,nt_intronic,sum_intronic" << std::endl;
    for(auto& t : tpm_fracs){
        frac_ss << std::get<0>(t) << "," << std::get<1>(t) << "," << std::get<2>(t) << "," << std::get<3>(t) << "," << std::get<4>(t) << "," << std::get<5>(t) << std::endl;
    }
    frac_ss.close();
}

void TrackingStats::save_tpm_joined(std::string base_out_fname){
    std::string joined_fname = base_out_fname+".tpms_joined";
    std::ofstream joined_ss(joined_fname.c_str());
    joined_ss << "nt_real,nt_splicing,nt_intronic,real,splicing,intronic" << std::endl;
    for(auto& t : tpms_joined){
        joined_ss << std::get<0>(t).size() << "," << std::get<1>(t).size() << "," << std::get<2>(t).size() << ",";
        if(std::get<0>(t).size() == 0){
            joined_ss << "0;";
        }
        else{
            for(auto& s : std::get<0>(t)){
                joined_ss << s <<";";
            }
        }
        joined_ss.seekp(-1, std::ios_base::end);
        joined_ss << ",";
        if(std::get<1>(t).size() == 0){
            joined_ss << "0;";
        }
        else{
            for(auto& s : std::get<1>(t)){
                joined_ss << s <<";";
            }
        }
        joined_ss.seekp(-1, std::ios_base::end);
        joined_ss << ",";
        if(std::get<2>(t).size() == 0){
            joined_ss << "0;";
        }
        else{
            for(auto& s : std::get<2>(t)){
                joined_ss << s <<";";
            }
        }
        joined_ss.seekp(-1, std::ios_base::end);
        joined_ss << std::endl;
    }
    joined_ss.close();
}

void TrackingStats::save_sample_cov(std::string base_out_fname){
    std::string cov_sample_real_fname = base_out_fname+".cov_sample_real";
    std::ofstream cov_sample_real_ss(cov_sample_real_fname.c_str());
    cov_sample_real_ss<<"cov,fpkm,tpm"<<std::endl;
    for(int i=0;i<this->cov_sample_real.size();i++){
        cov_sample_real_ss<<cov_sample_real[i]<<","<<fpkm_sample_real[i]<<","<<tpm_sample_real[i]<<std::endl;
    }
    cov_sample_real_ss.close();

    std::string cov_sample_intronic_fname = base_out_fname+".cov_sample_intronic";
    std::ofstream cov_sample_intronic_ss(cov_sample_intronic_fname.c_str());
    cov_sample_intronic_ss<<"cov,fpkm,tpm"<<std::endl;
    for(int i=0;i<this->cov_sample_intronic.size();i++){
        cov_sample_intronic_ss<<cov_sample_intronic[i]<<","<<fpkm_sample_intronic[i]<<","<<tpm_sample_intronic[i]<<std::endl;
    }
    cov_sample_intronic_ss.close();

    std::string cov_sample_splicing_fname = base_out_fname+".cov_sample_splicing";
    std::ofstream cov_sample_splicing_ss(cov_sample_splicing_fname.c_str());
    cov_sample_splicing_ss<<"cov,fpkm,tpm"<<std::endl;
    for(int i=0;i<this->cov_sample_splicing.size();i++){
        cov_sample_splicing_ss<<cov_sample_splicing[i]<<","<<fpkm_sample_splicing[i]<<","<<tpm_sample_splicing[i]<<std::endl;
    }
    cov_sample_splicing_ss.close();

    std::string cov_sample_intergenic_fname = base_out_fname+".cov_sample_intergenic";
    std::ofstream cov_sample_intergenic_ss(cov_sample_intergenic_fname.c_str());
    cov_sample_intergenic_ss<<"cov,fpkm,tpm"<<std::endl;
    for(int i=0;i<this->cov_sample_intergenic.size();i++){
        cov_sample_intergenic_ss<<cov_sample_intergenic[i]<<","<<fpkm_sample_intergenic[i]<<","<<tpm_sample_intergenic[i]<<std::endl;
    }
    cov_sample_intergenic_ss.close();
}

void TrackingStats::save_sample_loc(std::string base_out_fname){
    std::string num_locs_fname = base_out_fname+".num_locs_sample";
    std::ofstream num_locs_ss(num_locs_fname);
    num_locs_ss<<"real,splicing,intronic,real_splicing,real_intronic,splicing_intronic,real_splicing_intronic,all_real,intergenic"<<std::endl;
    for(auto& num : num_locs_sample){
        num_locs_ss << num[0]<<","<<num[1] <<","<<num[2]<<","<<num[3]<<","<<num[4]<<","<<num[5] <<","<<num[6]<<","<<num[7]<<","<<num[8]<< std::endl;
    }
    num_locs_ss.close();
}

// we need to have distribution of loci per tissue and per sample as well to choose from
void TrackingStats::save_tissue_loc(std::string base_out_fname){
    std::string num_locs_fname = base_out_fname+".num_locs_tissue";
    std::ofstream num_locs_ss(num_locs_fname);
    num_locs_ss<<"real,intergenic"<<std::endl;
    for(auto& num : num_locs_tissue){
        num_locs_ss <<std::get<0>(num)<<","<< std::get<1>(num)<<","<<std::get<2>(num) << std::endl;
    }
    num_locs_ss.close();
}

void TrackingStats::save_sample_txs(std::string base_out_fname){
    std::string num_tx_per_sample_fname = base_out_fname+".num_tx_per_sample";
    std::ofstream num_tx_per_sample_ss(num_tx_per_sample_fname.c_str());
    num_tx_per_sample_ss<<"sample,real,splicing,intronic,intergenic"<<std::endl;

    for(auto& v : sample_txs){
        num_tx_per_sample_ss<<v.first<<","<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
    }

    num_tx_per_sample_ss.close();
}

void TrackingStats::save_tissue_tx_stats3(std::string base_out_fname){
    std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc3";
    std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
    num_tx_per_tissue_loc_ss<<"real,splicing,intronic,intergenic"<<std::endl;

    for(auto& v : num_tx_per_tissue_loc3){
        num_tx_per_tissue_loc_ss<<std::get<0>(v)<<","<<std::get<1>(v)<<","<<std::get<2>(v)<<","<<std::get<3>(v)<<std::endl;
    }

    num_tx_per_tissue_loc_ss.close();
}

void TrackingStats::save_sample_tx_stats3(std::string base_out_fname){
    std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc3";
    std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
    num_tx_per_sample_loc_ss<<"real,splicing,intronic,intergenic"<<std::endl;

    for(auto& v : num_tx_per_sample_loc3){
        num_tx_per_sample_loc_ss<<std::get<0>(v)<<","<<std::get<1>(v)<<","<<std::get<2>(v)<<","<<std::get<3>(v)<<std::endl;
    }

    num_tx_per_sample_loc_ss.close();
}

void TrackingStats::save_sample_tx_stats4(std::string base_out_fname){
    std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc4";
    std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
    num_tx_per_sample_loc_ss<<"real,splicing,intronic,intergenic,total_real,total_splicing,total_intronic,total_intergenic"<<std::endl;

    for(auto& v : num_tx_per_sample_loc4){
        sample2tissue_loc_txs_it.first = sample2tissue_loc_txs.find(std::get<4>(v.second));
        if(sample2tissue_loc_txs_it.first==sample2tissue_loc_txs.end()){
            std::cerr<<"non-existing locus"<<std::endl;
            exit(-1);
        }
        num_tx_per_sample_loc_ss<<std::get<0>(v.second)<<","
                                <<std::get<1>(v.second)<<","
                                <<std::get<2>(v.second)<<","
                                <<std::get<3>(v.second)<<","
                                <<sample2tissue_loc_txs_it.first->second[0]<<","
                                <<sample2tissue_loc_txs_it.first->second[1]<<","
                                <<sample2tissue_loc_txs_it.first->second[2]<<","
                                <<sample2tissue_loc_txs_it.first->second[3]<<std::endl;
    }

    num_tx_per_sample_loc_ss.close();
}

void TrackingStats::save_tissue_tx_stats5(std::string base_out_fname){
    std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc5";
    std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
    num_tx_per_tissue_loc_ss<<"real,splicing,intronic,intergenic"<<std::endl;
    for(auto& v : num_tx_per_tissue_loc5){
        num_tx_per_tissue_loc_ss<<v.first<<","<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
    }
    num_tx_per_tissue_loc_ss.close();
}

void TrackingStats::save_sample_tpm_gauss(std::string base_out_fname){
    std::string gauss_sample_tpm_per_tissue_loc_fname = base_out_fname+".tpm_sample_gauss";
    std::ofstream gauss_sample_tpm_per_tissue_loc_ss(gauss_sample_tpm_per_tissue_loc_fname.c_str());
    gauss_sample_tpm_per_tissue_loc_ss<<"real_mean,real_sd,real_num,splicing_mean,splicing_sd,splicing_num,intronic_mean,intronic_sd,intronic_num,intergenic_mean,intergenic_sd,intergenic_num";
    for(int i=0;i<gauss_sample_tpm_per_tissue_loc.size();i++){
        if(i%4==0){ // seen real,splice,int and pol and can output CR
            gauss_sample_tpm_per_tissue_loc_ss<<std::endl;
        }
        else{
            gauss_sample_tpm_per_tissue_loc_ss<<",";
        }
        int total_num = 0;
        for(auto& t : gauss_sample_tpm_per_tissue_loc[i]){
            if(t.first>0){ // is not empty
                total_num++;
            }
            gauss_sample_tpm_per_tissue_loc_ss<<t.first<<";";
        }
        gauss_sample_tpm_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
        gauss_sample_tpm_per_tissue_loc_ss <<",";
        for(auto& t : gauss_sample_tpm_per_tissue_loc[i]){
            gauss_sample_tpm_per_tissue_loc_ss<<t.second<<";";
        }
        gauss_sample_tpm_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
        gauss_sample_tpm_per_tissue_loc_ss <<","<<total_num;
    }
    gauss_sample_tpm_per_tissue_loc_ss<<std::endl;

    gauss_sample_tpm_per_tissue_loc_ss.close();
}

void TrackingStats::save_sample_tx_gauss(std::string base_out_fname){
    std::string gauss_sample_tx_per_tissue_loc_fname = base_out_fname+".num_tx_sample_gauss";
    std::ofstream gauss_sample_tx_per_tissue_loc_ss(gauss_sample_tx_per_tissue_loc_fname.c_str());
    gauss_sample_tx_per_tissue_loc_ss<<"real_mean,real_sd,real_num,splicing_mean,splicing_sd,splicing_num,intronic_mean,intronic_sd,intronic_num,intergenic_mean,intergenic_sd,intergenic_num";
    for(int i=0;i<gauss_sample_tx_per_tissue_loc.size();i++){
        if(i%4==0){ // seen real,splice,int and pol and can output CR
            gauss_sample_tx_per_tissue_loc_ss<<std::endl;
        }
        else{
            gauss_sample_tx_per_tissue_loc_ss<<",";
        }
        gauss_sample_tx_per_tissue_loc_ss<<std::get<0>(gauss_sample_tx_per_tissue_loc[i])<<","<<std::get<1>(gauss_sample_tx_per_tissue_loc[i])<<","<<std::get<2>(gauss_sample_tx_per_tissue_loc[i]);
    }
    gauss_sample_tx_per_tissue_loc_ss<<std::endl;

    gauss_sample_tx_per_tissue_loc_ss.close();
}

void TrackingStats::save_sample_gauss(std::string base_out_fname){
    std::string gauss_sample_per_tissue_loc_fname = base_out_fname+".sample_gauss";
    std::ofstream gauss_sample_per_tissue_loc_ss(gauss_sample_per_tissue_loc_fname.c_str());
    gauss_sample_per_tissue_loc_ss<<"real_mean,real_sd,real_num,total_real_num,real_mean_tpm,real_sd_tpm,real_tpms,splicing_mean,splicing_sd,splicing_num,total_splicing_num,splicing_mean_tpm,splicing_sd_tpm,splicing_tpms,intronic_mean,intronic_sd,intronic_num,total_intronic_num,intronic_mean_tpm,intronic_sd_tpm,intronic_tpms,intergenic_mean,intergenic_sd,intergenic_num,total_intergenic_num,intergenic_mean_tpm,intergenic_sd_tpm,intergenic_tpms";
    for(int i=0;i<gauss_sample_per_tissue_loc.size();i++){
        if(i%4==0){ // seen real,splice,int and pol and can output CR
            gauss_sample_per_tissue_loc_ss<<std::endl;
        }
        else{
            gauss_sample_per_tissue_loc_ss<<",";
        }
        gauss_sample_per_tissue_loc_ss<<std::get<0>(gauss_sample_per_tissue_loc[i].first)<<","<<std::get<1>(gauss_sample_per_tissue_loc[i].first)<<","<<std::get<2>(gauss_sample_per_tissue_loc[i].first)<<","<<std::get<3>(gauss_sample_per_tissue_loc[i].first)<<",";

        for(auto& t : gauss_sample_per_tissue_loc[i].second){
            gauss_sample_per_tissue_loc_ss<<std::get<0>(t)<<";";
        }
        gauss_sample_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
        gauss_sample_per_tissue_loc_ss <<",";
        for(auto& t : gauss_sample_per_tissue_loc[i].second){
            gauss_sample_per_tissue_loc_ss<<std::get<1>(t)<<";";
        }
        gauss_sample_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
        gauss_sample_per_tissue_loc_ss <<",";
        for(auto& t : gauss_sample_per_tissue_loc[i].second){
            if(std::get<2>(t).size()==0){
                gauss_sample_per_tissue_loc_ss<<0<<":";
            }
            else{
                for(auto& tpm_val : std::get<2>(t)){
                    if(tpm_val==0){
                        gauss_sample_per_tissue_loc_ss<<0<<":";
                    }
                    else{
                        gauss_sample_per_tissue_loc_ss<<FIXED_FLOAT(tpm_val)<<":";
                    }
                }
                gauss_sample_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
                gauss_sample_per_tissue_loc_ss<<";";
            }
        }
        gauss_sample_per_tissue_loc_ss.unsetf(std::ios_base::fixed);
        gauss_sample_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
    }
    gauss_sample_per_tissue_loc_ss<<std::endl;

    gauss_sample_per_tissue_loc_ss.close();
}

void TrackingStats::save_sample_tx_stats6(std::string base_out_fname) {
    std::string num_tx_per_sample_loc_fname = base_out_fname + ".num_tx_per_sample_loc6";
    std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
    num_tx_per_sample_loc_ss<<"total_real,total_splicing,total_intronic,total_intergenic,"<<
                            "tissue_real,tissue_splicing,tissue_intronic,tissue_intergenic,"<<
                            "tpms"<<std::endl;

    for(auto& tl : num_tx_per_sample_loc6){ // tissue + locus
        num_tx_per_sample_loc_ss<<std::get<0>(tl.second)[0]<<","<<std::get<0>(tl.second)[1]<<","<<std::get<0>(tl.second)[2]<<","<<std::get<0>(tl.second)[3]<<","; // save total number of transcript per ALL locus
        num_tx_per_sample_loc_ss<<std::get<1>(tl.second)[0].size()<<","<<std::get<1>(tl.second)[1].size()<<","<<std::get<1>(tl.second)[2].size()<<","<<std::get<1>(tl.second)[3].size()<<","; // save the number of transcript per current tissue locus
        for(auto& sample : std::get<2>(tl.second)){ // for each sample
            if(sample.second[0].size()==0){
                num_tx_per_sample_loc_ss<<"0"<<"-";
            }
            for(auto& tpm : sample.second[0]){ // for each real tpm
                if(tpm==0){
                    num_tx_per_sample_loc_ss<<0<<"-";
                }
                else{
                    num_tx_per_sample_loc_ss<<FIXED_FLOAT(tpm)<<"-";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss<<":"; // delimiter between types within sample

            if(sample.second[1].size()==0){
                num_tx_per_sample_loc_ss<<"0"<<"-";
            }
            for(auto& tpm : sample.second[1]){ // for each splicing tpm
                if(tpm==0){
                    num_tx_per_sample_loc_ss<<0<<"-";
                }
                else{
                    num_tx_per_sample_loc_ss<<FIXED_FLOAT(tpm)<<"-";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss<<":"; // delimiter between types within sample

            if(sample.second[2].size()==0){
                num_tx_per_sample_loc_ss<<"0"<<"-";
            }
            for(auto& tpm : sample.second[2]){ // for each intronic tpm
                if(tpm==0){
                    num_tx_per_sample_loc_ss<<0<<"-";
                }
                else{
                    num_tx_per_sample_loc_ss<<FIXED_FLOAT(tpm)<<"-";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss<<":"; // delimiter between types within sample

            if(sample.second[3].size()==0){
                num_tx_per_sample_loc_ss<<"0"<<"-";
            }
            for(auto& tpm : sample.second[3]){ // for each intergenic tpm
                if(tpm==0){
                    num_tx_per_sample_loc_ss<<0<<"-";
                }
                else{
                    num_tx_per_sample_loc_ss<<FIXED_FLOAT(tpm)<<"-";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss<<";"; // delimiter between samples
        }
        num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
        num_tx_per_sample_loc_ss<<std::endl; // new line
    }
    num_tx_per_sample_loc_ss.close();
}

void TrackingStats::save_sample_tx_stats5(std::string base_out_fname){
    std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc5";
    std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
    num_tx_per_sample_loc_ss<<"sample,real,splicing,intronic,intergenic,tpms_real,total_tpm_real,tpms_splicing,total_tpm_splicing,tpms_intronic,total_tpm_intronic,tpms_intergenic,total_tpm_intergenic"<<std::endl;

    for(auto& v : num_tx_per_sample_loc5){
        num_tx_per_sample_loc_ss<<v.first.first<<","<<std::get<0>(v.second).size()<<","<<std::get<1>(v.second).size()<<","<<std::get<2>(v.second).size()<<","<<std::get<3>(v.second).size()<<",";
        float sum_real = 0;
        if(std::get<0>(v.second).size() == 0){
            num_tx_per_sample_loc_ss << "0;";
        }
        else{
            for(auto& s : std::get<0>(v.second)){
                sum_real+=s;
                num_tx_per_sample_loc_ss << s <<";";
            }
        }
        num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
        num_tx_per_sample_loc_ss << ","<<sum_real<<",";
        float sum_splicing = 0;
        if(std::get<1>(v.second).size() == 0){
            num_tx_per_sample_loc_ss << "0;";
        }
        else{
            for(auto& s : std::get<1>(v.second)){
                sum_splicing+=s;
                num_tx_per_sample_loc_ss << s <<";";
            }
        }
        num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
        num_tx_per_sample_loc_ss << ","<<sum_splicing<<",";
        float sum_intronic = 0;
        if(std::get<2>(v.second).size() == 0){
            num_tx_per_sample_loc_ss << "0;";
        }
        else{
            for(auto& s : std::get<2>(v.second)){
                sum_intronic+=s;
                num_tx_per_sample_loc_ss << s <<";";
            }
        }
        num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
        num_tx_per_sample_loc_ss << ","<<sum_intronic<<",";
        float sum_intergenic = 0;
        if(std::get<3>(v.second).size() == 0){
            num_tx_per_sample_loc_ss << "0;";
        }
        else{
            for(auto& s : std::get<3>(v.second)){
                sum_intergenic+=s;
                num_tx_per_sample_loc_ss << s <<";";
            }
        }
        num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
        num_tx_per_sample_loc_ss<<","<<sum_intergenic;
        num_tx_per_sample_loc_ss << std::endl;
    }

    num_tx_per_sample_loc_ss.close();
}

void TrackingStats::save_loc_stats(std::string base_out_fname){
    std::string loc_fname = base_out_fname+".loc_stats";
    std::ofstream loc_ss(loc_fname.c_str());
    loc_ss<<"real: "<<num_avg_real_locs_per_tissue<<std::endl;
    loc_ss<<"noise: "<<num_avg_noise_locs_per_tissue<<std::endl;
    loc_ss<<"undefined: "<<num_avg_undef_locs_per_tissue<<std::endl;
    loc_ss.close();
}

void TrackingStats::save_num_tx_tissue(std::string base_out_fname){
    std::string num_tx_per_tissue_fname = base_out_fname+".num_tx_per_tissue";
    std::ofstream num_tx_per_tissue_ss(num_tx_per_tissue_fname.c_str());
    num_tx_per_tissue_ss<<"real,splicing,intronic,intergenic"<<std::endl;
    for(auto& v : num_tx_per_tissue){
        num_tx_per_tissue_ss<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
    }
    num_tx_per_tissue_ss.close();
}

void TrackingStats::save_tx_freq(std::string base_out_fname){
    std::string tissue_tx_freq_fname = base_out_fname+".tissue_tx_freq";
    std::ofstream tissue_tx_freq_ss(tissue_tx_freq_fname.c_str());
    tissue_tx_freq_ss<<"type,freq"<<std::endl;
    for(auto& v : this->tissue_tx_freq_real){
        tissue_tx_freq_ss<<"real"<<","<<v<<std::endl;
    }
    for(auto& v : this->tissue_tx_freq_intronic){
        tissue_tx_freq_ss<<"intronic"<<","<<v<<std::endl;
    }
    for(auto& v : this->tissue_tx_freq_splicing){
        tissue_tx_freq_ss<<"splicing"<<","<<v<<std::endl;
    }
    for(auto& v : this->tissue_tx_freq_intergenic){
        tissue_tx_freq_ss<<"intergenic"<<","<<v<<std::endl;
    }
    tissue_tx_freq_ss.close();
}
