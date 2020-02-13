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

#include "gff.h"
#include "GFaSeqGet.h"
#include "arg_parse.h"

// Types of transcripts
enum TYPE {REAL_TX        = 1,
           SPLICING_TX = 2,
           INTRONIC_TX    = 3,
           POLYMERASE_TX  = 4,
    UNDEFINED_TX      = -1};

// Minimum Transcript
struct MT{
    float cov; // coverage
    float tpm; // tpm
    float fpkm; // fpkm
    std::string tid; // transcript id
    std::string locus; // locus id
    std::string tissue;
    std::string sample;
    int type = TYPE::UNDEFINED_TX;
    MT(std::string tid,std::string locus,float fpkm,float tpm,float cov,std::string tissue,std::string sample):tid(tid), locus(locus), fpkm(fpkm), tpm(tpm), cov(cov), tissue(tissue),sample(sample) { }
    void set_cov(int cov){cov=cov;}
    void set_tpm(int tpm){tpm=tpm;}
    std::string get_strg(){
        return tissue+"\t"+sample+"\t"+tid+":"+locus+"|"+std::to_string(fpkm)+"|"+std::to_string(tpm)+"|"+std::to_string(cov);
    }
    void set_type(int type){
        if(this->type!=TYPE::UNDEFINED_TX && this->type!=type){
            std::cerr<<"repeating types detected on sample level"<<std::endl;
            exit(-1);
        }
        this->type=type;
    }
};

typedef std::map<std::string,MT> MTM; // map of transcripts per sample
typedef std::pair<MTM::iterator,bool> MTM_IT;
typedef std::map<std::string,std::vector<MTM_IT>> Loci_Sample;
typedef std::pair<Loci_Sample::iterator,bool> Loci_Sample_IT;

// Minimum transcript to tissue
struct MTT{
    std::string tid;
    std::string locus;
    std::string tissue;
    std::set<std::string> samples;
    std::vector<MTM_IT> txs;
    int type = TYPE::UNDEFINED_TX;
    MTT(std::string tid,std::string locus,std::string tissue):tid(tid),locus(locus),tissue(tissue){}
    void add_tx(MTM_IT mtm_it){
        this->txs.push_back(mtm_it);
        this->samples.insert(mtm_it.first->second.sample);
    }

    std::string get_strg(){
        std::string res = "";
        for(auto& v : txs){
            res+=tid+"\t"+locus+"\t"+v.first->second.get_strg()+"\n";
        }
        return res;
    }

    int get_num_txs(){
        return txs.size();
    }

    float get_sum_tpms(){
        float res=0;
        for(auto& t : txs){
            res+=t.first->second.tpm;
        }

        return res;
    }

    void set_type(int type){
        if(this->type!=TYPE::UNDEFINED_TX && this->type!=type){
            std::cerr<<"repeating types detected on tissue level"<<std::endl;
            exit(-1);
        }
        this->type=type;
        // propagate down
        for(auto& t : txs){
            t.first->second.set_type(type);
        }
    }
};

typedef std::map<std::string,MTT> MTTM; // map of tissue merged transcripts to sample transcripts
typedef std::pair<MTTM::iterator,bool> MTTM_IT;
typedef std::set<std::string> Tissues;
typedef std::pair<std::set<std::string>::iterator,bool> Tissues_IT;
typedef std::map<std::string,std::vector<MTTM_IT>> Loci_Tissue;
typedef std::pair<Loci_Tissue::iterator,bool> Loci_Tissue_IT;

struct MATT{
    std::string tid;
    std::string locus;
    std::vector<MTTM_IT> txs;
    std::set<std::string> tissues;
    int type = TYPE::UNDEFINED_TX;
    MATT(std::string tid,std::string loc):tid(tid),locus(loc){}
    void add_tx(MTTM_IT mttm_it){
        this->txs.push_back(mttm_it);
        this->tissues.insert(mttm_it.first->second.tissue);
    }
    int num_tissues(){
        return this->tissues.size();
    }
    std::string get_strg(){
        std::string res = "";
        for(auto& v : txs){
            res+=tid+"\t"+locus+"\t"+v.first->second.get_strg()+"\n";
        }

        return res;
    }

    int get_num_sample_txs(){
        int res=0;
        for(auto& v : txs){
            res += v.first->second.get_num_txs();
        }

        return res;
    }

    float get_sum_sample_tpms(){
        float res=0;
        for(auto& t : txs){
            res += t.first->second.get_sum_tpms();
        }

        return res;
    }

    int set_type(int type){
        if(this->type!=TYPE::UNDEFINED_TX){
            std::cerr<<"Type was already set for: "<<tid<<std::endl;
            exit(-1);
        }
        this->type=type;
        // now propagate down to the tissue and sample levels
        for(auto& v : txs){
            v.first->second.set_type(type);
        }
    }
};

typedef std::map<std::string,MATT> MATTM; // map of all merged transcripts to tissue merged transcripts
typedef std::pair<MATTM::iterator,bool> MATTM_IT;
typedef std::map<std::string,std::pair<std::vector<MATTM_IT>,int>> Loci; // value contains 1. iterators to transcripts; 2. type of locus (real,noise)
typedef std::pair<Loci::iterator,bool> Loci_IT;

struct TR{
    std::string new_tid,new_loc,old_tid,old_loc,qj,fpkm,tpm,cov;
    TR() = default;
    TR(std::string new_tid,std::string new_loc,std::string old_tid,std::string old_loc,std::string qj,std::string fpkm,std::string tpm,std::string cov):new_tid(new_tid),new_loc(new_loc),old_tid(old_tid),old_loc(old_loc),qj(qj),fpkm(fpkm),tpm(tpm),cov(cov){}
};

struct Stats{
    std::vector<int> num_tx_per_all_loc,num_tx_per_all_loc_real,num_tx_per_all_loc_intronic,num_tx_per_all_loc_splicing,num_tx_per_all_loc_polymerase,num_tx_per_all_loc_left;
    std::vector<float> sum_tx_per_all_loc,sum_tx_per_all_loc_real,sum_tx_per_all_loc_intronic,sum_tx_per_all_loc_splicing,sum_tx_per_all_loc_polymerase,sum_tx_per_all_loc_left;
    void save_all_tx_stats(std::string base_out_fname){
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
        else if(num_tx_per_all_loc_splicing.size() != num_tx_per_all_loc_polymerase.size()) {
            std::cerr<<"all loc vectors not equal size"<<std::endl;
            std::cerr << "num_tx_per_all_loc_splicing: " << num_tx_per_all_loc_splicing.size() << std::endl;
            std::cerr<<"num_tx_per_all_loc_polymerase: "<<num_tx_per_all_loc_polymerase.size()<<std::endl;
            exit(-1);
        }
        else if(num_tx_per_all_loc_polymerase.size()!=num_tx_per_all_loc_left.size()) {
            std::cerr<<"all loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_all_loc_polymerase: "<<num_tx_per_all_loc_polymerase.size()<<std::endl;
            std::cerr<<"num_tx_per_all_loc_left: "<<num_tx_per_all_loc_left.size()<<std::endl;
            exit(-1);
        }
        else{
            std::string num_tx_per_all_loc_fname = base_out_fname+".num_tx_per_all_loc";
            std::ofstream num_tx_per_all_loc_ss(num_tx_per_all_loc_fname.c_str());
            num_tx_per_all_loc_ss<<"all,real,intronic,splicing,polymerase,left"<<std::endl;

            for(int i=0;i<this->num_tx_per_all_loc.size();i++){
                num_tx_per_all_loc_ss << num_tx_per_all_loc[i] << "," << num_tx_per_all_loc_real[i] << "," << num_tx_per_all_loc_intronic[i] << "," << num_tx_per_all_loc_splicing[i] << "," << num_tx_per_all_loc_polymerase[i] << "," << num_tx_per_all_loc_left[i] << std::endl;
            }

            num_tx_per_all_loc_ss.close();

            std::string sum_tx_per_all_loc_fname = base_out_fname+".sum_tx_per_all_loc";
            std::ofstream sum_tx_per_all_loc_ss(sum_tx_per_all_loc_fname.c_str());
            sum_tx_per_all_loc_ss<<"all,real,intronic,splicing,polymerase,left"<<std::endl;

            for(int i=0;i<this->sum_tx_per_all_loc.size();i++){
                sum_tx_per_all_loc_ss << sum_tx_per_all_loc[i] << "," << sum_tx_per_all_loc_real[i] << "," << sum_tx_per_all_loc_intronic[i] << "," << sum_tx_per_all_loc_splicing[i] << "," << sum_tx_per_all_loc_polymerase[i] << "," << sum_tx_per_all_loc_left[i] << std::endl;
            }

            sum_tx_per_all_loc_ss.close();
        }
    }

    std::vector<int> num_tx_per_tissue_loc,num_tx_per_tissue_loc_real,num_tx_per_tissue_loc_intronic,num_tx_per_tissue_loc_splicing,num_tx_per_tissue_loc_polymerase,num_tx_per_tissue_loc_left;
    std::vector<float> sum_tx_per_tissue_loc,sum_tx_per_tissue_loc_real,sum_tx_per_tissue_loc_intronic,sum_tx_per_tissue_loc_splicing,sum_tx_per_tissue_loc_polymerase,sum_tx_per_tissue_loc_left;
    void save_tissue_tx_stats(std::string base_out_fname){
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
        else if(num_tx_per_tissue_loc_splicing.size() != num_tx_per_tissue_loc_polymerase.size()) {
            std::cerr<<"tissue loc vectors not equal size"<<std::endl;
            std::cerr << "num_tx_per_tissue_loc_splicing: " << num_tx_per_tissue_loc_splicing.size() << std::endl;
            std::cerr<<"num_tx_per_tissue_loc_polymerase: "<<num_tx_per_tissue_loc_polymerase.size()<<std::endl;
            exit(-1);
        }
        else if(num_tx_per_tissue_loc_polymerase.size()!=num_tx_per_tissue_loc_left.size()) {
            std::cerr<<"tissue loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_tissue_loc_polymerase: "<<num_tx_per_tissue_loc_polymerase.size()<<std::endl;
            std::cerr<<"num_tx_per_tissue_loc_left: "<<num_tx_per_tissue_loc_left.size()<<std::endl;
            exit(-1);
        }
        else{
            std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc";
            std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
            num_tx_per_tissue_loc_ss<<"all,real,intronic,splicing,polymerase,left"<<std::endl;

            for(int i=0;i<this->num_tx_per_tissue_loc.size();i++){
                num_tx_per_tissue_loc_ss << num_tx_per_tissue_loc[i] << "," << num_tx_per_tissue_loc_real[i] << "," << num_tx_per_tissue_loc_intronic[i] << "," << num_tx_per_tissue_loc_splicing[i] << "," << num_tx_per_tissue_loc_polymerase[i] << "," << num_tx_per_tissue_loc_left[i] << std::endl;
            }

            num_tx_per_tissue_loc_ss.close();

            std::string sum_tx_per_tissue_loc_fname = base_out_fname+".sum_tx_per_tissue_loc";
            std::ofstream sum_tx_per_tissue_loc_ss(sum_tx_per_tissue_loc_fname.c_str());
            sum_tx_per_tissue_loc_ss<<"all,real,intronic,splicing,polymerase,left"<<std::endl;

            for(int i=0;i<this->sum_tx_per_tissue_loc.size();i++){
                sum_tx_per_tissue_loc_ss << sum_tx_per_tissue_loc[i] << "," << sum_tx_per_tissue_loc_real[i] << "," << sum_tx_per_tissue_loc_intronic[i] << "," << sum_tx_per_tissue_loc_splicing[i] << "," << sum_tx_per_tissue_loc_polymerase[i] << "," << sum_tx_per_tissue_loc_left[i] << std::endl;
            }

            sum_tx_per_tissue_loc_ss.close();
        }
    }

    std::vector<int> num_tx_per_sample_loc,num_tx_per_sample_loc_real,num_tx_per_sample_loc_intronic,num_tx_per_sample_loc_splicing,num_tx_per_sample_loc_polymerase,num_tx_per_sample_loc_left;
    std::vector<float> sum_tx_per_sample_loc,sum_tx_per_sample_loc_real,sum_tx_per_sample_loc_intronic,sum_tx_per_sample_loc_splicing,sum_tx_per_sample_loc_polymerase,sum_tx_per_sample_loc_left;
    void save_sample_tx_stats(std::string base_out_fname){
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
        else if(num_tx_per_sample_loc_splicing.size() != num_tx_per_sample_loc_polymerase.size()) {
            std::cerr<<"sample loc vectors not equal size"<<std::endl;
            std::cerr << "num_tx_per_sample_loc_splicing: " << num_tx_per_sample_loc_splicing.size() << std::endl;
            std::cerr<<"num_tx_per_sample_loc_polymerase: "<<num_tx_per_sample_loc_polymerase.size()<<std::endl;
            exit(-1);
        }
        else if(num_tx_per_sample_loc_polymerase.size()!=num_tx_per_sample_loc_left.size()) {
            std::cerr<<"sample loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_sample_loc_polymerase: "<<num_tx_per_sample_loc_polymerase.size()<<std::endl;
            std::cerr<<"num_tx_per_sample_loc_left: "<<num_tx_per_sample_loc_left.size()<<std::endl;
            exit(-1);
        }
        else{
            std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc";
            std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
            num_tx_per_sample_loc_ss<<"all,real,intronic,splicing,polymerase,left"<<std::endl;

            for(int i=0;i<this->num_tx_per_sample_loc.size();i++){
                num_tx_per_sample_loc_ss << num_tx_per_sample_loc[i] << "," << num_tx_per_sample_loc_real[i] << "," << num_tx_per_sample_loc_intronic[i] << "," << num_tx_per_sample_loc_splicing[i] << "," << num_tx_per_sample_loc_polymerase[i] << "," << num_tx_per_sample_loc_left[i] << std::endl;
            }

            num_tx_per_sample_loc_ss.close();

            std::string sum_tx_per_sample_loc_fname = base_out_fname+".sum_tx_per_sample_loc";
            std::ofstream sum_tx_per_sample_loc_ss(sum_tx_per_sample_loc_fname.c_str());
            sum_tx_per_sample_loc_ss<<"all,real,intronic,splicing,polymerase,left"<<std::endl;

            for(int i=0;i<this->sum_tx_per_sample_loc.size();i++){
                sum_tx_per_sample_loc_ss << sum_tx_per_sample_loc[i] << "," << sum_tx_per_sample_loc_real[i] << "," << sum_tx_per_sample_loc_intronic[i] << "," << sum_tx_per_sample_loc_splicing[i] << "," << sum_tx_per_sample_loc_polymerase[i] << "," << sum_tx_per_sample_loc_left[i] << std::endl;
            }

            sum_tx_per_sample_loc_ss.close();
        }
    }

    std::vector<int> num_tx_per_tissue_loc2,num_tx_per_tissue_loc_real2,num_tx_per_tissue_loc_intronic2,num_tx_per_tissue_loc_splicing2,num_tx_per_tissue_loc_polymerase2,num_tx_per_tissue_loc_left2;
    void save_tissue_tx_stats2(std::string base_out_fname){
        std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc2";
        std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
        num_tx_per_tissue_loc_ss<<"total,real,intronic,splicing,polymerase,left"<<std::endl;

        for(int i=0;i<this->num_tx_per_tissue_loc2.size();i++){
            num_tx_per_tissue_loc_ss<<num_tx_per_tissue_loc2[i]<<","<<num_tx_per_tissue_loc_real2[i]<<","<<num_tx_per_tissue_loc_intronic2[i]<<","<<num_tx_per_tissue_loc_splicing2[i]<<","<<num_tx_per_tissue_loc_polymerase2[i]<<","<<num_tx_per_tissue_loc_left2[i]<<std::endl;
        }

        num_tx_per_tissue_loc_ss.close();
    }

    std::vector<int> num_tx_per_sample_loc2,num_tx_per_sample_loc_real2,num_tx_per_sample_loc_intronic2,num_tx_per_sample_loc_splicing2,num_tx_per_sample_loc_polymerase2,num_tx_per_sample_loc_left2;
    void save_sample_tx_stats2(std::string base_out_fname){
        std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc2";
        std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
        num_tx_per_sample_loc_ss<<"total,real,intronic,splicing,polymerase,left"<<std::endl;

        for(int i=0;i<this->num_tx_per_sample_loc2.size();i++){
            num_tx_per_sample_loc_ss<<num_tx_per_sample_loc2[i]<<","<<num_tx_per_sample_loc_real2[i]<<","<<num_tx_per_sample_loc_intronic2[i]<<","<<num_tx_per_sample_loc_splicing2[i]<<","<<num_tx_per_sample_loc_polymerase2[i]<<","<<num_tx_per_sample_loc_left2[i]<<std::endl;
        }

        num_tx_per_sample_loc_ss.close();
    }

    // for real loci only - used to compute fraction of real to splicing transcripts per locus
    std::vector<std::tuple<int,float,int,float,int,float>> cov_fracs; // number of real transcripts, sum of tpms of real, number of splicing transcripts, sum of tpms of splicing, num of intronic transcirpts, sum of intronic
    void save_cov_fracs(std::string base_out_fname){
        std::string frac_fname = base_out_fname+".cov_fracs";
        std::ofstream frac_ss(frac_fname.c_str());
        frac_ss << "nt_real,sum_real,nt_nonint,sum_nonint,nt_int,sum_int" << std::endl;
        for(auto& t : cov_fracs){
            frac_ss << std::get<0>(t) << "," << std::get<1>(t) << "," << std::get<2>(t) << "," << std::get<3>(t) << "," << std::get<4>(t) << "," << std::get<5>(t) << std::endl;
        }
        frac_ss.close();
    }

    std::vector<std::tuple<std::vector<float>,std::vector<float>,std::vector<float>>> covs_joined;
    void save_cov_joined(std::string base_out_fname){
        std::string joined_fname = base_out_fname+".covs_joined";
        std::ofstream joined_ss(joined_fname.c_str());
        joined_ss << "nt_real,nt_nonint,nt_int,real,nonint,int" << std::endl;
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

    std::vector<std::tuple<int,float,int,float,int,float>> tpm_fracs; // number of real transcripts, sum of tpms of real, number of splicing transcripts, sum of tpms of splicing, num intronic, sum intronic
    void save_tpm_fracs(std::string base_out_fname){
        std::string frac_fname = base_out_fname+".tpm_fracs";
        std::ofstream frac_ss(frac_fname.c_str());
        frac_ss << "nt_real,sum_real,nt_nonint,sum_nonint,nt_int,sum_int" << std::endl;
        for(auto& t : tpm_fracs){
            frac_ss << std::get<0>(t) << "," << std::get<1>(t) << "," << std::get<2>(t) << "," << std::get<3>(t) << "," << std::get<4>(t) << "," << std::get<5>(t) << std::endl;
        }
        frac_ss.close();
    }

    std::vector<std::tuple<std::vector<float>,std::vector<float>,std::vector<float>>> tpms_joined;
    void save_tpm_joined(std::string base_out_fname){
        std::string joined_fname = base_out_fname+".tpms_joined";
        std::ofstream joined_ss(joined_fname.c_str());
        joined_ss << "nt_real,nt_nonint,nt_int,real,nonint,int" << std::endl;
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

    std::vector<float> cov_sample_real,cov_sample_intronic,cov_sample_splicing,cov_sample_polymerase,cov_sample_left;
    std::vector<float> fpkm_sample_real,fpkm_sample_intronic,fpkm_sample_splicing,fpkm_sample_polymerase,fpkm_sample_left;
    std::vector<float> tpm_sample_real,tpm_sample_intronic,tpm_sample_splicing,tpm_sample_polymerase,tpm_sample_left;
    void save_sample_cov(std::string base_out_fname){
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

        std::string cov_sample_polymerase_fname = base_out_fname+".cov_sample_polymerase";
        std::ofstream cov_sample_polymerase_ss(cov_sample_polymerase_fname.c_str());
        cov_sample_polymerase_ss<<"cov,fpkm,tpm"<<std::endl;
        for(int i=0;i<this->cov_sample_polymerase.size();i++){
            cov_sample_polymerase_ss<<cov_sample_polymerase[i]<<","<<fpkm_sample_polymerase[i]<<","<<tpm_sample_polymerase[i]<<std::endl;
        }
        cov_sample_polymerase_ss.close();
    }

    std::vector<int> num_real_locs_sample,num_noise_locs_sample;
    void save_sample_loc(std::string base_out_fname){
        std::string num_real_locs_fname = base_out_fname+".num_real_locs_sample";
        std::ofstream num_real_locs_ss(num_real_locs_fname);
        for(auto& num : num_real_locs_sample){
            num_real_locs_ss << num << std::endl;
        }
        num_real_locs_ss.close();

        std::string num_noise_locs_fname = base_out_fname+".num_noise_locs_sample";
        std::ofstream num_noise_locs_ss(num_noise_locs_fname);
        for(auto& num : num_noise_locs_sample){
            num_noise_locs_ss << num << std::endl;
        }
        num_noise_locs_ss.close();
    }

    // we need to have distribution of loci per tissue and per sample as well to choose from
    std::vector<int> num_real_locs_tissue,num_noise_locs_tissue;
    void save_tissue_loc(std::string base_out_fname){
        std::string num_real_locs_fname = base_out_fname+".num_real_locs_tissue";
        std::ofstream num_real_locs_ss(num_real_locs_fname);
        for(auto& num : num_real_locs_tissue){
            num_real_locs_ss << num << std::endl;
        }
        num_real_locs_ss.close();

        std::string num_noise_locs_fname = base_out_fname+".num_noise_locs_tissue";
        std::ofstream num_noise_locs_ss(num_noise_locs_fname);
        for(auto& num : num_noise_locs_tissue){
            num_noise_locs_ss << num << std::endl;
        }
        num_noise_locs_ss.close();
    }

    std::map<std::string,std::tuple<int,int,int,int>> sample_txs;
    void save_sample_txs(std::string base_out_fname){
        std::string num_tx_per_sample_fname = base_out_fname+".num_tx_per_sample";
        std::ofstream num_tx_per_sample_ss(num_tx_per_sample_fname.c_str());
        num_tx_per_sample_ss<<"sample,real,intronic,splicing,polymerase"<<std::endl;

        for(auto& v : sample_txs){
            num_tx_per_sample_ss<<v.first<<","<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
        }

        num_tx_per_sample_ss.close();
    }

    std::vector<std::tuple<int,int,int,int>> num_tx_per_tissue_loc3;
    void save_tissue_tx_stats3(std::string base_out_fname){
        std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc3";
        std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
        num_tx_per_tissue_loc_ss<<"real,intronic,splicing,polymerase"<<std::endl;

        for(auto& v : num_tx_per_tissue_loc3){
            num_tx_per_tissue_loc_ss<<std::get<0>(v)<<","<<std::get<1>(v)<<","<<std::get<2>(v)<<","<<std::get<3>(v)<<std::endl;
        }

        num_tx_per_tissue_loc_ss.close();
    }

    std::vector<std::tuple<int,int,int,int>> num_tx_per_sample_loc3;
    void save_sample_tx_stats3(std::string base_out_fname){
        std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc3";
        std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
        num_tx_per_sample_loc_ss<<"real,intronic,splicing,polymerase"<<std::endl;

        for(auto& v : num_tx_per_sample_loc3){
            num_tx_per_sample_loc_ss<<std::get<0>(v)<<","<<std::get<1>(v)<<","<<std::get<2>(v)<<","<<std::get<3>(v)<<std::endl;
        }

        num_tx_per_sample_loc_ss.close();
    }

    std::map<std::string,std::tuple<int,int,int,int>> num_tx_per_sample_loc4;
    void save_sample_tx_stats4(std::string base_out_fname){
        std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc4";
        std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
        num_tx_per_sample_loc_ss<<"real,intronic,splicing,polymerase"<<std::endl;

        for(auto& v : num_tx_per_sample_loc4){
            num_tx_per_sample_loc_ss<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
        }

        num_tx_per_sample_loc_ss.close();
    }

    std::map<std::string,std::tuple<int,int,int,int>> num_tx_per_tissue_loc5;
    void save_tissue_tx_stats5(std::string base_out_fname){
        std::string num_tx_per_tissue_loc_fname = base_out_fname+".num_tx_per_tissue_loc5";
        std::ofstream num_tx_per_tissue_loc_ss(num_tx_per_tissue_loc_fname.c_str());
        num_tx_per_tissue_loc_ss<<"real,intronic,splicing,polymerase"<<std::endl;
        for(auto& v : num_tx_per_tissue_loc5){
            num_tx_per_tissue_loc_ss<<v.first<<","<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
        }
        num_tx_per_tissue_loc_ss.close();
    }

    std::vector<std::vector<std::pair<float,float>>> gauss_sample_tpm_per_tissue_loc;
    void save_sample_tpm_gauss(std::string base_out_fname){
        std::string gauss_sample_tpm_per_tissue_loc_fname = base_out_fname+".tpm_sample_gauss";
        std::ofstream gauss_sample_tpm_per_tissue_loc_ss(gauss_sample_tpm_per_tissue_loc_fname.c_str());
        gauss_sample_tpm_per_tissue_loc_ss<<"real_mean,real_sd,intronic_mean,intronic_sd,splicing_mean,splicing_sd,polymerase_mean,polymerase_sd";
        for(int i=0;i<gauss_sample_tpm_per_tissue_loc.size();i++){
            if(i%4==0){ // seen real,splice,int and pol and can output CR
                gauss_sample_tpm_per_tissue_loc_ss<<std::endl;
            }
            else{
                gauss_sample_tpm_per_tissue_loc_ss<<",";
            }
            for(auto& t : gauss_sample_tpm_per_tissue_loc[i]){
                gauss_sample_tpm_per_tissue_loc_ss<<t.first<<";";
            }
            gauss_sample_tpm_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
            gauss_sample_tpm_per_tissue_loc_ss << ",";
            for(auto& t : gauss_sample_tpm_per_tissue_loc[i]){
                gauss_sample_tpm_per_tissue_loc_ss<<t.second<<";";
            }
            gauss_sample_tpm_per_tissue_loc_ss.seekp(-1, std::ios_base::end);
        }
        gauss_sample_tpm_per_tissue_loc_ss<<std::endl;

        gauss_sample_tpm_per_tissue_loc_ss.close();
    }

    std::vector<std::pair<float,float>> gauss_sample_tx_per_tissue_loc;
    void save_sample_tx_gauss(std::string base_out_fname){
        std::string gauss_sample_tx_per_tissue_loc_fname = base_out_fname+".num_tx_sample_gauss";
        std::ofstream gauss_sample_tx_per_tissue_loc_ss(gauss_sample_tx_per_tissue_loc_fname.c_str());
        gauss_sample_tx_per_tissue_loc_ss<<"real_mean,real_sd,intronic_mean,intronic_sd,splicing_mean,splicing_sd,polymerase_mean,polymerase_sd";
        for(int i=0;i<gauss_sample_tx_per_tissue_loc.size();i++){
            if(i%4==0){ // seen real,splice,int and pol and can output CR
                gauss_sample_tx_per_tissue_loc_ss<<std::endl;
            }
            else{
                gauss_sample_tx_per_tissue_loc_ss<<",";
            }
            gauss_sample_tx_per_tissue_loc_ss<<gauss_sample_tx_per_tissue_loc[i].first<<","<<gauss_sample_tx_per_tissue_loc[i].second;
        }
        gauss_sample_tx_per_tissue_loc_ss<<std::endl;

        gauss_sample_tx_per_tissue_loc_ss.close();
    }

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
    //     <7> - float - tpms observed in polymerase transcripts
    std::map<std::string,std::tuple<std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>>> num_tx_per_sample_loc5;
    void save_sample_tx_stats5(std::string base_out_fname){
        std::string num_tx_per_sample_loc_fname = base_out_fname+".num_tx_per_sample_loc5";
        std::ofstream num_tx_per_sample_loc_ss(num_tx_per_sample_loc_fname.c_str());
        num_tx_per_sample_loc_ss<<"real,intronic,splicing,polymerase,tpms_real,tpms_intronic,tpms_splicing,tpms_polymerase"<<std::endl;

        for(auto& v : num_tx_per_sample_loc5){
            num_tx_per_sample_loc_ss<<std::get<0>(v.second).size()<<","<<std::get<1>(v.second).size()<<","<<std::get<2>(v.second).size()<<","<<std::get<3>(v.second).size()<<",";
            if(std::get<0>(v.second).size() == 0){
                num_tx_per_sample_loc_ss << "0;";
            }
            else{
                for(auto& s : std::get<0>(v.second)){
                    num_tx_per_sample_loc_ss << s <<";";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss << ",";
            if(std::get<1>(v.second).size() == 0){
                num_tx_per_sample_loc_ss << "0;";
            }
            else{
                for(auto& s : std::get<1>(v.second)){
                    num_tx_per_sample_loc_ss << s <<";";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss << ",";
            if(std::get<2>(v.second).size() == 0){
                num_tx_per_sample_loc_ss << "0;";
            }
            else{
                for(auto& s : std::get<2>(v.second)){
                    num_tx_per_sample_loc_ss << s <<";";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss << ",";
            if(std::get<3>(v.second).size() == 0){
                num_tx_per_sample_loc_ss << "0;";
            }
            else{
                for(auto& s : std::get<3>(v.second)){
                    num_tx_per_sample_loc_ss << s <<";";
                }
            }
            num_tx_per_sample_loc_ss.seekp(-1, std::ios_base::end);
            num_tx_per_sample_loc_ss << std::endl;
        }

        num_tx_per_sample_loc_ss.close();
    }

    int num_avg_real_locs_per_tissue=0,num_avg_noise_locs_per_tissue=0,num_avg_undef_locs_per_tissue=0;
    void save_loc_stats(std::string base_out_fname){
        std::string loc_fname = base_out_fname+".loc_stats";
        std::ofstream loc_ss(loc_fname.c_str());
        loc_ss<<"real: "<<num_avg_real_locs_per_tissue<<std::endl;
        loc_ss<<"noise: "<<num_avg_noise_locs_per_tissue<<std::endl;
        loc_ss<<"undefined: "<<num_avg_undef_locs_per_tissue<<std::endl;
        loc_ss.close();
    }

    std::map<std::string,std::tuple<int,int,int,int>> num_tx_per_tissue;
    void save_num_tx_tissue(std::string base_out_fname){
        std::string num_tx_per_tissue_fname = base_out_fname+".num_tx_per_tissue";
        std::ofstream num_tx_per_tissue_ss(num_tx_per_tissue_fname.c_str());
        num_tx_per_tissue_ss<<"real,splicing,intronic,polymerase"<<std::endl;
        for(auto& v : num_tx_per_tissue){
            num_tx_per_tissue_ss<<std::get<0>(v.second)<<","<<std::get<1>(v.second)<<","<<std::get<2>(v.second)<<","<<std::get<3>(v.second)<<std::endl;
        }
        num_tx_per_tissue_ss.close();
    }

    void save_stats(std::string base_out_fname){
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
//        save_tissue_loc(base_out_fname);
//        save_sample_loc(base_out_fname);
//        save_tissue_tx_stats3(base_out_fname);
//        save_sample_tx_stats3(base_out_fname);
//        save_sample_txs(base_out_fname);
//        save_sample_tx_stats4(base_out_fname);
//        save_sample_tx_stats5(base_out_fname);
//        save_tissue_tx_stats5(base_out_fname);
//        save_num_tx_tissue(base_out_fname);
//        save_sample_tx_gauss(base_out_fname);
        save_sample_tpm_gauss(base_out_fname);
        std::cout<<"done saving stats"<<std::endl;
    }
};

class TrackingTree{
public:
    TrackingTree(std::string ttf,std::string atf,std::string agf,std::string tr_tf,std::string sp_tf,std::string in_tf,std::string pl_tf):tissue_tracking_fname(ttf),all_tracking_fname(atf),all_gtf_fname(agf),true_gff_fname(tr_tf),
                                                                        splicing_gff_fname(sp_tf),intronic_gff_fname(in_tf),polymerase_gff_fname(pl_tf){}
    ~TrackingTree() = default;

    void load(){
        // load tissue_tracking
        std::cout<<"<<<loading tissue tracking info"<<std::endl;
        load_tt();
        // TODO: need to add some tracking to be able to compute the number of loci/transcripts/transcripts_per_locus of each type per sample/tissue
        // What if we output averages across samples?
        // How do we compute the average number of transcripts/loci per sample/tissue?

        // why not load it with the probabilities? for a given locus, what is the probability on a tissue level of observing this locus?
        // that's fine for separate modeling of loci and transcripts, but in order to model the relationship between each transcript and each locus there has to be some kind of conditional probability?
        // or we just write a gtf-like file of probabilities for each locus and exon for an average sample, and an average tissue, and then select things that way

        // anoter question then becomes, how do you model samples, such that they resemble an average tissue?
        // meaning what should sample-like probabilities look like? how are they conditioned on a sample?
        // same goes for coverages and expression

        // What we need
        //  1. probability of each locus being found in an average tissue (separatefor each type of locus)
        //  2. number of loci per average tissue
        // to do this we need
        //  1.

        // Do we really need that though?
        // Aren't we just going to use distributions to perform sampling?
        // we might still need an average number of loci and transcripts of each type, but otherwise, it should be possible to simply use distributions we have already computed

        // add all tracking
        std::cout<<"<<<loading all tracking info"<<std::endl;
        load_at();

        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading true transcript ids"<<std::endl;
        load_type(this->true_gff_fname,TYPE::REAL_TX);
        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading non-intronic/splicing noise transcript ids"<<std::endl;
        load_type(this->splicing_gff_fname,TYPE::SPLICING_TX);
        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading intronic noise transcript ids"<<std::endl;
        load_type(this->intronic_gff_fname,TYPE::INTRONIC_TX);
        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading poymerase noise transcript ids"<<std::endl;
        load_type(this->polymerase_gff_fname,TYPE::POLYMERASE_TX);
    }

    void get_num_tx_per_sample_locus5(Stats& stats){
        std::cout<<"computing the number of transcripts per locus per sample5"<<std::endl;
        std::pair<std::map<std::string,std::tuple<std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>>>::iterator,bool> stloc5_it;
        for(auto &atx : this->mattm){ // iterate over ALL transcripts
            int tx_type = atx.second.type;
            std::string tx_locus = atx.second.locus;
            for(auto& ttx : atx.second.txs){
                for(auto& stx : ttx.first->second.txs){
                    stloc5_it = stats.num_tx_per_sample_loc5.insert(std::make_pair(stx.first->second.sample+tx_locus,std::make_tuple(std::vector<float>{},std::vector<float>{},std::vector<float>{},std::vector<float>{})));
                    if(tx_type == TYPE::REAL_TX){
                        std::get<0>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                    }
                    else if(tx_type == TYPE::SPLICING_TX){
                        std::get<1>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                    }
                    else if(tx_type == TYPE::INTRONIC_TX){
                        std::get<2>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                    }
                    else if(tx_type == TYPE::POLYMERASE_TX){
                        std::get<3>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                    }
                    else{
                        continue;
                    }
                }
            }
        }
    }

    void get_num_tx_tissue(Stats& stats){
        // also need to write a number of transcripts per tissue to gauge whether the stage1 tissue simmulation is correct
        std::cout<<"computing the number of transcripts of each type in a tissue"<<std::endl;
        std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> tt_it;
        for(auto &atx : this->mattm){ // iterate over ALL transcripts
            int tx_type = atx.second.type;
            for(auto& ttx : atx.second.txs){
                tt_it = stats.num_tx_per_tissue.insert(std::make_pair(ttx.first->second.tissue,std::make_tuple(0,0,0,0)));
                if(tx_type == TYPE::REAL_TX){
                    std::get<0>(tt_it.first->second)++;
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    std::get<1>(tt_it.first->second)++;
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    std::get<2>(tt_it.first->second)++;
                }
                else if(tx_type == TYPE::POLYMERASE_TX){
                    std::get<3>(tt_it.first->second)++;
                }
                else{
                    continue;
                }
            }
        }
    }

    void get_num_tx_per_tissue_locus5(Stats& stats){
        std::cout<<"computing the number of transcripts of each type per tissue locus"<<std::endl;
        std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> ttloc5_it;
        for(auto &atx : this->mattm){ // iterate over ALL transcripts
            int tx_type = atx.second.type;
            std::string tx_locus = atx.second.locus;
            for(auto& ttx : atx.second.txs){
                ttloc5_it = stats.num_tx_per_tissue_loc5.insert(std::make_pair(ttx.first->second.tissue+tx_locus,std::make_tuple(0,0,0,0)));
                if(tx_type == TYPE::REAL_TX){
                    std::get<0>(ttloc5_it.first->second)++;
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    std::get<1>(ttloc5_it.first->second)++;
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    std::get<2>(ttloc5_it.first->second)++;
                }
                else if(tx_type == TYPE::POLYMERASE_TX){
                    std::get<3>(ttloc5_it.first->second)++;
                }
                else{
                    continue;
                }
            }
        }
    }

    void get_num_tx_per_sample_locus4(Stats& stats){
        std::cout<<"computing the number of transcripts per locus per sample4"<<std::endl;
        std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> stloc_it;
        for(auto &atx : this->mattm){ // iterate over ALL loci
            int tx_type = atx.second.type;
            std::string tx_locus = atx.second.locus;
            for(auto& ttx : atx.second.txs){
                for(auto& sample : ttx.first->second.samples){
                    stloc_it = stats.num_tx_per_sample_loc4.insert(std::make_pair(sample+tx_locus,std::make_tuple(0,0,0,0)));
                    if(tx_type == TYPE::REAL_TX){
                        std::get<0>(stloc_it.first->second) += 1;
                    }
                    else if(tx_type == TYPE::SPLICING_TX){
                        std::get<1>(stloc_it.first->second) += 1;
                    }
                    else if(tx_type == TYPE::INTRONIC_TX){
                        std::get<2>(stloc_it.first->second) += 1;
                    }
                    else if(tx_type == TYPE::POLYMERASE_TX){
                        std::get<3>(stloc_it.first->second) += 1;
                    }
                    else{
                        continue;
                    }
                }
            }
        }
    }

    void get_num_tx_sample(Stats& stats){
        // now need to do the same for sample level
        // namely we need to find how many ALL transcripts does a sample contibute to an ALL locus
        std::cout<<"computing the number of transcripts per sample"<<std::endl;
        std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> st_it;
        for(auto &atx : this->mattm){ // iterate over ALL loci
            int tx_type = atx.second.type;
            for(auto& ttx : atx.second.txs){
                for(auto& sample : ttx.first->second.samples){
                    st_it = stats.sample_txs.insert(std::make_pair(sample,std::make_tuple(0,0,0,0)));
                    if(tx_type == TYPE::REAL_TX){
                        std::get<0>(st_it.first->second) += 1;
                    }
                    else if(tx_type == TYPE::SPLICING_TX){
                        std::get<1>(st_it.first->second) += 1;
                    }
                    else if(tx_type == TYPE::INTRONIC_TX){
                        std::get<2>(st_it.first->second) += 1;
                    }
                    else if(tx_type == TYPE::POLYMERASE_TX){
                        std::get<3>(st_it.first->second) += 1;
                    }
                    else{
                        continue;
                    }
                }
            }
        }
    }

    void get_num_tx_per_tissue_locus3(Stats& stats){
        // Now can we also compute the number of loci per tissue and per sample
        //    as well as the number of transcripts per tissue and per sample?
        // 1. how many true loci?
        // 2. how many false loci?
        // 3. how many of the true loci have noise transcripts

        // need to get averages across tissues: how?
        std::cout<<"computing the number of transcripts per locus per tissue3"<<std::endl;
        std::map<std::string,int> num_avg_real_locs,num_avg_noise_locs,num_avg_undef_locs;
        std::pair<std::map<std::string,int>::iterator,bool> num_it;
        std::cout<<"number of loci: "<<this->loci.size()<<std::endl;
        // the question that we should be asking instead is how many ALL transcripts does a given tissue contain within an ALL locus (instead of _loc and _loc2 files)
        for(auto& loc : this->loci){ // begin iterating over all loci
            // we only care about the composition of loci for tissues where a locus is expressed
            std::map<std::string,std::tuple<int,int,int,int>> found_tissues; // tuple contains real,nonint,int,pol
            std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> ft_it;
            for(auto& tx : loc.second.first){
                for(auto& tissue : tx.first->second.tissues){
                    ft_it = found_tissues.insert(std::make_pair(tissue,std::make_tuple(0,0,0,0)));
                    if(tx.first->second.type == TYPE::REAL_TX){
                        std::get<0>(ft_it.first->second) += 1;
                    }
                    else if(tx.first->second.type == TYPE::SPLICING_TX){
                        std::get<1>(ft_it.first->second) += 1;
                    }
                    else if(tx.first->second.type == TYPE::INTRONIC_TX){
                        std::get<2>(ft_it.first->second) += 1;
                    }
                    else if(tx.first->second.type == TYPE::POLYMERASE_TX){
                        std::get<3>(ft_it.first->second) += 1;
                    }
                    else{
                        continue;
                    }
                }
            }
            for(auto& tissue : found_tissues){
                stats.num_tx_per_tissue_loc3.push_back(tissue.second);
            }
        }
    }

    void get_real_noise_locs_sample(Stats& stats){
        std::cout<<"computing the number of real and noise ALL loci as identified in each samples"<<std::endl;
        std::unordered_map<std::string,int> sample_locs_real,sample_locs_noise,sample_locs_other;
        std::pair<std::unordered_map<std::string,int>::iterator,bool> sit;
        for(auto& loc : this->loci){ // for now let's just output everything
            int loc_type = loc.second.second;
            std::set<std::string> all_samples;
            for(auto& tx : loc.second.first){
                for(auto& sx : tx.first->second.txs){
                    all_samples.insert(sx.first->second.samples.begin(),sx.first->second.samples.end());
                }
            }
            for(auto& sample : all_samples){
                if(loc_type==1){
                    sit = sample_locs_real.insert(std::make_pair(sample,0));
                    sit.first->second++;
                }
                else if(loc_type==0){ // is 0
                    sit = sample_locs_noise.insert(std::make_pair(sample,0));
                    sit.first->second++;
                }
                else{ // is undefined
                    continue;
                }
            }
        }
        for(auto& tt : sample_locs_real){
            stats.num_real_locs_sample.push_back(tt.second);
        }
        for(auto& tt : sample_locs_noise){
            stats.num_noise_locs_sample.push_back(tt.second);
        }
    }

    void get_real_noise_locs_tissue(Stats& stats){
        // compute number of loci per sample - how?
        // do we look at the ones in main merge and trace it back?
        // how do we tell the type of locus?
        // if it contains any real transcripts - it is real


        // what should we output?
        // 1. number of noise (only RNApol transcripts in locus)
        // 2. real loci - anything that contains nonint/int/real transcripts and not RNApol
        // 2. other - whatever does not fit these definitions


        // --------
        // alternatively follow the same solution as below
        // but instead of getting the total average - need to get a distribution instead
        std::cout<<"computing the number of real and noise ALL loci per identified in each tissue"<<std::endl;
        std::unordered_map<std::string,int> tissue_locs_real,tissue_locs_noise,tissue_locs_other;
        std::pair<std::unordered_map<std::string,int>::iterator,bool> tit;
        for(auto& loc : this->loci){ // for now let's just output everything
            int loc_type = loc.second.second;
            std::set<std::string> all_tissues;
            for(auto& tx : loc.second.first){
                all_tissues.insert(tx.first->second.tissues.begin(),tx.first->second.tissues.end());
            }
            for(auto& tissue : all_tissues){
                if(loc_type==1){
                    tit = tissue_locs_real.insert(std::make_pair(tissue,0));
                    tit.first->second++;
                }
                else if(loc_type==0){ // is 0
                    tit = tissue_locs_noise.insert(std::make_pair(tissue,0));
                    tit.first->second++;
                }
                else{ // is undefined
                    continue;
                }
            }
        }
        for(auto& tt : tissue_locs_real){
            stats.num_real_locs_tissue.push_back(tt.second);
        }
        for(auto& tt : tissue_locs_noise){
            stats.num_noise_locs_tissue.push_back(tt.second);
        }
    }

    void get_cov_sample(Stats& stats){
        // Now to figure out how to report coverages and the rest
        std::cout<<"computing sample transcript coverages"<<std::endl;
        for(auto& sit : this->mtm){
            switch(sit.second.type){
                case TYPE::REAL_TX:
                    stats.cov_sample_real.push_back(sit.second.cov);
                    stats.fpkm_sample_real.push_back(sit.second.fpkm);
                    stats.tpm_sample_real.push_back(sit.second.tpm);
                    break;
                case TYPE::INTRONIC_TX:
                    stats.cov_sample_intronic.push_back(sit.second.cov);
                    stats.fpkm_sample_intronic.push_back(sit.second.fpkm);
                    stats.tpm_sample_intronic.push_back(sit.second.tpm);
                    break;
                case TYPE::SPLICING_TX:
                    stats.cov_sample_splicing.push_back(sit.second.cov);
                    stats.fpkm_sample_splicing.push_back(sit.second.fpkm);
                    stats.tpm_sample_splicing.push_back(sit.second.tpm);
                    break;
                case TYPE::POLYMERASE_TX:
                    stats.cov_sample_polymerase.push_back(sit.second.cov);
                    stats.fpkm_sample_polymerase.push_back(sit.second.fpkm);
                    stats.tpm_sample_polymerase.push_back(sit.second.tpm);
                    break;
                default:
                    stats.cov_sample_left.push_back(sit.second.cov);
                    stats.fpkm_sample_left.push_back(sit.second.fpkm);
                    stats.tpm_sample_left.push_back(sit.second.tpm);
                    break;
            }
        }
    }

    void get_tpm_fracs_joined(Stats& stats){
        // same for TPMs
        std::cout<<"computing the same for the TPMs"<<std::endl;
        for(auto& v : this->loci_sample) {
            std::vector<float> sums_reals,sums_splicings,sums_intronics;
            float sum_real = 0,sum_splicing = 0,sum_intronic = 0;
            int num_real = 0,num_splicing = 0,num_intronic = 0;
            for (auto &sit : v.second) {
                switch(sit.first->second.type){
                    case TYPE::REAL_TX:
                        num_real += 1;
                        sum_real += sit.first->second.tpm;
                        sums_reals.push_back(sit.first->second.tpm);
                        break;
                    case TYPE::SPLICING_TX:
                        num_splicing += 1;
                        sum_splicing += sit.first->second.tpm;
                        sums_splicings.push_back(sit.first->second.tpm);
                        break;
                    case TYPE::INTRONIC_TX:
                        num_intronic += 1;
                        sum_intronic += sit.first->second.tpm;
                        sums_intronics.push_back(sit.first->second.tpm);
                        break;
                    default:
                        break;
                }
            }
            if(num_real>0){ // only do for the real transcripts
                stats.tpm_fracs.push_back(std::tuple<int,float,int,float,int,float>(num_real,sum_real,num_splicing,sum_splicing,num_intronic,sum_intronic));
                stats.tpms_joined.push_back(std::tuple<std::vector<float>,std::vector<float>,std::vector<float>>(sums_reals,sums_splicings,sums_intronics));
            }
        }
    }

    void get_cov_fracs_joined(Stats& stats){
        // what we can do is:
        // 1. for each locus save the mean real, mean intronic, mean splicing, mean polymerase coverages/tpms (sample level)
        // 2. when simulating - select the fraction from the tissue
        // 3. how will this work when the number of transcripts varries?

        // also, probably shouldn't be using coverages, since they are normalized across samples with varying sequencing depth... - TPMs would be better, but then there is no way of getting a coverage/nreads out of it...

        // also we should probably begin by computing distributions based on the simulation parameters, so that we know if they look similar to those observed in GTEx

        // however, something like this needs to be done on the tissue level and the propagated onto the sample level - question is how??? - does it though?
        //   on tissue level, we only get the number of samples - no need for abundance - that can be done at sample level

        // question is how do we compute the fraction information. What we need:
        // 1. TPM
        //    - real
        //    - nonint
        // 2. Number of Transcripts
        //    - real
        //    - nonint
        std::cout<<"computing coverage fracs and coverage joined"<<std::endl;
        for(auto& v : this->loci_sample) {
            std::vector<float> sums_reals,sums_splicings,sums_intronics;
            float sum_real = 0,sum_splicing = 0,sum_intronic = 0;
            int num_real = 0,num_splicing = 0,num_intronic = 0;
            for (auto &sit : v.second) {
                switch(sit.first->second.type){
                    case TYPE::REAL_TX:
                        num_real += 1;
                        sum_real += sit.first->second.cov;
                        sums_reals.push_back(sit.first->second.cov);
                        break;
                    case TYPE::SPLICING_TX:
                        num_splicing += 1;
                        sum_splicing += sit.first->second.cov;
                        sums_splicings.push_back(sit.first->second.cov);
                        break;
                    case TYPE::INTRONIC_TX:
                        num_intronic += 1;
                        sum_intronic += sit.first->second.cov;
                        sums_intronics.push_back(sit.first->second.cov);
                        break;
                    default:
                        break;
                }
            }
            if(num_real>0){ // only do for the real transcripts
                stats.cov_fracs.push_back(std::tuple<int,float,int,float,int,float>(num_real,sum_real,num_splicing,sum_splicing,num_intronic,sum_intronic));
                stats.covs_joined.push_back(std::tuple<std::vector<float>,std::vector<float>,std::vector<float>>(sums_reals,sums_splicings,sums_intronics));
            }
        }
    }

    void get_num_tx_per_sample_locus3(Stats& stats){
        std::cout<<"computing the number of transcripts per locus per sample3"<<std::endl;
        for(auto& loc : this->loci){ // begin iterating over all loci
            // we only care about the composition of loci for samples where a locus is expressed
            std::map<std::string,std::tuple<int,int,int,int>> found_samples; // tuple contains real,nonint,int,pol
            std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> ft_it;
            for(auto& tx : loc.second.first){
                for(auto& stx : tx.first->second.txs){ // now to descend onto the sample level
                    for(auto& sample : stx.first->second.samples){
                        ft_it = found_samples.insert(std::make_pair(sample,std::make_tuple(0,0,0,0)));
                        if(tx.first->second.type == TYPE::REAL_TX){
                            std::get<0>(ft_it.first->second) += 1;
                        }
                        else if(tx.first->second.type == TYPE::SPLICING_TX){
                            std::get<1>(ft_it.first->second) += 1;
                        }
                        else if(tx.first->second.type == TYPE::INTRONIC_TX){
                            std::get<2>(ft_it.first->second) += 1;
                        }
                        else if(tx.first->second.type == TYPE::POLYMERASE_TX){
                            std::get<3>(ft_it.first->second) += 1;
                        }
                        else{
                            continue;
                        }
                    }
                }
            }
            for(auto& tissue : found_samples){
                stats.num_tx_per_sample_loc3.push_back(tissue.second);
            }
        }
    }

    void get_num_tx_per_sample_locus2(Stats& stats){
        // now to figure out how to report the average number of transcripts which a sample contributes to an tissue locus
        std::cout<<"computing the number of transcripts per sample loc2"<<std::endl;
        for(auto& lit : this->loci_tissue){
            int num_txs_total=lit.second.size();
            std::map<std::string,int> num_txs_real,num_txs_intronic,num_txs_splicing,num_txs_polymerase,num_txs_left;
            std::pair<std::map<std::string,int>::iterator,bool> nit;
            for(auto& tit : lit.second){
                // each tissue should only appear once in her - make sure that is so
                for(auto& ts : tit.first->second.samples){
                    switch(tit.first->second.type){
                        case TYPE::REAL_TX:
                            nit = num_txs_real.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        case TYPE::INTRONIC_TX:
                            nit = num_txs_intronic.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        case TYPE::SPLICING_TX:
                            nit = num_txs_splicing.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        case TYPE::POLYMERASE_TX:
                            nit = num_txs_polymerase.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        default:
                            nit = num_txs_left.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                    }
                }
            }
            // compute averages across tissues now
            int sum_real=0,sum_int=0,sum_nonint=0,sum_pol=0,sum_left=0;
            for(auto& nt : num_txs_real){sum_real+=nt.second;}
            int nt_real = 0;
            if(sum_real>0){
                nt_real = (int)((float)sum_real/(float)num_txs_real.size());
            }

            for(auto& nt : num_txs_intronic){sum_int+=nt.second;}
            int nt_int = 0;
            if(sum_int>0){
                nt_int = (int)((float)sum_int/(float)num_txs_intronic.size());
            }

            for(auto& nt : num_txs_splicing){sum_nonint+=nt.second;}
            int nt_nonint = 0;
            if(sum_nonint>0){
                nt_nonint = (int)((float)sum_nonint/(float)num_txs_splicing.size());
            }

            for(auto& nt : num_txs_polymerase){sum_pol+=nt.second;}
            int nt_pol = 0;
            if(sum_pol>0){
                nt_pol = (int)((float)sum_pol/(float)num_txs_polymerase.size());
            }

            for(auto& nt : num_txs_left){sum_left+=nt.second;}
            int nt_left = 0;
            if(sum_left>0){
                nt_left = (int)((float)sum_left/(float)num_txs_left.size());
            }

            // save the average number of transcripts per ALL locus for each type of transcripts
            stats.num_tx_per_sample_loc_real2.push_back(nt_real);
            stats.num_tx_per_sample_loc_intronic2.push_back(nt_int);
            stats.num_tx_per_sample_loc_splicing2.push_back(nt_nonint);
            stats.num_tx_per_sample_loc_polymerase2.push_back(nt_pol);
            stats.num_tx_per_sample_loc_left2.push_back(nt_left);
            stats.num_tx_per_sample_loc2.push_back(num_txs_total);
        }
    }

    void get_num_tx_per_tissue_locus2(Stats& stats){
        // now to figure out how to report the average number of transcripts which a tissue contributes to an ALL locus
        std::cout<<"computing the number of transcripts per tissue loc2"<<std::endl;
        for(auto& lit : this->loci){
            int num_txs_total=lit.second.first.size();
            std::map<std::string,int> num_txs_real,num_txs_intronic,num_txs_splicing,num_txs_polymerase,num_txs_left;
            std::pair<std::map<std::string,int>::iterator,bool> nit;
            for(auto& tit : lit.second.first){
                // each tissue should only appear once in her - make sure that is so
                for(auto& ts : tit.first->second.tissues){
                    switch(tit.first->second.type){
                        case TYPE::REAL_TX:
                            nit = num_txs_real.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        case TYPE::INTRONIC_TX:
                            nit = num_txs_intronic.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        case TYPE::SPLICING_TX:
                            nit = num_txs_splicing.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        case TYPE::POLYMERASE_TX:
                            nit = num_txs_polymerase.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                        default:
                            nit = num_txs_left.insert(std::make_pair(ts,0));
                            nit.first->second++;
                            break;
                    }
                }
            }
            // compute averages across tissues now
            int sum_real=0,sum_int=0,sum_nonint=0,sum_pol=0,sum_left=0;
            for(auto& nt : num_txs_real){sum_real+=nt.second;}
            int nt_real = 0;
            if(sum_real>0){
                nt_real = (int)((float)sum_real/(float)num_txs_real.size());
            }

            for(auto& nt : num_txs_intronic){sum_int+=nt.second;}
            int nt_int = 0;
            if(sum_int>0){
                nt_int = (int)((float)sum_int/(float)num_txs_intronic.size());
            }

            for(auto& nt : num_txs_splicing){sum_nonint+=nt.second;}
            int nt_nonint = 0;
            if(sum_nonint>0){
                nt_nonint = (int)((float)sum_nonint/(float)num_txs_splicing.size());
            }

            for(auto& nt : num_txs_polymerase){sum_pol+=nt.second;}
            int nt_pol = 0;
            if(sum_pol>0){
                nt_pol = (int)((float)sum_pol/(float)num_txs_polymerase.size());
            }

            for(auto& nt : num_txs_left){sum_left+=nt.second;}
            int nt_left = 0;
            if(sum_left>0){
                nt_left = (int)((float)sum_left/(float)num_txs_left.size());
            }

            // save the average number of transcripts per ALL locus for each type of transcripts
            stats.num_tx_per_tissue_loc_real2.push_back(nt_real);
            stats.num_tx_per_tissue_loc_intronic2.push_back(nt_int);
            stats.num_tx_per_tissue_loc_splicing2.push_back(nt_nonint);
            stats.num_tx_per_tissue_loc_polymerase2.push_back(nt_pol);
            stats.num_tx_per_tissue_loc_left2.push_back(nt_left);
            stats.num_tx_per_tissue_loc2.push_back(num_txs_total);
        }
    }

    void get_num_tx_per_sample_locus(Stats& stats){
        // Now do the same on a sample level
        std::cout<<"computing the number of transcripts per sample locus"<<std::endl;
        for(auto& v : this->loci_sample){
            int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_splicing=0,num_txs_polymerase=0,num_txs_left=0;
            float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_splicing=0,sum_txs_polymerase=0,sum_txs_left=0;
            for(auto& sit : v.second){
                num_txs_all += 1;
                switch(sit.first->second.type){
                    case TYPE::REAL_TX:
                        num_txs_real += 1;
                        sum_txs_real += sit.first->second.tpm;
                        break;
                    case TYPE::INTRONIC_TX:
                        num_txs_intronic += 1;
                        sum_txs_intronic += sit.first->second.tpm;
                        break;
                    case TYPE::SPLICING_TX:
                        num_txs_splicing += 1;
                        sum_txs_splicing += sit.first->second.tpm;
                        break;
                    case TYPE::POLYMERASE_TX:
                        num_txs_polymerase += 1;
                        sum_txs_polymerase += sit.first->second.tpm;
                        break;
                    default:
                        num_txs_left += 1;
                        sum_txs_left += sit.first->second.tpm;
                        break;
                }
            }
            stats.num_tx_per_sample_loc_real.push_back(num_txs_real);
            stats.num_tx_per_sample_loc_intronic.push_back(num_txs_intronic);
            stats.num_tx_per_sample_loc_splicing.push_back(num_txs_splicing);
            stats.num_tx_per_sample_loc_polymerase.push_back(num_txs_polymerase);
            stats.num_tx_per_sample_loc_left.push_back(num_txs_left);
            stats.num_tx_per_sample_loc.push_back(num_txs_all);

            stats.sum_tx_per_sample_loc_real.push_back(sum_txs_real);
            stats.sum_tx_per_sample_loc_intronic.push_back(sum_txs_intronic);
            stats.sum_tx_per_sample_loc_splicing.push_back(sum_txs_splicing);
            stats.sum_tx_per_sample_loc_polymerase.push_back(sum_txs_polymerase);
            stats.sum_tx_per_sample_loc_left.push_back(sum_txs_left);
            stats.sum_tx_per_sample_loc.push_back(sum_txs_all);
        }
    }

    void get_num_tx_per_tissue_locus(Stats& stats){
        // Now do the same on a tissue level
        std::cout<<"computing the number of transcripts per tissue locus"<<std::endl;
        for(auto& v : this->loci_tissue){
            int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_splicing=0,num_txs_polymerase=0,num_txs_left=0;
            float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_splicing=0,sum_txs_polymerase=0,sum_txs_left=0;
            for(auto& sit : v.second){
                num_txs_all += sit.first->second.get_num_txs();
                switch(sit.first->second.type){
                    case TYPE::REAL_TX:
                        num_txs_real += sit.first->second.get_num_txs();
                        sum_txs_real += sit.first->second.get_sum_tpms();
                        break;
                    case TYPE::INTRONIC_TX:
                        num_txs_intronic += sit.first->second.get_num_txs();
                        sum_txs_intronic += sit.first->second.get_sum_tpms();
                        break;
                    case TYPE::SPLICING_TX:
                        num_txs_splicing += sit.first->second.get_num_txs();
                        sum_txs_splicing += sit.first->second.get_sum_tpms();
                        break;
                    case TYPE::POLYMERASE_TX:
                        num_txs_polymerase += sit.first->second.get_num_txs();
                        sum_txs_polymerase += sit.first->second.get_sum_tpms();
                        break;
                    default:
                        num_txs_left += sit.first->second.get_num_txs();
                        sum_txs_left += sit.first->second.get_sum_tpms();
                        break;
                }
            }
            stats.num_tx_per_tissue_loc_real.push_back(num_txs_real);
            stats.num_tx_per_tissue_loc_intronic.push_back(num_txs_intronic);
            stats.num_tx_per_tissue_loc_splicing.push_back(num_txs_splicing);
            stats.num_tx_per_tissue_loc_polymerase.push_back(num_txs_polymerase);
            stats.num_tx_per_tissue_loc_left.push_back(num_txs_left);
            stats.num_tx_per_tissue_loc.push_back(num_txs_all);

            stats.sum_tx_per_tissue_loc_real.push_back(sum_txs_real);
            stats.sum_tx_per_tissue_loc_intronic.push_back(sum_txs_intronic);
            stats.sum_tx_per_tissue_loc_splicing.push_back(sum_txs_splicing);
            stats.sum_tx_per_tissue_loc_polymerase.push_back(sum_txs_polymerase);
            stats.sum_tx_per_tissue_loc_left.push_back(sum_txs_left);
            stats.sum_tx_per_tissue_loc.push_back(sum_txs_all);
        }
    }

    void get_loc_stats(Stats& stats){
        std::map<std::string,int> num_avg_real_locs,num_avg_noise_locs,num_avg_undef_locs;
        std::pair<std::map<std::string,int>::iterator,bool> num_it;
        std::cout<<"computing average locus stats"<<std::endl;
        for(auto& loc : this->loci){
            // make sure all come from the same tissue!
            int loc_type = loc.second.second;
            std::set<std::string> all_tissues;
            for(auto& tx : loc.second.first){
                all_tissues.insert(tx.first->second.tissues.begin(),tx.first->second.tissues.end());
            }
            for(auto& tissue : all_tissues){
                if(loc_type==1){
                    num_it = num_avg_real_locs.insert(std::make_pair(tissue,0));
                    num_it.first->second++;
                }
                else if(loc_type==0){ // is 0
                    num_it = num_avg_noise_locs.insert(std::make_pair(tissue,0));
                    num_it.first->second++;
                }
                else{ // is undefined
                    num_it = num_avg_undef_locs.insert(std::make_pair(tissue,0));
                    num_it.first->second++;
                }
            }
        }
        // get average numbers here
        int sum_num_locs_real=0,sum_num_locs_noise=0,sum_num_locs_undef=0;
        for(auto& t : num_avg_noise_locs){
            sum_num_locs_noise+=t.second;
        }
        for(auto& t : num_avg_real_locs){
            sum_num_locs_real+=t.second;
        }
        for(auto& t : num_avg_undef_locs){
            sum_num_locs_undef+=t.second;
        }
        stats.num_avg_noise_locs_per_tissue = (int)((float)sum_num_locs_noise/(float)num_avg_noise_locs.size());
        stats.num_avg_real_locs_per_tissue = (int)((float)sum_num_locs_real/(float)num_avg_real_locs.size());
        stats.num_avg_undef_locs_per_tissue = (int)((float)sum_num_locs_undef/(float)num_avg_undef_locs.size());
    }

    void get_gauss_sample_tpm_per_tissue_loc(Stats& stats){
        // need a function which computes the parameters of a normal distribution for the number of transcripts per locus in a tissue
        // and for the TPMs per locus in a tissue
        // we can use this to do the following:
        // 1. selecting total number of transcripts in a tissue (and selecting them)
        // 2. selecting the TPMs and the number of transcripts for each sample
        // The problem is that we will invaitable loose information about the relationship etwee real/splicing/intronic transcripts per single locus
        //  Is there a way to construct the parameters for the joint normal distribution, which accounts for all three?

        // can we turn this number into a multivariate distribution instead?
        // for each locus, we want to know mean and covariance matrix between the 3 variants (real,splicing,intronic)
        // however, we want to do this on a tissue level
        // meaning that for a given tissue and a given ALL locus, what is the mean expression of each group and the covariance between them?

        // for polymerase we can simply compute the mean and std

        // to do this we need to store the following
        // the problem here is that we have a different number of transcripts for each
        // one way to circumvent this would be to build separate 1-D gausians for each type

        // NEW PLAN
        // create a distribution of means and stds to assign to each transcript in a tissue
        // when creating samples, each tissue will be broken down into samples
        // where each sample receives a number of each type of transcript
        // the number which each receives can also be computed by the same distribution with respect to number of transcripts instead of TPMs
        // and then when splitting into samples, we can simply sample from those distributions

        // also initialize for the TPMS
        // however here we might need to do on a transcript level - each transcript stored and reported individually
//        typedef std::map<std::string,int> ST; // holds a map from transcript ID to all expressions of that transcript
//        std::pair<std::map<std::string,int>::iterator,bool> st_it;
        typedef std::map<std::string,std::vector<float>> STPM; // holds a map of sample to all transcripts
        std::pair<STPM::iterator,bool> stpm_it;
        std::map<std::string,std::array<STPM,4>> tpm_txs;
        std::pair<std::map<std::string,std::array<STPM,4>>::iterator,bool> tp_it;

        for(auto &atx : this->mattm){ // iterate over ALL transcripts
            int tx_type = atx.second.type;
            std::string tx_locus = atx.second.locus;
            for(auto& ttx : atx.second.txs){
                for(auto& stx : ttx.first->second.txs){
                    tp_it = tpm_txs.insert(std::make_pair(stx.first->second.tissue+tx_locus,std::array<STPM,4>{STPM{},STPM{},STPM{},STPM{}}));
                    if(tx_type == TYPE::REAL_TX){
                        stpm_it = tp_it.first->second[0].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                        stpm_it.first->second.push_back(stx.first->second.tpm);
                    }
                    else if(tx_type == TYPE::SPLICING_TX){
                        stpm_it = tp_it.first->second[1].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                        stpm_it.first->second.push_back(stx.first->second.tpm);
                    }
                    else if(tx_type == TYPE::INTRONIC_TX){
                        stpm_it = tp_it.first->second[2].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                        stpm_it.first->second.push_back(stx.first->second.tpm);
                    }
                    else if(tx_type == TYPE::POLYMERASE_TX){
                        stpm_it = tp_it.first->second[3].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                        stpm_it.first->second.push_back(stx.first->second.tpm);
                    }
                    else{
                        continue;
                    }
                }
            }
        }

        // now we need to process these to get std and mean
        for(auto& v : tpm_txs){ // iterate over loci at the tissue level
            for(int i=0;i<4;i++){ // iterate over the real,splice,intr,pol
                stats.gauss_sample_tpm_per_tissue_loc.push_back(std::vector<std::pair<float,float>>{});
                if(v.second[i].size()==0){ // no transcripts of the given type found for the tissue ALL locus
                    stats.gauss_sample_tpm_per_tissue_loc.back().push_back(std::make_pair(0,0));
                    continue;
                }
                for(auto& a : v.second[i]){ // iterate over the transcripts
                    float total_num = 0;
                    float sum = 0;
                    for(auto& t : a.second){
                        sum+=t;
                        total_num++;
                    }
                    float mean = 0;
                    if(total_num>0){
                        mean = sum/total_num;
                    }
                    // compute SD
                    float total = 0;
                    for(auto& t : a.second){
                        total += (t - mean)*(t - mean);
                    }
                    float sd = sqrt(total / total_num);
                    stats.gauss_sample_tpm_per_tissue_loc.back().push_back(std::make_pair(mean,sd));
                }
            }
        }
    }

    void get_gauss_sample_tx_per_tissue_loc(Stats& stats){
        // what we need for this is a mean and a std for bth the number of transcripts and the corresponding tpm generators
        std::cout<<"computing gaussians of the number of transcripts on tissue level"<<std::endl;
        typedef std::map<std::string,int> STX;
        std::pair<STX::iterator,bool> stx_it;
        std::map<std::string,std::array<STX,4>> num_txs;
        std::pair<std::map<std::string,std::array<STX,4>>::iterator,bool> nt_it;

        for(auto &atx : this->mattm){ // iterate over ALL transcripts
            int tx_type = atx.second.type;
            std::string tx_locus = atx.second.locus;
            for(auto& ttx : atx.second.txs){
                for(auto& stx : ttx.first->second.txs){
                    nt_it = num_txs.insert(std::make_pair(stx.first->second.tissue+tx_locus,std::array<STX,4>{STX{},STX{},STX{},STX{}}));
                    if(tx_type == TYPE::REAL_TX){
                        stx_it = nt_it.first->second[0].insert(std::make_pair(stx.first->second.sample,0));
                        stx_it.first->second++;
                    }
                    else if(tx_type == TYPE::SPLICING_TX){
                        stx_it = nt_it.first->second[1].insert(std::make_pair(stx.first->second.sample,0));
                        stx_it.first->second++;
                    }
                    else if(tx_type == TYPE::INTRONIC_TX){
                        stx_it = nt_it.first->second[2].insert(std::make_pair(stx.first->second.sample,0));
                        stx_it.first->second++;
                    }
                    else if(tx_type == TYPE::POLYMERASE_TX){
                        stx_it = nt_it.first->second[3].insert(std::make_pair(stx.first->second.sample,0));
                        stx_it.first->second++;
                    }
                    else{
                        continue;
                    }
                }
            }
        }

        // now we need to process these to get std and mean
        for(auto& v : num_txs){ // iterate over loci at the tissue level
            for(int i=0;i<4;i++){ // iterate over the real,splice,intr,pol
                float total_num = 0;
                float sum = 0;
                for(auto& a : v.second[i]){ // iterate over the samples -> for each sample map inside STX
                    sum+=a.second;
                    total_num++;
                }
                if(total_num==0){ // not present
                    stats.gauss_sample_tx_per_tissue_loc.push_back(std::make_pair(0,0));
                    continue;
                }
                float mean = sum/total_num;
                float total = 0;
                for(auto& a : v.second[i]){
                    total += ((float)a.second - mean)*((float)a.second - mean);
                }
                float sd = sqrt(total / total_num);
                stats.gauss_sample_tx_per_tissue_loc.push_back(std::make_pair(mean,sd));
            }
        }
    }

    void get_stats(Stats& stats){
        std::cout<<"<<<computing stats"<<std::endl;

//        get_num_tx_per_sample_locus(stats);
//        get_num_tx_per_tissue_locus(stats);
//
//        get_num_tx_per_tissue_locus2(stats);
//        get_num_tx_per_sample_locus2(stats);
//
//        get_num_tx_per_tissue_locus3(stats);
//        get_num_tx_per_sample_locus3(stats);
//
//        get_cov_fracs_joined(stats);
//        get_tpm_fracs_joined(stats);
//
//        get_cov_sample(stats);
//
//        get_real_noise_locs_tissue(stats);
//        get_real_noise_locs_sample(stats);
//
//        get_num_tx_sample(stats);
//        get_num_tx_tissue(stats);
//
//        get_num_tx_per_sample_locus4(stats);
//
//        get_num_tx_per_tissue_locus5(stats);
//        get_num_tx_per_sample_locus5(stats);
//
//        get_loc_stats(stats);

//        get_gauss_sample_tx_per_tissue_loc(stats);
        get_gauss_sample_tpm_per_tissue_loc(stats);
    }

private:
    MTM mtm;
    MTM_IT mtm_it;
    MTTM mttm;
    MTTM_IT mttm_it;
    MATTM mattm;
    MATTM_IT mattm_it;
    Loci loci;
    Loci_IT loci_it;
    Tissues tissues;
    Tissues_IT tissues_it;
    Loci_Tissue loci_tissue;
    Loci_Tissue_IT loci_tissue_it;
    Loci_Sample loci_sample;
    Loci_Sample_IT loci_sample_it;

    int loc_is_real(std::vector<MATTM_IT>& txs){
        bool found_real = false;
        bool found_nonint = false;
        bool found_int = false;
        bool found_pol = false;
        // TODO: should we instead add intronic noise to the real? since such loci will automatical
        for(auto& tx : txs){
            if(tx.first->second.type == TYPE::REAL_TX){
                found_real = true;
            }
            if(tx.first->second.type == TYPE::SPLICING_TX){
                found_nonint = true;
            }
            if(tx.first->second.type == TYPE::INTRONIC_TX){
                found_int = true;
            }
            if(tx.first->second.type == TYPE::POLYMERASE_TX){
                found_pol = true;
            }
        }
        if(!found_real && !found_nonint && !found_int && !found_pol){
            return TYPE::UNDEFINED_TX;
        }
        if((found_real || found_nonint) && (!found_int && !found_pol)){
            return 1;
        }
        else if((found_int || found_pol) && (!found_real && !found_nonint)){
            return 0;
        }
        else {
            return 2;
        }
    }

    bool loc_tissue_is_real(std::vector<MTTM_IT>& txs){
        for(auto& tx : txs){
            if(tx.first->second.type != TYPE::REAL_TX){
                return true;
            }
        }
        return false;
    }

    bool loc_sample_is_real(std::vector<MTM_IT>& txs){
        for(auto& tx : txs){
            if(tx.first->second.type != TYPE::REAL_TX){
                return true;
            }
        }
        return false;
    }

    std::string get_tissue_name(std::string fname){
        std::size_t found_sample = fname.rfind("/");
        if (found_sample!=std::string::npos){

            std::size_t found_ext = fname.rfind(".");
            if (found_ext!=std::string::npos){
                return fname.substr(found_sample+1,found_ext-found_sample-1);
            }
            else{
                std::cerr<<"no extension separator found in the path: "<<fname<<std::endl;
                exit(-1);
            }
        }
        else{
            std::cerr<<"no sample separator found in the path: "<<fname<<std::endl;
            exit(-1);
        }
    }

    void break_tracking(std::vector<TR>& trs,std::string& tline){
        std::string new_tid,new_loc;
        // get tid by finding first tab
        std::size_t found_tid = tline.find("\t");
        if (found_tid!=std::string::npos){
            new_tid = tline.substr(0,found_tid);
        }
        else{
            std::cerr<<"ERROR #1"<<std::endl;
            exit(-1);
        }
        // get locus
        std::size_t found_gid = tline.find("\t",found_tid+1);
        if(found_gid!=std::string::npos){
            new_loc = tline.substr(found_tid+1,found_gid-found_tid - 1);
        }
        // get individual transcripts
        std::size_t found_name = tline.find("\t",found_gid+1);
        std::size_t found_code = tline.find("\t",found_name+1);

        std::stringstream *tx_stream = new std::stringstream(tline.substr(found_code+1,tline.size()-found_code));
        std::string tx;
        while(std::getline(*tx_stream,tx,'\t')) {
            if (std::strcmp(tx.c_str(), "-") == 0) { // transcript does not exist in a sample
                continue;
            }
            TR tr;
            tr.new_loc = new_loc;
            tr.new_tid = new_tid;
            // get qJ
            std::size_t found_qj = tx.find(":");
            if (found_qj != std::string::npos) {
                tr.qj = tx.substr(0, found_qj);
            } else {
                std::cerr << "ERROR #2" << std::endl;
                exit(-1);
            }
            // get old_gid
            std::size_t found_old_gid = tx.find("|");
            if (found_old_gid != std::string::npos) {
                tr.old_loc = tx.substr(found_qj + 1, found_old_gid - found_qj - 1);
            } else {
                std::cerr << "ERROR #3" << std::endl;
                exit(-1);
            }
            // get old_tid
            std::size_t found_old_tid = tx.find("|", found_old_gid + 1);
            if (found_old_tid != std::string::npos) {
                tr.old_tid = tx.substr(found_old_gid + 1, found_old_tid - found_old_gid - 1);
            } else {
                std::cerr << "ERROR #4" << std::endl;
                exit(-1);
            }
            // skip num_exons
            std::size_t found_ne = tx.find("|", found_old_tid + 1);
            // get fpkm
            std::size_t found_fpkm = tx.find("|", found_ne + 1);
            if (found_fpkm != std::string::npos) {
                if((found_fpkm - found_ne - 1)==0){
                    tr.fpkm = "0";
                }
                else{
                    tr.fpkm = tx.substr(found_ne + 1, found_fpkm - found_ne - 1);
                }
            } else {
                std::cerr << "ERROR #5" << std::endl;
                exit(-1);
            }
            // get tpm
            std::size_t found_tpm = tx.find("|", found_fpkm + 1);
            if (found_tpm != std::string::npos) {
                if((found_fpkm - found_ne - 1)==0){
                    tr.tpm = "0";
                }
                else {
                    tr.tpm = tx.substr(found_fpkm + 1, found_tpm - found_fpkm - 1);
                }
            } else {
                std::cerr << "ERROR #6" << std::endl;
                exit(-1);
            }
            // get cov
            std::size_t found_cov = tx.find("|", found_tpm + 1);
            if (found_cov != std::string::npos) {
                if((found_fpkm - found_ne - 1)==0){
                    tr.cov = "0";
                }
                else{
                    tr.cov = tx.substr(found_tpm + 1, found_cov - found_tpm - 1);
                }
            } else {
                std::cerr << "ERROR #7" << std::endl;
                exit(-1);
            }
            trs.push_back(tr);
        }
        delete tx_stream;
    }

    void add_tracking(std::string fname){
        std::string tissue = get_tissue_name(fname);

        // parse lines
        std::ifstream tissue_track_stream;
        tissue_track_stream.open(fname.c_str(),std::ios::in);
        if (!tissue_track_stream.good()){
            std::cerr<<"@ERROR::Couldn't open tracking file: "<<fname<<std::endl;
            exit(1);
        }
        std::ios::sync_with_stdio(false);

        std::string tline,uid,ulid,sid;
        std::vector<TR> trs;
        int line_n = 0;
        while (std::getline(tissue_track_stream,tline)) {
            break_tracking(trs,tline);
            line_n+=1;
//            std::cout<<line_n<<std::endl;
            for(auto tr : trs){
                uid = tissue+"_"+tr.qj+"_"+tr.old_tid;
                ulid = tissue+"_"+tr.qj+"_"+tr.old_loc;
                sid = tissue+"_"+tr.qj; // sample id
                // now that we have this information - we need to store it in a realtionship
//                std::cout<<uid<<"\t"<<tr.fpkm<<"\t"<<tr.tpm<<"\t"<<tr.cov<<std::endl;
                this->mtm_it = this->mtm.insert(std::make_pair(uid,MT(uid,tr.old_loc,std::stof(tr.fpkm),std::stof(tr.tpm),std::stof(tr.cov),tissue,sid)));
                // create an entry about the locus information here
                this->loci_sample_it = loci_sample.insert(std::make_pair(ulid,std::vector<MTM_IT>{}));
                this->loci_sample_it.first->second.push_back(this->mtm_it);
                // now to create/update entry in the tissue to sample map
                this->mttm_it = this->mttm.insert(std::make_pair(tr.new_tid,MTT(tr.new_tid,tr.new_loc,tissue)));
                this->mttm_it.first->second.add_tx(this->mtm_it);
                // add to tissue locus map for downstream stats computation
                this->loci_tissue_it = loci_tissue.insert(std::make_pair(tissue+"_"+tr.new_loc,std::vector<MTTM_IT>{}));
                this->loci_tissue_it.first->second.push_back(this->mttm_it);
                this->tissues.insert(tissue);
            }
            trs.clear();
        }
        tissue_track_stream.close();
    }

    void load_tt(){
        std::ifstream tissue_track_stream;
        tissue_track_stream.open(tissue_tracking_fname.c_str(),std::ios::in);
        if (!tissue_track_stream.good()){
            std::cerr<<"@ERROR::Couldn't open file with the list of sample paths: "<<tissue_tracking_fname<<std::endl;
            exit(1);
        }
        std::ios::sync_with_stdio(false);

        std::string aline,tissue;
        while (std::getline(tissue_track_stream,aline)) {
            std::cout<<aline<<std::endl;
            add_tracking(aline);
        }
        tissue_track_stream.close();
    }

    void load_at(){
        std::ifstream all_track_stream;
        all_track_stream.open(all_tracking_fname.c_str(),std::ios::in);
        if (!all_track_stream.good()){
            std::cerr<<"@ERROR::Couldn't open all tracking file: "<<all_tracking_fname<<std::endl;
            exit(1);
        }
        std::ios::sync_with_stdio(false);

        std::string aline;
        std::vector<TR> trs;
        while (std::getline(all_track_stream,aline)) {
            break_tracking(trs,aline);
            for(auto tr : trs){
                // find transcript
                this->mttm_it.first = this->mttm.find(tr.old_tid);
                if(this->mttm_it.first == this->mttm.end()){
                    std::cerr<<"transcript not found: "<<aline<<std::endl;
                    exit(-1);
                }
                else{
                    // add entry to the map linking all transcripts and
                    this->mattm_it = mattm.insert(std::make_pair(tr.new_tid,MATT(tr.new_tid,tr.new_loc)));
                    this->mattm_it.first->second.add_tx(this->mttm_it);

                    // also add to locus map linking all loci to all transcripts
                    this->loci_it = loci.insert(std::make_pair(tr.new_loc,std::make_pair(std::vector<MATTM_IT>{},-1)));
//                    std::cout<<tr.new_loc<<std::endl;
                    this->loci_it.first->second.first.push_back(this->mattm_it);
                }
            }
            trs.clear();
//            return;
        }
        all_track_stream.close();
    }

    void load_type(std::string gff_fname,int type){
        // these are gtf files, wo can simply use the gtf/gff mosules
        FILE* gff_file = fopen(gff_fname.c_str(), "r");
        if (gff_file == nullptr)
        {
            std::cerr << "@ERROR::Couldn't open annotation: " << gff_fname<< std::endl;
            exit(1);
        }
        GffReader gffReader;
        gffReader.init(gff_file,true);
        gffReader.readAll();

        GffObj *p_gffObj;

        for (int i = 0; i < gffReader.gflst.Count(); ++i){
            p_gffObj = gffReader.gflst.Get(i);
            if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0){
                continue;
            }

            // find the transcript in the index
            this->mattm_it.first = this->mattm.find(p_gffObj->getID());
            if(this->mattm_it.first==this->mattm.end()){
                std::cerr<<"transcript not found: "<<p_gffObj->getID()<<" in: "<<gff_fname<<std::endl;
                exit(-1);
            }
            this->mattm_it.first->second.set_type(type);

            // now also find the corresponding locus in the index and report cheange its type according to what's been precomputed
            this->loci_it.first = loci.find(p_gffObj->getGeneID());
            if(this->loci_it.first != loci.end()){
                if(type == TYPE::REAL_TX){
                    this->loci_it.first->second.second = 1;
                }
                else if(type == TYPE::SPLICING_TX){
                    this->loci_it.first->second.second = 1;
                }
                else if(type == TYPE::INTRONIC_TX){
                    this->loci_it.first->second.second = 1;
                }
                else if(type == TYPE::POLYMERASE_TX){
                    this->loci_it.first->second.second = 0;
                }
                else{
                    std::cerr<<"unknown type passed"<<std::endl;
                    exit(-1);
                }
            }
            else{ // danger zone
                std::cerr<<"locus: "<<p_gffObj->getGeneID()<<" from type is not found in index"<<std::endl;
                exit(-1);
            }
        }
    }

    std::string tissue_tracking_fname,all_tracking_fname,all_gtf_fname,true_gff_fname,polymerase_gff_fname,intronic_gff_fname,splicing_gff_fname;
};

enum Opt {TISSUES   = 't',
          GTF       = 'g',
          ALL       = 'a',
          OUTPUT    = 'o',
          REAL      = 'r',
          SPLICING  = 's',
          INTRONIC  = 'i',
          POLYMERASE= 'p'};

int main(int argc, char** argv) {

    ArgParse args("gtex_cor");
    args.add_string(Opt::TISSUES,"tissues","","File containing a list of paths to tracking files per tissue",true);
    args.add_string(Opt::GTF,"gtf","","File containing the merged gtf of all tissues",true);
    args.add_string(Opt::ALL,"all","","Path to the tracking for all samples across all tissues",true);
    args.add_string(Opt::OUTPUT,"output","","Basename for the output files",true);
    args.add_string(Opt::REAL,"real","","Real annotation for reporting stats on real transcripts",true);
    args.add_string(Opt::SPLICING,"splice","","Splicing noise gtf file",true);
    args.add_string(Opt::INTRONIC,"intron","","Intronic noise gtf file",true);
    args.add_string(Opt::POLYMERASE,"poly","","Polymerase noise gtf file",true);

    // we probably only the GTFs to get effective lengths of the transcripts, since we can get coverage and TPM information from the tissue tracking files

    if(argc <= 1 || strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    // first create the execution string
    std::string cl="gtex_stats ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    // now can load the tissue to sample transcript relationship table
    TrackingTree tt(args.get_string(Opt::TISSUES),args.get_string(Opt::ALL),args.get_string(Opt::GTF),args.get_string(Opt::REAL),args.get_string(Opt::SPLICING),args.get_string(Opt::INTRONIC),args.get_string(Opt::POLYMERASE));
    tt.load();

    Stats stats;
    tt.get_stats(stats);
    stats.save_stats(args.get_string(Opt::OUTPUT));

    return 0;
}

// another strategy - do not read individual gtfs

// start with tissue trackings from which you can get tpms, coverarages, tids, gids and lengths

// then proceed to load and link ALL.tracking

// eventually, read in the ALL.combined.gtf which can be used to get total intron length per transcript
// and the computed intron length can be subtracted from sample-specific transcripts to compute the effective length

// What questions are we trying to answer
// 1. distribution of the number of transcripts per locus
// 2. distribution of the number of noise transcripts per locus
// 3. distribution of the number of true transcripts per locus
// 4. distribution of coverage per true transcripts
// 5. distribution of coverage per noise transcripts
// 6. distribution of tpm per noise transcripts
// 7. distribution of tpm per true transcripts
// 8. total expression contribution for a given locus at each level
// 9. number of samples that contribute transcripts to a given locus (all and tissue levels)
// 10. number of tissues that contribute transcripts to a given locus (all level)


// Length is not always present in the .tracking files - seems like it is only there for class code "="
// instead let's forgo getting the number of reads (no need for length and effective length)
// and focus on coverage only. We can then target the same coverage distribution instead of reads distribution