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
           NONINTRONIC_TX = 2,
           INTRONIC_TX    = 3,
           POLYMERASE_TX  = 4};

// Minimum Transcript
struct MT{
    float cov; // coverage
    float tpm; // tpm
    float fpkm; // fpkm
    std::string tid; // transcript id
    std::string locus; // locus id
    std::string tissue;
    std::string sample;
    int type = -1;
    MT(std::string tid,std::string locus,float fpkm,float tpm,float cov,std::string tissue,std::string sample):tid(tid), locus(locus), fpkm(fpkm), tpm(tpm), cov(cov), tissue(tissue),sample(sample) { }
    void set_cov(int cov){cov=cov;}
    void set_tpm(int tpm){tpm=tpm;}
    std::string get_strg(){
        return tissue+"\t"+sample+"\t"+tid+":"+locus+"|"+std::to_string(fpkm)+"|"+std::to_string(tpm)+"|"+std::to_string(cov);
    }
    void set_type(int type){
        if(this->type!=-1 && this->type!=type){
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
    std::vector<MTM_IT> txs;
    int type = -1;
    MTT(std::string tid,std::string locus):tid(tid),locus(locus){}
    void add_tx(MTM_IT mtm_it){
        this->txs.push_back(mtm_it);
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
        if(this->type!=-1 && this->type!=type){
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
    int type=-1;
    MATT(std::string tid,std::string loc):tid(tid),locus(loc){}
    void add_tx(MTTM_IT mttm_it){
        this->txs.push_back(mttm_it);
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
        if(this->type!=-1){
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
typedef std::map<std::string,std::vector<MATTM_IT>> Loci;
typedef std::pair<Loci::iterator,bool> Loci_IT;

struct TR{
    std::string new_tid,new_loc,old_tid,old_loc,qj,fpkm,tpm,cov;
    TR() = default;
    TR(std::string new_tid,std::string new_loc,std::string old_tid,std::string old_loc,std::string qj,std::string fpkm,std::string tpm,std::string cov):new_tid(new_tid),new_loc(new_loc),old_tid(old_tid),old_loc(old_loc),qj(qj),fpkm(fpkm),tpm(tpm),cov(cov){}
};

struct Stats{
    std::vector<int> num_tx_per_all_loc,num_tx_per_all_loc_real,num_tx_per_all_loc_intronic,num_tx_per_all_loc_nonintronic,num_tx_per_all_loc_polymerase,num_tx_per_all_loc_left;
    std::vector<float> sum_tx_per_all_loc,sum_tx_per_all_loc_real,sum_tx_per_all_loc_intronic,sum_tx_per_all_loc_nonintronic,sum_tx_per_all_loc_polymerase,sum_tx_per_all_loc_left;
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
        else if(num_tx_per_all_loc_intronic.size()!=num_tx_per_all_loc_nonintronic.size()) {
            std::cerr<<"all loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_all_loc_intronic: "<<num_tx_per_all_loc_intronic.size()<<std::endl;
            std::cerr<<"num_tx_per_all_loc_nonintronic: "<<num_tx_per_all_loc_nonintronic.size()<<std::endl;
            exit(-1);
        }
        else if(num_tx_per_all_loc_nonintronic.size()!=num_tx_per_all_loc_polymerase.size()) {
            std::cerr<<"all loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_all_loc_nonintronic: "<<num_tx_per_all_loc_nonintronic.size()<<std::endl;
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
            num_tx_per_all_loc_ss<<"all,real,intronic,nonintronic,polymerase,left"<<std::endl;

            for(int i=0;i<this->num_tx_per_all_loc.size();i++){
                num_tx_per_all_loc_ss<<num_tx_per_all_loc[i]<<","<<num_tx_per_all_loc_real[i]<<","<<num_tx_per_all_loc_intronic[i]<<","<<num_tx_per_all_loc_nonintronic[i]<<","<<num_tx_per_all_loc_polymerase[i]<<","<<num_tx_per_all_loc_left[i]<<std::endl;
            }

            num_tx_per_all_loc_ss.close();

            std::string sum_tx_per_all_loc_fname = base_out_fname+".sum_tx_per_all_loc";
            std::ofstream sum_tx_per_all_loc_ss(sum_tx_per_all_loc_fname.c_str());
            sum_tx_per_all_loc_ss<<"all,real,intronic,nonintronic,polymerase,left"<<std::endl;

            for(int i=0;i<this->sum_tx_per_all_loc.size();i++){
                sum_tx_per_all_loc_ss<<sum_tx_per_all_loc[i]<<","<<sum_tx_per_all_loc_real[i]<<","<<sum_tx_per_all_loc_intronic[i]<<","<<sum_tx_per_all_loc_nonintronic[i]<<","<<sum_tx_per_all_loc_polymerase[i]<<","<<sum_tx_per_all_loc_left[i]<<std::endl;
            }

            sum_tx_per_all_loc_ss.close();
        }
    }

    std::vector<int> num_tx_per_tissue_loc,num_tx_per_tissue_loc_real,num_tx_per_tissue_loc_intronic,num_tx_per_tissue_loc_nonintronic,num_tx_per_tissue_loc_polymerase,num_tx_per_tissue_loc_left;
    std::vector<float> sum_tx_per_tissue_loc,sum_tx_per_tissue_loc_real,sum_tx_per_tissue_loc_intronic,sum_tx_per_tissue_loc_nonintronic,sum_tx_per_tissue_loc_polymerase,sum_tx_per_tissue_loc_left;
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
        else if(num_tx_per_tissue_loc_intronic.size()!=num_tx_per_tissue_loc_nonintronic.size()) {
            std::cerr<<"tissue loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_tissue_loc_intronic: "<<num_tx_per_tissue_loc_intronic.size()<<std::endl;
            std::cerr<<"num_tx_per_tissue_loc_nonintronic: "<<num_tx_per_tissue_loc_nonintronic.size()<<std::endl;
            exit(-1);
        }
        else if(num_tx_per_tissue_loc_nonintronic.size()!=num_tx_per_tissue_loc_polymerase.size()) {
            std::cerr<<"tissue loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_tissue_loc_nonintronic: "<<num_tx_per_tissue_loc_nonintronic.size()<<std::endl;
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
            num_tx_per_tissue_loc_ss<<"all,real,intronic,nonintronic,polymerase,left"<<std::endl;

            for(int i=0;i<this->num_tx_per_tissue_loc.size();i++){
                num_tx_per_tissue_loc_ss<<num_tx_per_tissue_loc[i]<<","<<num_tx_per_tissue_loc_real[i]<<","<<num_tx_per_tissue_loc_intronic[i]<<","<<num_tx_per_tissue_loc_nonintronic[i]<<","<<num_tx_per_tissue_loc_polymerase[i]<<","<<num_tx_per_tissue_loc_left[i]<<std::endl;
            }

            num_tx_per_tissue_loc_ss.close();

            std::string sum_tx_per_tissue_loc_fname = base_out_fname+".sum_tx_per_tissue_loc";
            std::ofstream sum_tx_per_tissue_loc_ss(sum_tx_per_tissue_loc_fname.c_str());
            sum_tx_per_tissue_loc_ss<<"all,real,intronic,nonintronic,polymerase,left"<<std::endl;

            for(int i=0;i<this->sum_tx_per_tissue_loc.size();i++){
                sum_tx_per_tissue_loc_ss<<sum_tx_per_tissue_loc[i]<<","<<sum_tx_per_tissue_loc_real[i]<<","<<sum_tx_per_tissue_loc_intronic[i]<<","<<sum_tx_per_tissue_loc_nonintronic[i]<<","<<sum_tx_per_tissue_loc_polymerase[i]<<","<<sum_tx_per_tissue_loc_left[i]<<std::endl;
            }

            sum_tx_per_tissue_loc_ss.close();
        }
    }

    std::vector<int> num_tx_per_sample_loc,num_tx_per_sample_loc_real,num_tx_per_sample_loc_intronic,num_tx_per_sample_loc_nonintronic,num_tx_per_sample_loc_polymerase,num_tx_per_sample_loc_left;
    std::vector<float> sum_tx_per_sample_loc,sum_tx_per_sample_loc_real,sum_tx_per_sample_loc_intronic,sum_tx_per_sample_loc_nonintronic,sum_tx_per_sample_loc_polymerase,sum_tx_per_sample_loc_left;
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
        else if(num_tx_per_sample_loc_intronic.size()!=num_tx_per_sample_loc_nonintronic.size()) {
            std::cerr<<"sample loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_sample_loc_intronic: "<<num_tx_per_sample_loc_intronic.size()<<std::endl;
            std::cerr<<"num_tx_per_sample_loc_nonintronic: "<<num_tx_per_sample_loc_nonintronic.size()<<std::endl;
            exit(-1);
        }
        else if(num_tx_per_sample_loc_nonintronic.size()!=num_tx_per_sample_loc_polymerase.size()) {
            std::cerr<<"sample loc vectors not equal size"<<std::endl;
            std::cerr<<"num_tx_per_sample_loc_nonintronic: "<<num_tx_per_sample_loc_nonintronic.size()<<std::endl;
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
            num_tx_per_sample_loc_ss<<"all,real,intronic,nonintronic,polymerase,left"<<std::endl;

            for(int i=0;i<this->num_tx_per_sample_loc.size();i++){
                num_tx_per_sample_loc_ss<<num_tx_per_sample_loc[i]<<","<<num_tx_per_sample_loc_real[i]<<","<<num_tx_per_sample_loc_intronic[i]<<","<<num_tx_per_sample_loc_nonintronic[i]<<","<<num_tx_per_sample_loc_polymerase[i]<<","<<num_tx_per_sample_loc_left[i]<<std::endl;
            }

            num_tx_per_sample_loc_ss.close();

            std::string sum_tx_per_sample_loc_fname = base_out_fname+".sum_tx_per_sample_loc";
            std::ofstream sum_tx_per_sample_loc_ss(sum_tx_per_sample_loc_fname.c_str());
            sum_tx_per_sample_loc_ss<<"all,real,intronic,nonintronic,polymerase,left"<<std::endl;

            for(int i=0;i<this->sum_tx_per_sample_loc.size();i++){
                sum_tx_per_sample_loc_ss<<sum_tx_per_sample_loc[i]<<","<<sum_tx_per_sample_loc_real[i]<<","<<sum_tx_per_sample_loc_intronic[i]<<","<<sum_tx_per_sample_loc_nonintronic[i]<<","<<sum_tx_per_sample_loc_polymerase[i]<<","<<sum_tx_per_sample_loc_left[i]<<std::endl;
            }

            sum_tx_per_sample_loc_ss.close();
        }
    }

    std::vector<float> cov_sample_real,cov_sample_intronic,cov_sample_nonintronic,cov_sample_polymerase,cov_sample_left;
    std::vector<float> fpkm_sample_real,fpkm_sample_intronic,fpkm_sample_nonintronic,fpkm_sample_polymerase,fpkm_sample_left;
    std::vector<float> tpm_sample_real,tpm_sample_intronic,tpm_sample_nonintronic,tpm_sample_polymerase,tpm_sample_left;
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

        std::string cov_sample_nonintronic_fname = base_out_fname+".cov_sample_nonintronic";
        std::ofstream cov_sample_nonintronic_ss(cov_sample_nonintronic_fname.c_str());
        cov_sample_nonintronic_ss<<"cov,fpkm,tpm"<<std::endl;
        for(int i=0;i<this->cov_sample_nonintronic.size();i++){
            cov_sample_nonintronic_ss<<cov_sample_nonintronic[i]<<","<<fpkm_sample_nonintronic[i]<<","<<tpm_sample_nonintronic[i]<<std::endl;
        }
        cov_sample_nonintronic_ss.close();

        std::string cov_sample_polymerase_fname = base_out_fname+".cov_sample_polymerase";
        std::ofstream cov_sample_polymerase_ss(cov_sample_polymerase_fname.c_str());
        cov_sample_polymerase_ss<<"cov,fpkm,tpm"<<std::endl;
        for(int i=0;i<this->cov_sample_polymerase.size();i++){
            cov_sample_polymerase_ss<<cov_sample_polymerase[i]<<","<<fpkm_sample_polymerase[i]<<","<<tpm_sample_polymerase[i]<<std::endl;
        }
        cov_sample_polymerase_ss.close();
    }

    void save_stats(std::string base_out_fname){
        save_all_tx_stats(base_out_fname);
        save_tissue_tx_stats(base_out_fname);
        save_sample_tx_stats(base_out_fname);
        save_sample_cov(base_out_fname);
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

        // add all tracking
        std::cout<<"<<<loading all tracking info"<<std::endl;
        load_at();

        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading true transcript ids"<<std::endl;
        load_type(this->true_gff_fname,TYPE::REAL_TX);
        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading intronic noise transcript ids"<<std::endl;
        load_type(this->intronic_gff_fname,TYPE::INTRONIC_TX);
        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading non-intronic/splicing noise transcript ids"<<std::endl;
        load_type(this->splicing_gff_fname,TYPE::NONINTRONIC_TX);
        // load transcript IDs for ALL that are true
        std::cout<<"<<<loading poymerase noise transcript ids"<<std::endl;
        load_type(this->polymerase_gff_fname,TYPE::POLYMERASE_TX);
    }

    void get_stats(Stats& stats){
        std::cout<<"<<<computing stats"<<std::endl;
        // use TYPE tags to report stats per each type of transcript
        for(auto& v : this->loci){
            int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_nonintronic=0,num_txs_polymerase=0,num_txs_left=0;
            float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_nonintronic=0,sum_txs_polymerase=0,sum_txs_left=0;
            for(auto& tit : v.second){
                num_txs_all += tit.first->second.get_num_sample_txs();
                switch(tit.first->second.type){
                    case TYPE::REAL_TX:
                        num_txs_real += tit.first->second.get_num_sample_txs();
                        sum_txs_real += tit.first->second.get_sum_sample_tpms();
                        break;
                    case TYPE::INTRONIC_TX:
                        num_txs_intronic += tit.first->second.get_num_sample_txs();
                        sum_txs_intronic += tit.first->second.get_sum_sample_tpms();
                        break;
                    case TYPE::NONINTRONIC_TX:
                        num_txs_nonintronic += tit.first->second.get_num_sample_txs();
                        sum_txs_nonintronic += tit.first->second.get_sum_sample_tpms();
                        break;
                    case TYPE::POLYMERASE_TX:
                        num_txs_polymerase += tit.first->second.get_num_sample_txs();
                        sum_txs_polymerase += tit.first->second.get_sum_sample_tpms();
                        break;
                    default:
                        num_txs_left += tit.first->second.get_num_sample_txs();
                        sum_txs_left += tit.first->second.get_sum_sample_tpms();
                        break;
                }
            }
            stats.num_tx_per_all_loc_real.push_back(num_txs_real);
            stats.num_tx_per_all_loc_intronic.push_back(num_txs_intronic);
            stats.num_tx_per_all_loc_nonintronic.push_back(num_txs_nonintronic);
            stats.num_tx_per_all_loc_polymerase.push_back(num_txs_polymerase);
            stats.num_tx_per_all_loc_left.push_back(num_txs_left);
            stats.num_tx_per_all_loc.push_back(num_txs_all);

            stats.sum_tx_per_all_loc_real.push_back(sum_txs_real);
            stats.sum_tx_per_all_loc_intronic.push_back(sum_txs_intronic);
            stats.sum_tx_per_all_loc_nonintronic.push_back(sum_txs_nonintronic);
            stats.sum_tx_per_all_loc_polymerase.push_back(sum_txs_polymerase);
            stats.sum_tx_per_all_loc_left.push_back(sum_txs_left);
            stats.sum_tx_per_all_loc.push_back(sum_txs_all);
        }

        // Now do the same on a tissue level
        for(auto& v : this->loci_tissue){
            int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_nonintronic=0,num_txs_polymerase=0,num_txs_left=0;
            float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_nonintronic=0,sum_txs_polymerase=0,sum_txs_left=0;
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
                    case TYPE::NONINTRONIC_TX:
                        num_txs_nonintronic += sit.first->second.get_num_txs();
                        sum_txs_nonintronic += sit.first->second.get_sum_tpms();
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
            stats.num_tx_per_tissue_loc_nonintronic.push_back(num_txs_nonintronic);
            stats.num_tx_per_tissue_loc_polymerase.push_back(num_txs_polymerase);
            stats.num_tx_per_tissue_loc_left.push_back(num_txs_left);
            stats.num_tx_per_tissue_loc.push_back(num_txs_all);

            stats.sum_tx_per_tissue_loc_real.push_back(sum_txs_real);
            stats.sum_tx_per_tissue_loc_intronic.push_back(sum_txs_intronic);
            stats.sum_tx_per_tissue_loc_nonintronic.push_back(sum_txs_nonintronic);
            stats.sum_tx_per_tissue_loc_polymerase.push_back(sum_txs_polymerase);
            stats.sum_tx_per_tissue_loc_left.push_back(sum_txs_left);
            stats.sum_tx_per_tissue_loc.push_back(sum_txs_all);
        }

        // Now do the same on a sample level
        for(auto& v : this->loci_sample){
            int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_nonintronic=0,num_txs_polymerase=0,num_txs_left=0;
            float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_nonintronic=0,sum_txs_polymerase=0,sum_txs_left=0;
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
                    case TYPE::NONINTRONIC_TX:
                        num_txs_nonintronic += 1;
                        sum_txs_nonintronic += sit.first->second.tpm;
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
            stats.num_tx_per_sample_loc_nonintronic.push_back(num_txs_nonintronic);
            stats.num_tx_per_sample_loc_polymerase.push_back(num_txs_polymerase);
            stats.num_tx_per_sample_loc_left.push_back(num_txs_left);
            stats.num_tx_per_sample_loc.push_back(num_txs_all);

            stats.sum_tx_per_sample_loc_real.push_back(sum_txs_real);
            stats.sum_tx_per_sample_loc_intronic.push_back(sum_txs_intronic);
            stats.sum_tx_per_sample_loc_nonintronic.push_back(sum_txs_nonintronic);
            stats.sum_tx_per_sample_loc_polymerase.push_back(sum_txs_polymerase);
            stats.sum_tx_per_sample_loc_left.push_back(sum_txs_left);
            stats.sum_tx_per_sample_loc.push_back(sum_txs_all);
        }

        // Now to figure out how to report coverages and the rest
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
                case TYPE::NONINTRONIC_TX:
                    stats.cov_sample_nonintronic.push_back(sit.second.cov);
                    stats.fpkm_sample_nonintronic.push_back(sit.second.fpkm);
                    stats.tpm_sample_nonintronic.push_back(sit.second.tpm);
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
                tr.fpkm = tx.substr(found_ne + 1, found_fpkm - found_ne - 1);
            } else {
                std::cerr << "ERROR #5" << std::endl;
                exit(-1);
            }
            // get tpm
            std::size_t found_tpm = tx.find("|", found_fpkm + 1);
            if (found_tpm != std::string::npos) {
                tr.tpm = tx.substr(found_fpkm + 1, found_tpm - found_fpkm - 1);
            } else {
                std::cerr << "ERROR #6" << std::endl;
                exit(-1);
            }
            // get cov
            std::size_t found_cov = tx.find("|", found_tpm + 1);
            if (found_cov != std::string::npos) {
                tr.cov = tx.substr(found_tpm + 1, found_cov - found_tpm - 1);
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

        std::string tline,uid,ulid;
        std::vector<TR> trs;
        while (std::getline(tissue_track_stream,tline)) {
            break_tracking(trs,tline);
            for(auto tr : trs){
                uid = tissue+"_"+tr.qj+"_"+tr.old_tid;
                ulid = tissue+"_"+tr.qj+"_"+tr.old_loc;
                // now that we have this information - we need to store it in a realtionship
                this->mtm_it = this->mtm.insert(std::make_pair(uid,MT(uid,tr.old_loc,std::stof(tr.fpkm),std::stof(tr.tpm),std::stof(tr.cov),tissue,tr.qj)));
                // create an entry about the locus information here
                this->loci_sample_it = loci_sample.insert(std::make_pair(ulid,std::vector<MTM_IT>{}));
                this->loci_sample_it.first->second.push_back(this->mtm_it);
                // now to create/update entry in the tissue to sample map
                this->mttm_it = this->mttm.insert(std::make_pair(tr.new_tid,MTT(tr.new_tid,tr.new_loc)));
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
                    this->loci_it = loci.insert(std::make_pair(tr.new_loc,std::vector<MATTM_IT>{}));
                    this->loci_it.first->second.push_back(this->mattm_it);
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
        }
    }

    std::string tissue_tracking_fname,all_tracking_fname,all_gtf_fname,true_gff_fname,polymerase_gff_fname,intronic_gff_fname,splicing_gff_fname;

    // how to identify old transcripts (sample specific)???
    // create unique ID from tissue tracking as:
    // 1. tissue ID (requires a lookup index for tissue to id
    // 2. qJ
    // 3. transcript
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