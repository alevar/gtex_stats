//
// Created by sparrow on 3/4/20.
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

#include "TrackingStats.h"
#include "gff.h"
#include "GFaSeqGet.h"

#ifndef GTEX_STATS_TRACKINGTREE_H
#define GTEX_STATS_TRACKINGTREE_H

// Types of transcripts
enum TYPE {REAL_TX        = 1,
    SPLICING_TX = 2,
    INTRONIC_TX    = 3,
    intergenic_TX  = 4,
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

class TrackingTree {
public:
    TrackingTree(std::string ttf,std::string atf,std::string agf,std::string tr_tf,std::string sp_tf,std::string in_tf,std::string pl_tf):tissue_tracking_fname(ttf),all_tracking_fname(atf),all_gtf_fname(agf),true_gff_fname(tr_tf),
                                                                                                                                          splicing_gff_fname(sp_tf),intronic_gff_fname(in_tf),intergenic_gff_fname(pl_tf){}
    ~TrackingTree() = default;

    void load();

    void get_num_tx_tissue(TrackingStats& stats);

    void get_num_tx_per_tissue_locus5(TrackingStats& stats);

    void get_num_tx_per_sample_locus4(TrackingStats& stats);

    void get_num_tx_per_sample_locus5(TrackingStats& stats);

    void get_num_tx_sample(TrackingStats& stats);

    void get_num_tx_per_tissue_locus3(TrackingStats& stats);

    void get_real_noise_locs_sample(TrackingStats& stats);

    void get_real_noise_locs_tissue(TrackingStats& stats);

    void get_cov_sample(TrackingStats& stats);

    void get_tpm_fracs_joined(TrackingStats& stats);

    void get_cov_fracs_joined(TrackingStats& stats);

    void get_num_tx_per_sample_locus3(TrackingStats& stats);

    void get_num_tx_per_sample_locus2(TrackingStats& stats);

    void get_num_tx_per_tissue_locus2(TrackingStats& stats);

    void get_num_tx_per_sample_locus(TrackingStats& stats);

    void get_num_tx_per_tissue_locus(TrackingStats& stats);

    void get_loc_stats(TrackingStats& stats);

    void get_gauss_sample_tpm_per_tissue_loc(TrackingStats& stats);

    void get_gauss_sample_tx_per_tissue_loc(TrackingStats& stats);

    void get_num_tx_per_sample_locus6(TrackingStats& stats);

    void get_gauss_sample_per_tissue_loc(TrackingStats& stats);

    void get_gauss_sample_per_tissue_loc2(TrackingStats& stats);

    void get_stats(TrackingStats& stats){
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
//        get_num_tx_tissue(stats);
//
//        get_num_tx_per_tissue_locus5(stats);

//        get_loc_stats(stats);

//        get_gauss_sample_tx_per_tissue_loc(stats);
//        get_gauss_sample_tpm_per_tissue_loc(stats);
//        get_gauss_sample_per_tissue_loc(stats);

        get_num_tx_per_sample_locus5(stats);
        get_num_tx_per_sample_locus4(stats);
        get_real_noise_locs_tissue(stats);
        get_real_noise_locs_sample(stats);
        get_num_tx_sample(stats);
        get_num_tx_per_sample_locus6(stats);
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

    std::map<std::string,std::set<std::string>> tts; // tissue to sample map; stores the sample names for each tissue
    std::pair<std::map<std::string,std::set<std::string>>::iterator,bool> tts_it;

    int loc_is_real(std::vector<MATTM_IT>& txs);

    bool loc_tissue_is_real(std::vector<MTTM_IT>& txs);

    bool loc_sample_is_real(std::vector<MTM_IT>& txs);

    std::string get_tissue_name(std::string fname);

    void break_tracking(std::vector<TR>& trs,std::string& tline);

    void add_tracking(std::string fname);

    void load_tt();

    void load_at();

    void load_type(std::string gff_fname,int type);

    std::string tissue_tracking_fname,all_tracking_fname,all_gtf_fname,true_gff_fname,intergenic_gff_fname,intronic_gff_fname,splicing_gff_fname;
};


#endif //GTEX_STATS_TRACKINGTREE_H
