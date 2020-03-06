//
// Created by sparrow on 3/4/20.
//

#include <cmath>
#include "TrackingTree.h"

void TrackingTree::load(){
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
    std::cout<<"<<<loading non-intronic/splicing noise transcript ids"<<std::endl;
    load_type(this->splicing_gff_fname,TYPE::SPLICING_TX);
    // load transcript IDs for ALL that are true
    std::cout<<"<<<loading intronic noise transcript ids"<<std::endl;
    load_type(this->intronic_gff_fname,TYPE::INTRONIC_TX);
    // load transcript IDs for ALL that are true
    std::cout<<"<<<loading intergenic noise transcript ids"<<std::endl;
    load_type(this->intergenic_gff_fname,TYPE::intergenic_TX);
}

int TrackingTree::loc_is_real(std::vector<MATTM_IT>& txs){
    bool found_real = false;
    bool found_nonint = false;
    bool found_int = false;
    bool found_pol = false;
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
        if(tx.first->second.type == TYPE::intergenic_TX){
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

bool TrackingTree::loc_tissue_is_real(std::vector<MTTM_IT>& txs){
    for(auto& tx : txs){
        if(tx.first->second.type != TYPE::REAL_TX){
            return true;
        }
    }
    return false;
}

bool TrackingTree::loc_sample_is_real(std::vector<MTM_IT>& txs){
    for(auto& tx : txs){
        if(tx.first->second.type != TYPE::REAL_TX){
            return true;
        }
    }
    return false;
}

std::string TrackingTree::get_tissue_name(std::string fname){
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

void TrackingTree::break_tracking(std::vector<TR>& trs,std::string& tline){
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

void TrackingTree::add_tracking(std::string fname){
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
            this->tts_it = this->tts.insert(std::make_pair(tissue,std::set<std::string>{})); // add to tissue_sample map
            this->tts_it.first->second.insert(tr.qj);

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

void TrackingTree::load_tt(){
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

void TrackingTree::load_at(){
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
                this->loci_it.first->second.first.push_back(this->mattm_it);
            }
        }
        trs.clear();
    }
    all_track_stream.close();
}

void TrackingTree::load_type(std::string gff_fname,int type){
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
            else if(type == TYPE::intergenic_TX){
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

void TrackingTree::get_num_tx_tissue(TrackingStats& stats){
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
            else if(tx_type == TYPE::intergenic_TX){
                std::get<3>(tt_it.first->second)++;
            }
            else{
                continue;
            }
        }
    }
}

void TrackingTree::get_num_tx_per_tissue_locus5(TrackingStats& stats){
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
            else if(tx_type == TYPE::intergenic_TX){
                std::get<3>(ttloc5_it.first->second)++;
            }
            else{
                continue;
            }
        }
    }
}

void TrackingTree::get_num_tx_per_sample_locus4(TrackingStats& stats){
    std::cout<<"computing the number of transcripts per locus per sample4"<<std::endl;
    std::pair<std::map<std::string,std::tuple<int,int,int,int,std::string>>::iterator,bool> stloc_it;
    std::pair<std::map<std::string,std::array<int,4>>::iterator,bool> s2t_it;
    for(auto &atx : this->mattm){ // iterate over ALL loci
        int tx_type = atx.second.type;
        std::string atx_locus = atx.second.locus;
        for(auto& ttx : atx.second.txs){
            s2t_it = stats.sample2tissue_loc_txs.insert(std::make_pair(ttx.first->second.tissue+atx_locus,std::array<int,4>{0,0,0,0}));
            if(tx_type == TYPE::REAL_TX){
                s2t_it.first->second[0]++;
            }
            else if(tx_type == TYPE::SPLICING_TX){
                s2t_it.first->second[1]++;
            }
            else if(tx_type == TYPE::INTRONIC_TX){
                s2t_it.first->second[2]++;
            }
            else if(tx_type == TYPE::intergenic_TX){
                s2t_it.first->second[3]++;
            }
            else{
                continue;
            }
            for(auto& stx : ttx.first->second.txs){
                stloc_it = stats.num_tx_per_sample_loc4.insert(std::make_pair(stx.first->second.sample+atx_locus,std::make_tuple(0,0,0,0,ttx.first->second.tissue+atx_locus)));
                if(tx_type == TYPE::REAL_TX){
                    std::get<0>(stloc_it.first->second) += 1;
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    std::get<1>(stloc_it.first->second) += 1;
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    std::get<2>(stloc_it.first->second) += 1;
                }
                else if(tx_type == TYPE::intergenic_TX){
                    std::get<3>(stloc_it.first->second) += 1;
                }
                else{
                    continue;
                }
            }
        }
    }
}

void TrackingTree::get_num_tx_per_sample_locus5(TrackingStats& stats){
    std::cout<<"computing the number of transcripts per locus per sample5"<<std::endl;
    std::pair<std::map<std::pair<std::string,std::string>,std::tuple<std::vector<float>,std::vector<float>,std::vector<float>,std::vector<float>>>::iterator,bool> stloc5_it;
    for(auto &atx : this->mattm){ // iterate over ALL transcripts
        int tx_type = atx.second.type;
        std::string tx_locus = atx.second.locus;
        for(auto& ttx : atx.second.txs){
            for(auto& stx : ttx.first->second.txs){
                stloc5_it = stats.num_tx_per_sample_loc5.insert(std::make_pair(std::make_pair(stx.first->second.sample,tx_locus),std::make_tuple(std::vector<float>{},std::vector<float>{},std::vector<float>{},std::vector<float>{})));
                if(tx_type == TYPE::REAL_TX){
                    std::get<0>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    std::get<1>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    std::get<2>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::intergenic_TX){
                    std::get<3>(stloc5_it.first->second).push_back(stx.first->second.tpm);
                }
                else{
                    continue;
                }
            }
        }
    }
}

void TrackingTree::get_num_tx_sample(TrackingStats& stats){
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
                else if(tx_type == TYPE::intergenic_TX){
                    std::get<3>(st_it.first->second) += 1;
                }
                else{
                    continue;
                }
            }
        }
    }
}

void TrackingTree::get_num_tx_per_tissue_locus3(TrackingStats& stats){
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
                else if(tx.first->second.type == TYPE::intergenic_TX){
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

void TrackingTree::get_real_noise_locs_sample(TrackingStats& stats){
    std::cout<<"computing the number of real and noise ALL loci as identified in each samples"<<std::endl;
    std::unordered_map<std::string,std::array<int,9>> sample_locs; // real,splice,int,all_real,pol
    std::pair<std::unordered_map<std::string,std::array<int,9>>::iterator,bool> sit;
    for(auto& loc : this->loci){ // for now let's just output everything
        int loc_type = loc.second.second;
        std::map<std::string,std::tuple<int,int,int,int>> all_samples;
        std::pair<std::map<std::string,std::tuple<int,int,int,int>>::iterator,bool> as_it;
        for(auto& tx : loc.second.first){
            for(auto& sx : tx.first->second.txs){
                for(auto& sample : sx.first->second.samples){
                    as_it = all_samples.insert(std::make_pair(sample,std::make_tuple(0,0,0,0)));
                    if(sx.first->second.type==TYPE::REAL_TX){
                        std::get<0>(as_it.first->second)=1;
                    }
                    else if(sx.first->second.type == TYPE::SPLICING_TX){
                        std::get<1>(as_it.first->second)=1;
                    }
                    else if(sx.first->second.type == TYPE::INTRONIC_TX){
                        std::get<2>(as_it.first->second)=1;
                    }
                    else if(sx.first->second.type == TYPE::intergenic_TX){
                        std::get<3>(as_it.first->second)=1;
                    }
                }
            }
        }
        for(auto& sample : all_samples){
            // 0 real_only
            // 1 splicing_only
            // 2 int_only
            // 3 real_splice
            // 4 real_int
            // 5 splice_int
            // 6 real_splice_int
            // 7 all real
            // 8 pol
            if(loc_type==1){
                sit = sample_locs.insert(std::make_pair(sample.first,std::array<int,9>{0,0,0,0,0,0,0,0,0}));
                std::get<7>(sit.first->second)++; // all real
                if(std::get<0>(sample.second)==1 &&
                   std::get<1>(sample.second)==0 &&
                   std::get<2>(sample.second)==0){ // has real transcripts only
                    std::get<0>(sit.first->second)++;
                }
                else if(std::get<0>(sample.second)==0 &&
                        std::get<1>(sample.second)==1 &&
                        std::get<2>(sample.second)==0){ // has splicing transcripts only
                    std::get<1>(sit.first->second)++;
                }
                else if(std::get<0>(sample.second)==0 &&
                        std::get<1>(sample.second)==0 &&
                        std::get<2>(sample.second)==1){ // has intronic transcripts only
                    std::get<2>(sit.first->second)++;
                }
                else if(std::get<0>(sample.second)==1 &&
                        std::get<1>(sample.second)==1 &&
                        std::get<2>(sample.second)==0){ // has real and splicing transcripts
                    std::get<3>(sit.first->second)++;
                }
                else if(std::get<0>(sample.second)==1 &&
                        std::get<1>(sample.second)==0 &&
                        std::get<2>(sample.second)==1){ // has real and intronic transcripts
                    std::get<4>(sit.first->second)++;
                }
                else if(std::get<0>(sample.second)==0 &&
                        std::get<1>(sample.second)==1 &&
                        std::get<2>(sample.second)==1){ // has splicing and intronic transcripts
                    std::get<5>(sit.first->second)++;
                }
                else if(std::get<0>(sample.second)==1 &&
                        std::get<1>(sample.second)==1 &&
                        std::get<2>(sample.second)==1){ // has real splicing and intronic transcripts
                    std::get<6>(sit.first->second)++;
                }
                else{
                    std::cerr<<"unknown code"<<std::endl;
                    exit(-1);
                }

            }
            else if(loc_type==0){ // is 0
                sit = sample_locs.insert(std::make_pair(sample.first,std::array<int,9>{0,0,0,0,0,0,0,0,0}));
                std::get<8>(sit.first->second)++;
            }
            else{ // is undefined
                continue;
            }
        }
    }
    for(auto& tt : sample_locs){
        stats.num_locs_sample.push_back(tt.second);
    }
}

void TrackingTree::get_real_noise_locs_tissue(TrackingStats& stats){
    std::cout<<"computing the number of real and noise ALL loci per identified in each tissue"<<std::endl;
    std::unordered_map<std::string,std::pair<int,int>> tissue_locs;
    std::pair<std::unordered_map<std::string,std::pair<int,int>>::iterator,bool> tit;
    for(auto& loc : this->loci){ // for now let's just output everything
        int loc_type = loc.second.second;
        std::set<std::string> all_tissues;
        for(auto& tx : loc.second.first){
            all_tissues.insert(tx.first->second.tissues.begin(),tx.first->second.tissues.end());
        }
        for(auto& tissue : all_tissues){
            if(loc_type==1){
                tit = tissue_locs.insert(std::make_pair(tissue,std::make_pair(0,0)));
                tit.first->second.first++;
            }
            else if(loc_type==0){ // is 0
                tit = tissue_locs.insert(std::make_pair(tissue,std::make_pair(0,0)));
                tit.first->second.second++;
            }
            else{ // is undefined
                continue;
            }
        }
    }
    for(auto& tt : tissue_locs){
        stats.num_locs_tissue.push_back(std::make_tuple(tt.first,tt.second.first,tt.second.second));
    }
}

void TrackingTree::get_cov_sample(TrackingStats& stats){
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
            case TYPE::intergenic_TX:
                stats.cov_sample_intergenic.push_back(sit.second.cov);
                stats.fpkm_sample_intergenic.push_back(sit.second.fpkm);
                stats.tpm_sample_intergenic.push_back(sit.second.tpm);
                break;
            default:
                stats.cov_sample_left.push_back(sit.second.cov);
                stats.fpkm_sample_left.push_back(sit.second.fpkm);
                stats.tpm_sample_left.push_back(sit.second.tpm);
                break;
        }
    }
}

void TrackingTree::get_tpm_fracs_joined(TrackingStats& stats){
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

void TrackingTree::get_cov_fracs_joined(TrackingStats& stats){
    // what we can do is:
    // 1. for each locus save the mean real, mean intronic, mean splicing, mean intergenic coverages/tpms (sample level)
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

void TrackingTree::get_num_tx_per_sample_locus3(TrackingStats& stats){
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
                    else if(tx.first->second.type == TYPE::intergenic_TX){
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

void TrackingTree::get_num_tx_per_sample_locus2(TrackingStats& stats){
    // now to figure out how to report the average number of transcripts which a sample contributes to an tissue locus
    std::cout<<"computing the number of transcripts per sample loc2"<<std::endl;
    for(auto& lit : this->loci_tissue){
        int num_txs_total=lit.second.size();
        std::map<std::string,int> num_txs_real,num_txs_intronic,num_txs_splicing,num_txs_intergenic,num_txs_left;
        std::pair<std::map<std::string,int>::iterator,bool> nit;
        for(auto& tit : lit.second){
            // each tissue should only appear once in here - make sure that is so
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
                    case TYPE::intergenic_TX:
                        nit = num_txs_intergenic.insert(std::make_pair(ts,0));
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

        for(auto& nt : num_txs_intergenic){sum_pol+=nt.second;}
        int nt_pol = 0;
        if(sum_pol>0){
            nt_pol = (int)((float)sum_pol/(float)num_txs_intergenic.size());
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
        stats.num_tx_per_sample_loc_intergenic2.push_back(nt_pol);
        stats.num_tx_per_sample_loc_left2.push_back(nt_left);
        stats.num_tx_per_sample_loc2.push_back(num_txs_total);
    }
}

void TrackingTree::get_num_tx_per_tissue_locus2(TrackingStats& stats){
    // now to figure out how to report the average number of transcripts which a tissue contributes to an ALL locus
    std::cout<<"computing the number of transcripts per tissue loc2"<<std::endl;
    for(auto& lit : this->loci){
        int num_txs_total=lit.second.first.size();
        std::map<std::string,int> num_txs_real,num_txs_intronic,num_txs_splicing,num_txs_intergenic,num_txs_left;
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
                    case TYPE::intergenic_TX:
                        nit = num_txs_intergenic.insert(std::make_pair(ts,0));
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

        for(auto& nt : num_txs_intergenic){sum_pol+=nt.second;}
        int nt_pol = 0;
        if(sum_pol>0){
            nt_pol = (int)((float)sum_pol/(float)num_txs_intergenic.size());
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
        stats.num_tx_per_tissue_loc_intergenic2.push_back(nt_pol);
        stats.num_tx_per_tissue_loc_left2.push_back(nt_left);
        stats.num_tx_per_tissue_loc2.push_back(num_txs_total);
    }
}

void TrackingTree::get_num_tx_per_sample_locus(TrackingStats& stats){
    // Now do the same on a sample level
    std::cout<<"computing the number of transcripts per sample locus"<<std::endl;
    for(auto& v : this->loci_sample){
        int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_splicing=0,num_txs_intergenic=0,num_txs_left=0;
        float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_splicing=0,sum_txs_intergenic=0,sum_txs_left=0;
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
                case TYPE::intergenic_TX:
                    num_txs_intergenic += 1;
                    sum_txs_intergenic += sit.first->second.tpm;
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
        stats.num_tx_per_sample_loc_intergenic.push_back(num_txs_intergenic);
        stats.num_tx_per_sample_loc_left.push_back(num_txs_left);
        stats.num_tx_per_sample_loc.push_back(num_txs_all);

        stats.sum_tx_per_sample_loc_real.push_back(sum_txs_real);
        stats.sum_tx_per_sample_loc_intronic.push_back(sum_txs_intronic);
        stats.sum_tx_per_sample_loc_splicing.push_back(sum_txs_splicing);
        stats.sum_tx_per_sample_loc_intergenic.push_back(sum_txs_intergenic);
        stats.sum_tx_per_sample_loc_left.push_back(sum_txs_left);
        stats.sum_tx_per_sample_loc.push_back(sum_txs_all);
    }
}

void TrackingTree::get_num_tx_per_tissue_locus(TrackingStats& stats){
    // Now do the same on a tissue level
    std::cout<<"computing the number of transcripts per tissue locus"<<std::endl;
    for(auto& v : this->loci_tissue){
        int num_txs_all=0,num_txs_real=0,num_txs_intronic=0,num_txs_splicing=0,num_txs_intergenic=0,num_txs_left=0;
        float sum_txs_all=0,sum_txs_real=0,sum_txs_intronic=0,sum_txs_splicing=0,sum_txs_intergenic=0,sum_txs_left=0;
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
                case TYPE::intergenic_TX:
                    num_txs_intergenic += sit.first->second.get_num_txs();
                    sum_txs_intergenic += sit.first->second.get_sum_tpms();
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
        stats.num_tx_per_tissue_loc_intergenic.push_back(num_txs_intergenic);
        stats.num_tx_per_tissue_loc_left.push_back(num_txs_left);
        stats.num_tx_per_tissue_loc.push_back(num_txs_all);

        stats.sum_tx_per_tissue_loc_real.push_back(sum_txs_real);
        stats.sum_tx_per_tissue_loc_intronic.push_back(sum_txs_intronic);
        stats.sum_tx_per_tissue_loc_splicing.push_back(sum_txs_splicing);
        stats.sum_tx_per_tissue_loc_intergenic.push_back(sum_txs_intergenic);
        stats.sum_tx_per_tissue_loc_left.push_back(sum_txs_left);
        stats.sum_tx_per_tissue_loc.push_back(sum_txs_all);
    }
}

void TrackingTree::get_loc_stats(TrackingStats& stats){
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

void TrackingTree::get_gauss_sample_tpm_per_tissue_loc(TrackingStats& stats){
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

    // for intergenic we can simply compute the mean and std

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
                else if(tx_type == TYPE::intergenic_TX){
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

void TrackingTree::get_gauss_sample_tx_per_tissue_loc(TrackingStats& stats){
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
                else if(tx_type == TYPE::intergenic_TX){
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
                stats.gauss_sample_tx_per_tissue_loc.push_back(std::make_tuple(0,0,0));
                continue;
            }
            float mean = sum/total_num;
            float total = 0;
            for(auto& a : v.second[i]){
                total += ((float)a.second - mean)*((float)a.second - mean);
            }
            float sd = sqrt(total / total_num);
            stats.gauss_sample_tx_per_tissue_loc.push_back(std::make_tuple(mean,sd,total_num));
        }
    }
}

void TrackingTree::get_num_tx_per_sample_locus6(TrackingStats& stats){
    std::cout<<"computing the number of transcripts per locus per sample5"<<std::endl;

    // need to add 0 tpms for each transcript in tissue but not sample - this will help with the ordering among transcripts
    std::map<std::pair<std::string,std::string>,std::map<std::string,std::array<std::map<std::string,float>,4>>> tstt; // map of tissue loc to sample to map of transcript to float
    std::pair<std::map<std::pair<std::string,std::string>,std::map<std::string,std::array<std::map<std::string,float>,4>>>::iterator,bool> tstt_it;
    std::pair<std::map<std::string,std::array<std::map<std::string,float>,4>>::iterator,bool> stt_it;
    std::pair<std::map<std::string,float>::iterator,bool> tt_it;

    std::pair<std::map<std::pair<std::string,std::string>, // tissue and locus
            std::tuple<std::array<int,4>, // 4 integers contain total number of transcripts in a locus
                    std::array<std::set<std::string>,4>, // 4 integers contain number of transcripts in a locus of this tissue specifically
                    std::map<std::string, // sample name
                            std::array<std::vector<float>,4> // tpms of each type
                    > > >::iterator,bool> stloc6_it;

    std::pair<std::map<std::string,std::array<std::vector<float>,4>>::iterator,bool> sample_it;
    for(auto &atx : this->mattm){ // iterate over ALL transcripts
        // get total number of transcripts
        std::set<std::string> total_num_real,total_num_splice,total_num_int,total_num_pol;
        this->loci_it.first = this->loci.find(atx.second.locus);
        if(this->loci_it.first==this->loci.end()){
            std::cerr<<"locus not found"<<std::endl;
            exit(-1);
        }
        for(auto& sub_tx : this->loci_it.first->second.first){ // get total number of transcripts in the current locus
            if(sub_tx.first->second.type == TYPE::REAL_TX){
                total_num_real.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::SPLICING_TX){
                total_num_splice.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::INTRONIC_TX){
                total_num_int.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::intergenic_TX){
                total_num_pol.insert(sub_tx.first->first);
            }
            else{
                continue;
            }
        }

        int tx_type = atx.second.type;
        std::string tx_locus = atx.second.locus;
        for(auto& ttx : atx.second.txs){
            stloc6_it = stats.num_tx_per_sample_loc6.insert(std::make_pair(std::make_pair(ttx.first->second.tissue,tx_locus), // tissue and locus key
                                                                           std::make_tuple(std::array<int,4>{static_cast<int>(total_num_real.size()),
                                                                                                             static_cast<int>(total_num_splice.size()),
                                                                                                             static_cast<int>(total_num_int.size()),
                                                                                                             static_cast<int>(total_num_pol.size())}, // total numbers of transcripts
                                                                                           std::array<std::set<std::string>,4>{std::set<std::string>{},
                                                                                                                               std::set<std::string>{},
                                                                                                                               std::set<std::string>{},
                                                                                                                               std::set<std::string>{}}, // number of transcripts in a locus of this tissue specifically
                                                                                           std::map<std::string, // sample name
                                                                                                   std::array<std::vector<float>,4>>{})));
            tstt_it = tstt.insert(std::make_pair(std::make_pair(ttx.first->second.tissue,tx_locus), // tissue and locus key
                                                 std::map<std::string,std::array<std::map<std::string,float>,4>>{}));
            for(auto& stx : ttx.first->second.txs){
                stt_it = tstt_it.first->second.insert(std::make_pair(stx.first->second.sample,std::array<std::map<std::string,float>,4>{}));
                if(tx_type == TYPE::REAL_TX){
                    std::get<1>(stloc6_it.first->second)[0].insert(atx.second.tid);
                    stt_it.first->second[0].insert(std::make_pair(atx.second.tid,stx.first->second.tpm));
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    std::get<1>(stloc6_it.first->second)[1].insert(atx.second.tid);
                    stt_it.first->second[1].insert(std::make_pair(atx.second.tid,stx.first->second.tpm));
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    std::get<1>(stloc6_it.first->second)[2].insert(atx.second.tid);
                    stt_it.first->second[2].insert(std::make_pair(atx.second.tid,stx.first->second.tpm));
                }
                else if(tx_type == TYPE::intergenic_TX){
                    std::get<1>(stloc6_it.first->second)[3].insert(atx.second.tid);
                    stt_it.first->second[3].insert(std::make_pair(atx.second.tid,stx.first->second.tpm));
                }
                else{
                    continue;
                }
            }
        }
    }
    // now can add the tpms, and fill in the 0 for missing transcripts in samples
    for(auto& tl : stats.num_tx_per_sample_loc6){
        stloc6_it.first = stats.num_tx_per_sample_loc6.find(tl.first);

        tstt_it.first = tstt.find(tl.first);
        if(tstt_it.first == tstt.end()){
            std::cerr<<"something went wrong"<<std::endl;
        }
        for(auto& real_tid : std::get<1>(stloc6_it.first->second)[0]){
            for(auto& stt : tstt_it.first->second){ // iterate over samples for the current tissue
                sample_it = std::get<2>(stloc6_it.first->second).insert(std::make_pair(stt.first,std::array<std::vector<float>,4>{}));

                tt_it.first = stt.second[0].find(real_tid);
                if(tt_it.first == stt.second[0].end()){ // transcript not expressed in this sample
                    sample_it.first->second[0].push_back(0);
                }
                else{ // transcript is expressed in this sample
                    sample_it.first->second[0].push_back(tt_it.first->second);
                }
            }
        }
        for(auto& splicing_tid : std::get<1>(stloc6_it.first->second)[1]){
            for(auto& stt : tstt_it.first->second){ // iterate over samples for the current tissue
                sample_it = std::get<2>(stloc6_it.first->second).insert(std::make_pair(stt.first,std::array<std::vector<float>,4>{}));

                tt_it.first = stt.second[1].find(splicing_tid);
                if(tt_it.first == stt.second[1].end()){ // transcript not expressed in this sample
                    sample_it.first->second[1].push_back(0);
                }
                else{ // transcript is expressed in this sample
                    sample_it.first->second[1].push_back(tt_it.first->second);
                }
            }
        }
        for(auto& intronic_tid : std::get<1>(stloc6_it.first->second)[2]){
            for(auto& stt : tstt_it.first->second){ // iterate over samples for the current tissue
                sample_it = std::get<2>(stloc6_it.first->second).insert(std::make_pair(stt.first,std::array<std::vector<float>,4>{}));

                tt_it.first = stt.second[2].find(intronic_tid);
                if(tt_it.first == stt.second[2].end()){ // transcript not expressed in this sample
                    sample_it.first->second[2].push_back(0);
                }
                else{ // transcript is expressed in this sample
                    sample_it.first->second[2].push_back(tt_it.first->second);
                }
            }
        }
        for(auto& intergenic_tid : std::get<1>(stloc6_it.first->second)[3]){
            for(auto& stt : tstt_it.first->second){ // iterate over samples for the current tissue
                sample_it = std::get<2>(stloc6_it.first->second).insert(std::make_pair(stt.first,std::array<std::vector<float>,4>{}));

                tt_it.first = stt.second[3].find(intergenic_tid);
                if(tt_it.first == stt.second[3].end()){ // transcript not expressed in this sample
                    sample_it.first->second[3].push_back(0);
                }
                else{ // transcript is expressed in this sample
                    sample_it.first->second[3].push_back(tt_it.first->second);
                }
            }
        }
    }
}

void TrackingTree::get_gauss_sample_per_tissue_loc(TrackingStats& stats){
    typedef std::map<std::string,std::vector<float>> STPM; // holds a map of sample to all transcripts
    std::pair<STPM::iterator,bool> stpm_it;
    std::map<std::string,std::array<STPM,4>> tpm_txs;
    std::pair<std::map<std::string,std::array<STPM,4>>::iterator,bool> tp_it;

    typedef std::map<std::string,int> STX;
    std::pair<STX::iterator,bool> stx_it;
    std::map<std::string,std::tuple<std::array<STX,4>,std::array<STPM,4>,std::string,std::array<int,4>>> num_txs;
    std::pair<std::map< std::string,std::tuple< std::array<STX,4>,
            std::array<STPM,4>,
            std::string,
            std::array<int,4>
    >
    >::iterator,bool> nt_it;

    for(auto &atx : this->mattm){ // iterate over ALL transcripts
        int tx_type = atx.second.type;
        // we need to know how many transcripts of each type there were in total for this locus
        this->loci_it.first = this->loci.find(atx.second.locus);
        if(this->loci_it.first==this->loci.end()){
            std::cerr<<"locus not found"<<std::endl;
            exit(-1);
        }

        std::set<std::string> total_num_real,total_num_splice,total_num_int,total_num_pol;
        for(auto& sub_tx : this->loci_it.first->second.first){
            if(sub_tx.first->second.type == TYPE::REAL_TX){
                total_num_real.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::SPLICING_TX){
                total_num_splice.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::INTRONIC_TX){
                total_num_int.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::intergenic_TX){
                total_num_pol.insert(sub_tx.first->first);
            }
            else{
                continue;
            }
        }
        std::string tx_locus = atx.second.locus;
        for(auto& ttx : atx.second.txs){
            for(auto& stx : ttx.first->second.txs){
                nt_it = num_txs.insert(std::make_pair(stx.first->second.tissue+tx_locus,std::tuple<std::array<STX,4>,
                        std::array<STPM,4>,
                        std::string,std::array<int,4>>{std::array<STX,4>{STX{},STX{},STX{},STX{}},
                                                       std::array<STPM,4>{STPM{},STPM{},STPM{},STPM{}},
                                                       stx.first->second.tissue,
                                                       std::array<int,4>{static_cast<int>(total_num_real.size()),
                                                                         static_cast<int>(total_num_splice.size()),
                                                                         static_cast<int>(total_num_int.size()),
                                                                         static_cast<int>(total_num_pol.size())}}));
                if(tx_type == TYPE::REAL_TX){
                    stx_it = std::get<0>(nt_it.first->second)[0].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[0].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    stx_it = std::get<0>(nt_it.first->second)[1].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[1].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    stx_it = std::get<0>(nt_it.first->second)[2].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[2].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::intergenic_TX){
                    stx_it = std::get<0>(nt_it.first->second)[3].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[3].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else{
                    continue;
                }
            }
        }
    }

    // now we need to process these to get std and mean
    for(auto& v : num_txs){ // iterate over loci at the tissue level
        int num_samples = this->tts[std::get<2>(v.second)].size();

        for(int i=0;i<4;i++){ // iterate over the real,splice,intr,pol
            stats.gauss_sample_per_tissue_loc.push_back(std::pair<std::tuple<float,float,int,int>,std::vector<std::tuple<float,float,std::vector<float>>>>{});
            if(std::get<1>(v.second)[i].size()==0){ // no transcripts of the given type found for the tissue ALL locus
                stats.gauss_sample_per_tissue_loc.back().second.push_back(std::make_tuple(0,0,std::vector<float>{}));
                stats.gauss_sample_per_tissue_loc.back().first = std::make_tuple(0,0,0,std::get<3>(v.second)[i]);
                continue;
            }

            float sum = 0;
            for(auto& a : std::get<0>(v.second)[i]){ // iterate over the samples -> for each sample map inside STX
                sum+=a.second;
            }
            float mean = sum/num_samples;
            float total = 0;
            for(auto& a : std::get<0>(v.second)[i]){
                total += ((float)a.second - mean)*((float)a.second - mean);
            }
            for(int idx=0;idx<num_samples-std::get<0>(v.second)[i].size();idx++){ // this along with all_samples is here to correct gaussians for the absence of transcript from each group in other groups (such as real not being accounted for in splice, etc)
                total += (0.0 - mean)*(0.0 - mean);
            }
            float sd = sqrt(total/num_samples);
            stats.gauss_sample_per_tissue_loc.back().first = std::make_tuple(mean,sd,std::get<1>(v.second)[i].size(),std::get<3>(v.second)[i]);

            for(auto& a : std::get<1>(v.second)[i]){ // iterate over the transcripts for tpms
                float total_num_tpm = 0;
                float sum_tpm = 0;
                for(auto& t : a.second){
                    sum_tpm+=t;
                    total_num_tpm++;
                }
                float mean_tpm = 0;
                if(total_num_tpm>0){
                    mean_tpm = sum_tpm/total_num_tpm;
                }
                else{
                    std::cerr<<"TPMS exist but no tpm"<<std::endl;
                }
                // compute SD
                float total_tpm = 0;
                for(auto& t : a.second){
                    total_tpm += (t - mean_tpm)*(t - mean_tpm);
                }
                float sd_tpm = sqrt(total_tpm / total_num_tpm);
                stats.gauss_sample_per_tissue_loc.back().second.push_back(std::make_tuple(mean_tpm,sd_tpm,a.second));
            }
        }
    }
}

void TrackingTree::get_gauss_sample_per_tissue_loc2(TrackingStats& stats){
    typedef std::map<std::string,std::vector<float>> STPM; // holds a map of sample to all transcripts
    std::pair<STPM::iterator,bool> stpm_it;
    std::map<std::string,std::array<STPM,4>> tpm_txs;
    std::pair<std::map<std::string,std::array<STPM,4>>::iterator,bool> tp_it;

    typedef std::map<std::string,int> STX;
    std::pair<STX::iterator,bool> stx_it;
    std::map<std::string,std::tuple<std::array<STX,4>,std::array<STPM,4>,std::string,std::array<int,4>>> num_txs;
    std::pair<std::map< std::string,std::tuple< std::array<STX,4>,
            std::array<STPM,4>,
            std::string,
            std::array<int,4>
    >
    >::iterator,bool> nt_it;

    for(auto &atx : this->mattm){ // iterate over ALL transcripts
        int tx_type = atx.second.type;
        // we need to know how many transcripts of each type there were in total for this locus
        this->loci_it.first = this->loci.find(atx.second.locus);
        if(this->loci_it.first==this->loci.end()){
            std::cerr<<"locus not found"<<std::endl;
            exit(-1);
        }

        std::set<std::string> total_num_real,total_num_splice,total_num_int,total_num_pol;
        for(auto& sub_tx : this->loci_it.first->second.first){
            if(sub_tx.first->second.type == TYPE::REAL_TX){
                total_num_real.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::SPLICING_TX){
                total_num_splice.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::INTRONIC_TX){
                total_num_int.insert(sub_tx.first->first);
            }
            else if(sub_tx.first->second.type == TYPE::intergenic_TX){
                total_num_pol.insert(sub_tx.first->first);
            }
            else{
                continue;
            }
        }
        std::string tx_locus = atx.second.locus;
        for(auto& ttx : atx.second.txs){
            for(auto& stx : ttx.first->second.txs){
                nt_it = num_txs.insert(std::make_pair(stx.first->second.tissue+tx_locus,std::tuple<std::array<STX,4>,
                        std::array<STPM,4>,
                        std::string,std::array<int,4>>{std::array<STX,4>{STX{},STX{},STX{},STX{}},
                                                       std::array<STPM,4>{STPM{},STPM{},STPM{},STPM{}},
                                                       stx.first->second.tissue,
                                                       std::array<int,4>{static_cast<int>(total_num_real.size()),
                                                                         static_cast<int>(total_num_splice.size()),
                                                                         static_cast<int>(total_num_int.size()),
                                                                         static_cast<int>(total_num_pol.size())}}));
                if(tx_type == TYPE::REAL_TX){
                    stx_it = std::get<0>(nt_it.first->second)[0].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[0].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::SPLICING_TX){
                    stx_it = std::get<0>(nt_it.first->second)[1].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[1].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::INTRONIC_TX){
                    stx_it = std::get<0>(nt_it.first->second)[2].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[2].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else if(tx_type == TYPE::intergenic_TX){
                    stx_it = std::get<0>(nt_it.first->second)[3].insert(std::make_pair(stx.first->second.sample,0));
                    stx_it.first->second++;

                    stpm_it = std::get<1>(nt_it.first->second)[3].insert(std::make_pair(atx.second.tid,std::vector<float>{}));
                    stpm_it.first->second.push_back(stx.first->second.tpm);
                }
                else{
                    continue;
                }
            }
        }
    }

    // now we need to process these to get std and mean
    for(auto& v : num_txs){ // iterate over loci at the tissue level
        int num_samples = this->tts[std::get<2>(v.second)].size();

        for(int i=0;i<4;i++){ // iterate over the real,splice,intr,pol
            stats.gauss_sample_per_tissue_loc.push_back(std::pair<std::tuple<float,float,int,int>,std::vector<std::tuple<float,float,std::vector<float>>>>{});
            if(std::get<1>(v.second)[i].size()==0){ // no transcripts of the given type found for the tissue ALL locus
                stats.gauss_sample_per_tissue_loc.back().second.push_back(std::make_tuple(0,0,std::vector<float>{}));
                stats.gauss_sample_per_tissue_loc.back().first = std::make_tuple(0,0,0,std::get<3>(v.second)[i]);
                continue;
            }

            float sum = 0;
            for(auto& a : std::get<0>(v.second)[i]){ // iterate over the samples -> for each sample map inside STX
                sum+=a.second;
            }
            float mean = sum/num_samples;
            float total = 0;
            for(auto& a : std::get<0>(v.second)[i]){
                total += ((float)a.second - mean)*((float)a.second - mean);
            }
            for(int idx=0;idx<num_samples-std::get<0>(v.second)[i].size();idx++){ // this along with all_samples is here to correct gaussians for the absence of transcript from each group in other groups (such as real not being accounted for in splice, etc)
                total += (0.0 - mean)*(0.0 - mean);
            }
            float sd = sqrt(total/num_samples);
            stats.gauss_sample_per_tissue_loc.back().first = std::make_tuple(mean,sd,std::get<1>(v.second)[i].size(),std::get<3>(v.second)[i]);

            for(auto& a : std::get<1>(v.second)[i]){ // iterate over the transcripts for tpms
                float total_num_tpm = 0;
                float sum_tpm = 0;
                for(auto& t : a.second){
                    sum_tpm+=t;
                    total_num_tpm++;
                }
                float mean_tpm = 0;
                if(total_num_tpm>0){
                    mean_tpm = sum_tpm/total_num_tpm;
                }
                else{
                    std::cerr<<"TPMS exist but no tpm"<<std::endl;
                }
                // compute SD
                float total_tpm = 0;
                for(auto& t : a.second){
                    total_tpm += (t - mean_tpm)*(t - mean_tpm);
                }
                float sd_tpm = sqrt(total_tpm / total_num_tpm);
                stats.gauss_sample_per_tissue_loc.back().second.push_back(std::make_tuple(mean_tpm,sd_tpm,a.second));
            }
        }
    }
}