#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "arg_parse.h"
#include "TrackingStats.h"
#include "TrackingTree.h"

enum Opt {TISSUES   = 't',
          GTF       = 'g',
          ALL       = 'a',
          OUTPUT    = 'o',
          REAL      = 'r',
          SPLICING  = 's',
          INTRONIC  = 'i',
          intergenic= 'p'};

int main(int argc, char** argv) {

    ArgParse args("gtex_cor");
    args.add_string(Opt::TISSUES,"tissues","","File containing a list of paths to tracking files per tissue",true);
    args.add_string(Opt::GTF,"gtf","","File containing the merged gtf of all tissues",true);
    args.add_string(Opt::ALL,"all","","Path to the tracking for all samples across all tissues",true);
    args.add_string(Opt::OUTPUT,"output","","Basename for the output files",true);
    args.add_string(Opt::REAL,"real","","Real annotation for reporting stats on real transcripts",true);
    args.add_string(Opt::SPLICING,"splice","","Splicing noise gtf file",true);
    args.add_string(Opt::INTRONIC,"intron","","Intronic noise gtf file",true);
    args.add_string(Opt::intergenic,"poly","","intergenic noise gtf file",true);

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
    TrackingTree tt(args.get_string(Opt::TISSUES),args.get_string(Opt::ALL),args.get_string(Opt::GTF),args.get_string(Opt::REAL),args.get_string(Opt::SPLICING),args.get_string(Opt::INTRONIC),args.get_string(Opt::intergenic));
    tt.load();

    TrackingStats stats;
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