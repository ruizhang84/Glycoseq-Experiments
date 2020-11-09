#include <algorithm>
#include <fstream>
#include <iostream>
#include <thread>  
#include <mutex> 
#include <chrono> 
#include <map>
#include <unordered_set>

#include <argp.h>

#include "search_parameter.h"
#include "search_dispatcher.h"
#include "search_helper.h"

#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../../util/io/train_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/glycan/glycan_builder.h"
#include "../../engine/spectrum/normalize.h"
#include "../../engine/search/precursor_match.h"
#include "../../engine/search/spectrum_search.h"
#include "../../engine/search/search_result.h"
#include "../../engine/analysis/multi_comparison.h"
#include "../../engine/learn/neural_network.h"


const char *argp_program_version =
  "glycoseq v2.0";
const char *argp_program_bug_address =
  "<rz20@iu.edu>";

static char doc[] =
  "Glycoseq -- a program to search glycopeptide from high thoughput LS-MS/MS";

static struct argp_option options[] = {
    {"tpath", 't',    "sidentified_scans.csv",  0,  "trianing dataset, file name with identified scan numbers" },
    {"fpath", 'f',    "protein.fasta",  0,  "fasta, Protein Sequence Input Path" },
    {"output",    'o',    "weight.txt",   0,  "optimized score weight" },
    {"pthread",   'p',  "6",  0,  "Number of Searching Threads" },
    {"digestion",   'd',  "TG",  0,  "The Digestion, Trypsin (T), Pepsin (P), Chymotrypsin (C), GluC (G)" }, 
    {"miss_cleavage",   's',  "2",  0,  "The Missing Cleavage Upto" },    
    {"HexNAc",   'x',  "12",  0,  "Search Up to Number of HexNAc" },
    {"HexNA",   'y',  "12",  0,  "Search Up to Number of Hex" },
    {"Fuc",   'z',  "5",  0,  "Search Up to Number of Fuc" },
    {"NeuAc",   'u',  "4",  0,  "Search Up to Number of NeuAc" },
    {"NeuGc",   'w',  "0",  0,  "Search Up to Number of NeuGc" },
    {"ms1_tol",   'm',  "10",  0,  "MS Tolereance" },
    {"ms2_tol",   'n',  "0.01",  0,  "MS2 Tolereance" },
    {"ms1_by",   'k',  "0",  0, "MS Tolereance By Int: PPM (0) or Dalton (1)" },
    {"ms2_by",   'l',  "1",  0, "MS2 Tolereance By Int: PPM (0) or Dalton (1)" },
    { 0 }
};

static std::string default_train_path = 
        "/home/ruiz/Documents/Glycoseq-Experiments/data/identified_scan.csv";
static std::string default_fasta_path = 
        "/home/ruiz/Documents/Glycoseq-Experiments/data/haptoglobin.fasta";
static std::string default_out_path = "weight.txt";
static std::string default_digestion = "TG";

struct arguments
{
    char * train_path = const_cast<char*> (default_train_path.c_str());
    char * fasta_path = const_cast<char*> (default_fasta_path.c_str());
    char * out_path = const_cast<char*> (default_out_path.c_str());
    //digestion
    int miss_cleavage = 2;
    char * digestion = const_cast<char*> (default_digestion.c_str());
    // upper bound of glycan seaerch
    int n_thread = 6;
    int hexNAc_upper_bound = 12;
    int hex_upper_bound = 12;
    int fuc_upper_bound = 5;
    int neuAc_upper_bound = 4;
    int neuGc_upper_bound = 0;
    // searching precision
    double ms1_tol = 10;
    double ms2_tol = 0.01;
    int ms1_by = 0;
    int ms2_by = 1;
};


static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
    error_t err = 0;
    struct arguments *arguments =  static_cast<struct arguments*>(state->input);

    switch (key)
    {
    case 'd':
        arguments->digestion = arg;
        break;

    case 'f':
        arguments->fasta_path = arg;
        break;

    case 'i':
        arguments->train_path = arg;
        break;

    case 'k':
        arguments->ms1_by = atoi(arg);
        break;

    case 'l':
        arguments->ms2_by = atoi(arg);
        break;

    case 'm':
        arguments->ms1_tol = atof(arg);
        break;

    case 'n':
        arguments->ms2_tol = atof(arg);
        break;

    case 'o':
        arguments->out_path = arg;
        break;
    
    case 'p':
        arguments->n_thread = atoi(arg);
        break;
    

    case 's':
        arguments->miss_cleavage = atoi(arg);
        break;

    case 'u':
        arguments->neuAc_upper_bound = atoi(arg);
        break;

    case 'w':
        arguments->neuGc_upper_bound = atoi(arg);
        break;

    case 'x':
        arguments->hexNAc_upper_bound = atoi(arg);
        break;

    case 'y':
        arguments->hex_upper_bound = atoi(arg);
        break;

    case 'z':
        arguments->fuc_upper_bound = atoi(arg);
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return err;
}

static struct argp argp = { options, parse_opt, 0, doc };


SearchParameter GetParameter(const struct arguments& arguments)
{
    SearchParameter parameter;
    parameter.n_thread = arguments.n_thread;
    parameter.miss_cleavage = arguments.miss_cleavage;
    parameter.hexNAc_upper_bound = arguments.hexNAc_upper_bound;
    parameter.hex_upper_bound = arguments.hex_upper_bound;
    parameter.neuAc_upper_bound = arguments.neuAc_upper_bound;
    parameter.neuGc_upper_bound = arguments.neuGc_upper_bound;
    parameter.ms1_tol = arguments.ms1_tol;
    parameter.ms1_by = arguments.ms1_by == 0 ?
        algorithm::search::ToleranceBy::PPM :
        algorithm::search::ToleranceBy::Dalton;
    parameter.ms2_tol = arguments.ms2_tol;
    parameter.ms2_by = arguments.ms2_by == 0 ?
        algorithm::search::ToleranceBy::PPM :
        algorithm::search::ToleranceBy::Dalton;
    std::string protease(arguments.digestion);
    for(const char& c : protease)
    {
        switch (c)
        {
        case 'T': case 't':
            parameter.proteases.push_back(engine::protein::Proteases::Trypsin);
            break;

        case 'G': case 'g':
            parameter.proteases.push_back(engine::protein::Proteases::GluC);
            break;

        case 'P': case 'p':
            parameter.proteases.push_back(engine::protein::Proteases::Pepsin);
            break;
        case 'C': case 'c':
            parameter.proteases.push_back(engine::protein::Proteases::Chymotrypsin);
            break;

        default:
            break;
        }
    }
    parameter.weights.assign(5, 1.0);
    parameter.bias = 1.0;
    return parameter;
}

void collect_data(const std::vector<engine::search::SearchResult>& targets, 
                  const std::unordered_set<int>& scan_set,
                  std::vector<std::vector<double>>& X, std::vector<int>& y)
{
    for(const auto& it : targets)
    {
        std::vector<double> scores = it.Score();
        X.push_back(scores);
        if (scan_set.find(it.Scan()) != scan_set.end())
        {
            y.push_back(1);
        }
        else
        {
            y.push_back(0);
        }
    }
}

int main(int argc, char *argv[])
{
    // parse arguments
    struct arguments arguments;
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    std::string train_path(arguments.train_path);
    std::string fasta_path(arguments.fasta_path);
    std::string out_path(arguments.out_path);
    SearchParameter parameter = GetParameter(arguments);

    // Get train Dataset
    util::io::TrainReader train_reader(train_path);
    train_reader.Init();
    std::map<std::string, std::vector<int>> dataset = train_reader.Dataset();

    // neural network
    engine::learn::Classifier classifier;
    classifier.set_weight(parameter.weights);
    classifier.set_bias(parameter.bias);

    // read fasta and build peptides
    std::vector<std::string> peptides, decoy_peptides;
    std::unordered_set<std::string> seqs = PeptidesDigestion(fasta_path, parameter);
    peptides.insert(peptides.end(), seqs.begin(), seqs.end());

    // // build glycans
    std::unique_ptr<engine::glycan::NGlycanBuilder> builder =
        std::make_unique<engine::glycan::NGlycanBuilder>(parameter.hexNAc_upper_bound, 
            parameter.hex_upper_bound, parameter.fuc_upper_bound, 
                parameter.neuAc_upper_bound, parameter.neuGc_upper_bound);
    builder->Build();

    // search
    std::cout << "Start to train\n"; 
    auto start = std::chrono::high_resolution_clock::now();

    // training
    std::vector<std::vector<double>> X_train, X_test;
    std::vector<int> y_train, y_test;

    int count = 0;
    int size = dataset.size();
    for(const auto& data : dataset)
    {
        // read spectrum
        std::string spectra_path = data.first;
        std::vector<int> scans = data.second;
        std::unordered_set<int> scan_set(scans.begin(), scans.end());
        count++;

        std::unique_ptr<util::io::SpectrumParser> parser = 
            std::make_unique<util::io::MGFParser>(spectra_path, util::io::SpectrumType::EThcD);
        std::unique_ptr<util::io::SpectrumReader> spectrum_reader
            = std::make_unique<util::io::SpectrumReader>(spectra_path, std::move(parser));
        spectrum_reader->Init();

        // seraching targets 
        SearchDispatcher target_searcher(spectrum_reader->GetSpectrum(), builder.get(), peptides, parameter);
        std::vector<engine::search::SearchResult> targets = target_searcher.Dispatch();

        // for (auto& it : targets)
        // {
        //     it.set_simple(true);
        // }

        std::cout << "Total target:" << targets.size() << std::endl;

        // collect data;
        if (count < 0.6 * size)
        {
            collect_data(targets, scan_set, X_train, y_train);
        }
        else
        {
            collect_data(targets, scan_set, X_test, y_test);
        }
        
    }

    classifier.Train(X_train, y_train, 10, 0.01);

    // testing
    double accu = classifier.Test(X_test, y_test);
    std::cout << "Accuracy: " << accu << std::endl;

    // output
    std::vector<double> weight = classifier.Weight();
    double bias = classifier.Bias();

    std::ofstream outfile;
    outfile.open (out_path);
    outfile<< "weight: ";

    for(auto it : weight)
    {
        outfile<<it<<",";
    }
    outfile<<"bias: "<<bias<<"\n";
    outfile.close();

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << "Total Time: " << duration.count() << std::endl; 

}