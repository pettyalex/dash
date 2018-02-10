#ifndef MAIN_H
#define MAIN_H

#include <fstream>
#include <iomanip>

using namespace std;

class Pedigree;

/// Parameter specifying the minimum density a cluster may be
extern float MIN_CLUSTER_DENSITY;
/// Parameter specifying the correlation at which clusters are printed
extern float MAX_DIFFERENCE_FOR_PRINT;
/// Parameter specifying minimum graph size to analyse
extern unsigned int MIN_GRAPH_SIZE;
/// Parameter specifying the window length
extern unsigned long WINDOW_SIZE;
/// Specifies current position in the window
extern unsigned long cur_pos;
/// Parameter specifying if haplotype mode (no splits) is on
extern bool HAP_MODE;
/// Parameter specifying if independent mode is on
extern bool INDEP_MODE;
/// Global access to Pedigree containing all individuals
extern Pedigree MAIN_PEDIGREE;
/// Global logging
extern ofstream of_log;
#endif
