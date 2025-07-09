#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>

// Check for OpenMP support
#ifdef _OPENMP
#include <omp.h>
#define PARALLEL_FOR _Pragma("omp parallel for")
#else
#define PARALLEL_FOR
#pragma message("OpenMP not supported. Compilation will proceed without parallel execution.")
#endif


// [[Rcpp::export]]
Rcpp::RObject allele_diff(std::vector<std::string> germs, 
                                           std::vector<std::string> inputs, 
                                           int X = 0, 
                                           bool parallel = false, 
                                           bool return_count = false) {
  // Set default non-mismatch characters
  std::unordered_set<char> non_mismatch_chars = {'N', '.', '-'};
  
 if (germs.size() != inputs.size()) {
   Rcpp::stop("The size of germs and inputs must be the same.");
 }
 
 size_t num_sequences = germs.size();
 auto pad_with_ns = [](std::string& seq, size_t target_length) {
   if (seq.size() < target_length) {
     seq.append(target_length - seq.size(), 'N');
   }
 };
 
 if (!parallel) {
   if (return_count) {
     std::vector<int> mutation_counts(num_sequences);
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       int count = 0;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           count++;
         }
       }
       mutation_counts[i] = count;
     }
     return Rcpp::wrap(mutation_counts);
   } else {
     Rcpp::List snp_list(num_sequences);
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       std::vector<int> snp_indices;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           snp_indices.push_back(j + 1);
         }
       }
       snp_list[i] = Rcpp::wrap(snp_indices);
     }
     return snp_list;
   }
 }
 
 if (parallel) {
   if (return_count) {
     std::vector<int> mutation_counts(num_sequences);
     PARALLEL_FOR
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       int count = 0;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           count++;
         }
       }
       mutation_counts[i] = count;
     }
     return Rcpp::wrap(mutation_counts);
   } else {
     Rcpp::List snp_list(num_sequences);
     PARALLEL_FOR
     for (size_t i = 0; i < num_sequences; ++i) {
       std::string germ = germs[i];
       std::string input = inputs[i];
       size_t max_length = std::max(germ.size(), input.size());
       pad_with_ns(germ, max_length);
       pad_with_ns(input, max_length);
       std::vector<int> snp_indices;
       for (size_t j = 0; j < max_length; ++j) {
         if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
             non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
             non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
           snp_indices.push_back(j + 1);
         }
       }
       snp_list[i] = Rcpp::wrap(snp_indices);
     }
     return snp_list;
   }
 }
 
 return Rcpp::wrap(R_NilValue);
}
