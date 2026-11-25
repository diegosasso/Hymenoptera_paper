#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// The same f+ * exp(H+) calculation as the plain-R check in run-shuffles.R.
// Integer codes identify whole 15-unit patterns; code zero is unaccelerated.
// This helper only speeds up repeated calculations. It does not randomize data.
// [[Rcpp::export]]
NumericVector mean_neff_by_group(IntegerMatrix patterns, List active_edges) {
  int n_maps = patterns.ncol();
  int n_edges = patterns.nrow();
  int n_groups = active_edges.size();
  NumericVector result(n_groups);
  std::vector<int> counts(32768, 0);
  std::vector<int> used(n_edges);
  std::vector<double> count_log_count(n_edges + 1, 0.0);
  for (int count = 1; count <= n_edges; ++count)
    count_log_count[count] = count * std::log(double(count));

  for (int group = 0; group < n_groups; ++group) {
    IntegerVector edges = active_edges[group];
    for (int map = 0; map < n_maps; ++map) {
      int positive = 0;
      int n_used = 0;
      for (int i = 0; i < edges.size(); ++i) {
        int pattern = patterns(edges[i] - 1, map);
        if (pattern == 0) continue;
        ++positive;
        if (counts[pattern] == 0) used[n_used++] = pattern;
        ++counts[pattern];
      }
      double sum_count_log_count = 0.0;
      for (int i = 0; i < n_used; ++i) {
        int pattern = used[i];
        sum_count_log_count += count_log_count[counts[pattern]];
        counts[pattern] = 0;
      }
      if (positive > 0) {
        double entropy = std::log(double(positive)) -
                         sum_count_log_count / positive;
        result[group] += double(positive) / edges.size() * std::exp(entropy);
      }
    }
    result[group] /= n_maps;
  }
  return result;
}
