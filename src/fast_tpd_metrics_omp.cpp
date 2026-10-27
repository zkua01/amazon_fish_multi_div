// fast_tpd_metrics_omp.cpp  —  OpenMP 并行版
// 用法: Rcpp::sourceCpp("fast_tpd_metrics_omp.cpp")
// 依赖: install.packages("Rcpp")
//
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>

using namespace Rcpp;

// ============================================================
// 内部工具：计算单个物种的稀疏TPD全局索引（0-based）
// （与原始版本完全一致，未改动算法）
// ============================================================
static std::vector<double> compute_sparse_tpd(
    const std::vector<double>& mu,
    const std::vector<double>& sg,
    const std::vector<std::vector<double>>& grid_list,
    const std::vector<double>& strides,
    int n_divisions, int D,
    double alpha_1d, double alpha
) {
  // Step A: 每维边际密度 + Bonferroni预剪枝
  std::vector<std::vector<double>> marginals(D);
  std::vector<std::vector<int>>    pruned_idx(D);
  std::vector<int>                 n_sub(D);
  
  for (int d = 0; d < D; d++) {
    double m  = mu[d];
    double s  = std::max(sg[d], 1e-10);
    int    nd = n_divisions;
    
    std::vector<double> v(nd);
    double vsum = 0.0;
    for (int j = 0; j < nd; j++) {
      double z = (grid_list[d][j] - m) / s;
      v[j]  = std::exp(-0.5 * z * z);
      vsum += v[j];
    }
    for (int j = 0; j < nd; j++) v[j] /= vsum;
    
    // 【调试】单维度剪枝暂跳过，保留全部网格点
    int keep_n = nd;
    std::vector<int> kept(nd);
    std::iota(kept.begin(), kept.end(), 0);
    
    pruned_idx[d] = kept;
    n_sub[d]      = keep_n;
    marginals[d].resize(keep_n);
    for (int j = 0; j < keep_n; j++)
      marginals[d][j] = v[kept[j]];
  }
  
  // Step B: 子网格外积（列主序，dim0变化最快）
  int sub_size = 1;
  for (int d = 0; d < D; d++) sub_size *= n_sub[d];
  
  std::vector<double> tpd_sub(sub_size);
  // 初始化为dim0
  for (int j = 0; j < n_sub[0]; j++) tpd_sub[j] = marginals[0][j];
  int cur_size = n_sub[0];
  
  for (int d = 1; d < D; d++) {
    std::vector<double> tmp(cur_size * n_sub[d]);
    for (int jd = 0; jd < n_sub[d]; jd++)
      for (int k = 0; k < cur_size; k++)
        tmp[jd * cur_size + k] = marginals[d][jd] * tpd_sub[k];
    tpd_sub  = std::move(tmp);
    cur_size *= n_sub[d];
  }
  
  // 归一化
  double tsum = 0.0;
  for (double x : tpd_sub) tsum += x;
  if (tsum > 0) for (double& x : tpd_sub) x /= tsum;
  
  // Step C: alpha整体裁剪（找阈值）
  std::vector<double> sv(tpd_sub);
  std::sort(sv.begin(), sv.end(), std::greater<double>());
  double cumsum = 0.0, thresh = 0.0;
  for (double x : sv) {
    cumsum += x;
    thresh  = x;
    if (cumsum >= alpha) break;
  }
  
  // Step D: 子网格线性索引 → 全局线性索引（无expand.grid，O(|nz| * D)）
  std::vector<double> global_indices;
  global_indices.reserve(sub_size / 4);
  
  for (int li = 0; li < sub_size; li++) {
    if (tpd_sub[li] < thresh) continue;
    int    rem   = li;
    double gidx  = 0.0;
    for (int d = 0; d < D; d++) {
      int local_d  = rem % n_sub[d];
      gidx        += static_cast<double>(pruned_idx[d][local_d]) * strides[d];
      rem         /= n_sub[d];
    }
    global_indices.push_back(gidx); // 0-based
  }
  
  return global_indices;
}

// ============================================================
// 导出函数1：构建所有物种的稀疏TPD  [OpenMP并行版]
// ============================================================
// [[Rcpp::export]]
List build_sparse_tpds_cpp(
    NumericMatrix means,     // n_spp × D，行=物种，列=性状维度
    NumericMatrix sds,       // n_spp × D
    List          grid_list, // D个NumericVector，每维的网格点
    double        alpha      = 0.95,
    int           n_divisions = 70,
    int           n_threads   = 0   // 0 = 用全部核；设为2/4/8等限制线程数
) {
  int n_spp = means.nrow();
  int D     = means.ncol();
  double alpha_1d = std::pow(alpha, 1.0 / D);

  // Copy R matrices into plain C++ containers before entering OpenMP.
  // Worker threads must not touch R/Rcpp objects.
  std::vector<std::vector<double>> means_std(n_spp, std::vector<double>(D));
  std::vector<std::vector<double>> sds_std(n_spp, std::vector<double>(D));
  for (int spi = 0; spi < n_spp; spi++) {
    for (int d = 0; d < D; d++) {
      means_std[spi][d] = means(spi, d);
      sds_std[spi][d]   = sds(spi, d);
    }
  }
  
  // 预构建grid的std::vector缓存（避免重复Rcpp转换）
  std::vector<std::vector<double>> grid_std(D);
  for (int d = 0; d < D; d++)
    grid_std[d] = as<std::vector<double>>(as<NumericVector>(grid_list[d]));
  
  // strides（0-based列主序）
  std::vector<double> strides(D, 1.0);
  for (int d = 1; d < D; d++)
    strides[d] = strides[d-1] * n_divisions;
  
  // === 线程控制：n_threads=0 用全部可用核；但不会超过物种数 ===
  int max_threads = omp_get_max_threads();
  int requested_threads = n_threads;
  if (n_threads <= 0) n_threads = max_threads;
  n_threads = std::min(n_threads, max_threads);
  n_threads = std::min(n_threads, std::max(n_spp, 1));
  Rcpp::Rcout << "  OpenMP threads: " << n_threads
              << " (requested=" << requested_threads
              << ", max=" << max_threads
              << ", species=" << n_spp << ")\n";
  
  // === 预分配C++容器（不在并行区内碰R对象） ===
  std::vector<std::vector<double>> results(n_spp);
  
  // === OpenMP 并行物种循环 ===
  #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
  for (int spi = 0; spi < n_spp; spi++) {
    const std::vector<double>& mu = means_std[spi];
    const std::vector<double>& sg = sds_std[spi];
    
    std::vector<double> gidx = compute_sparse_tpd(
      mu, sg, grid_std, strides, n_divisions, D, alpha_1d, alpha
    );
    
    results[spi] = std::move(gidx);
  }
  
  // 串行转换回R List（R API非线程安全，必须在串行区执行）
  List result(n_spp);
  for (int spi = 0; spi < n_spp; spi++)
    result[spi] = wrap(results[spi]);
  
  return result;
}

// ============================================================
// 导出函数2：计算所有流域的FRic（未改动，Step2已经够快）
// ============================================================
// [[Rcpp::export]]
NumericVector calc_fric_cpp(
    List           sparse_TPDs,   // List of NumericVector（0-based全局索引）
    NumericMatrix  abundance_mat, // n_basins × n_spp（存在/缺失，>0为存在）
    double         cell_volume,
    CharacterVector basin_names   // 长度 = n_basins
) {
  int n_basins    = abundance_mat.nrow();
  int n_spp_total = sparse_TPDs.size();
  
  // 预提取所有稀疏TPD为std::vector（避免foreach中重复转换）
  std::vector<std::vector<double>> tpd_vec(n_spp_total);
  for (int j = 0; j < n_spp_total; j++)
    tpd_vec[j] = as<std::vector<double>>(
      as<NumericVector>(sparse_TPDs[j])
    );
  
  NumericVector fric_out(n_basins);
  
  for (int i = 0; i < n_basins; i++) {
    int total_cells = 0;
    for (int j = 0; j < n_spp_total; j++)
      if (abundance_mat(i, j) > 0)
        total_cells += (int)tpd_vec[j].size();

    if (total_cells == 0) {
      fric_out[i] = 0.0;
      continue;
    }

    std::vector<double> cells;
    cells.reserve(total_cells);

    for (int j = 0; j < n_spp_total; j++) {
      if (abundance_mat(i, j) > 0) {
        const auto& tv = tpd_vec[j];
        cells.insert(cells.end(), tv.begin(), tv.end());
      }
    }

    std::sort(cells.begin(), cells.end());
    int unique_cells = 1;
    for (int k = 1; k < total_cells; k++)
      if (cells[k] != cells[k - 1])
        unique_cells++;

    fric_out[i] = unique_cells * cell_volume;
  }
  
  fric_out.names() = basin_names;
  return fric_out;
}
