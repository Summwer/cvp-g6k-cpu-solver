// #include "randomized_iterative_slicer.h"
#include <iostream>
#include <stack>
#include "siever.h"
// #include "utils.h"
#include <time.h>
#include "/usr/include/libalglib/stdafx.h"
#include "/usr/include/libalglib/specialfunctions.h"

using namespace alglib;

//recover vector from projected vector to find a close vector to t
//cvp extend left
void Siever::cvp_extend_left(Entry &e, unsigned int lp){
  // SFT mu,tmp;
  // initialize_local(ll, l - lp, r);
  CPUCOUNT(202);

  assert(lp <= l);
  if(lp==0) return;
  initialize_local(ll, l - lp, r);
  
  std::copy_backward(e.x.begin(), e.x.begin()+n-lp, e.x.begin()+n);
  std::fill(e.x.begin(), e.x.begin()+lp, 0);
  std::fill(e.x.begin()+n,e.x.end(),0);
  for (int i = lp - 1 ; i >= 0; --i){
    e.yr[i] = std::inner_product(e.x.cbegin()+i+1, e.x.cbegin()+n, muT[i].cbegin()+i+1, static_cast<FT>(0.));
    e.x[i] = - (ZT)std::round(e.yr[i] - yl[i+ll]);
    e.yr[i] += e.x[i];
    // e.yr[i] *= sqrt_rr[i];
  }
  for (unsigned int i = 0 ; i < n; i++){
    e.yr[i] = std::inner_product(e.x.cbegin()+i, e.x.cbegin()+n, muT[i].cbegin()+i,  static_cast<FT>(0.)); //* sqrt_rr[i];
  }
}


// //extend yr from projected lattice [l:r] to lattice [0:r], the coeffiencts on full GS-basis
// void Siever::left_recompute_yr(Entry &e, unsigned int lp){
//   CPUCOUNT(202);

//   assert(lp == ll);
//   if(lp==0) return;
//   initialize_local(0, 0, r);
//   std::copy_backward(e.x.begin(), e.x.begin()+n-lp, e.x.begin()+n);
//   std::fill(e.x.begin(), e.x.begin()+lp, 0);
//   std::fill(e.x.begin()+n,e.x.end(),0);
//   cout<<"e.yr: ";
//   for (unsigned int i = 0; i < n; ++i){
//     e.yr[i] = std::inner_product(e.x.cbegin(), e.x.cbegin()+n, muT[i].cbegin(),  static_cast<FT>(0.)) ; //* sqrt_rr[i]
//     cout<<"("<<e.yr[i]<<", "<<e.x[i]<<", ";
//   }
//   cout<<endl;
// }

//y: coefficients of projected close vector on projected gs 
//x: coefficients of projected close vector on lattice basis
void Siever::recover_vector_from_yr(double* y, long* x, Entry pe){
  //Recover the closest vector to a full-dimensional form.
  // vector<SFT> y = vector<SFT>(full_n,0.), z = vector<SFT>(full_n,0.);
  std::fill(y, &y[full_n], 0.);
  std::fill(x, &x[full_n], 0); 
  //cv: the full-dimensional closest vector in lattice to t.
  Entry cv;
  // cout<<"cv.yr:";
  // cout<<"pcv.x: [";
  for(unsigned int i = 0; i < n; i++){
    cv.yr[i] = pt.yr[i] - pe.yr[i];
    cv.x[i] = pe.x[i];
    // cout<<cv.x[i]<<", ";
    // cout<<cv.yr[i]<<" ";
    // cout<<pe.yr[i]<<" ";
  }
  // cout<<"]"<<endl;

  // cout<<endl;

  // cout<<"l:"<<l<<", ll: "<<ll<<endl;
  cvp_extend_left(cv, l-ll); //actual dim - sieve dim

  
  // left_recompute_yr(cv, ll); //extend left to ll.
  
  //Recover cv from coordinates on gso/gh.
  for(unsigned int i = 0; i < n; i++){
    y[i+ll] = cv.yr[i] / sqrt_rr[i]; 
  }

  // cout<<"y:";
  // for(unsigned int i = 0; i < full_n; i++)
  //   cout<<y[i]<<" ";
  // cout<<endl;

  for(unsigned int i = 0; i < n; i++){
    x[i+ll] = (int) cv.x[i];//sqrt_rr[i]; 
  }

  // cout<<"x:";
  // for(unsigned int i = 0; i < full_n; i++)
  //   cout<<x[i]<<" ";
  // cout<<endl;
}


// void Siever::initialize_target_vector(long* target_vector, fplll::MatGSO<SZT, SFT> M, vector<SFT> &yl){
//     vector<SFT> z = vector<SFT>(full_n,0.); //use yl to store the <t,bi*>/<bi*,bi*>
//     //convert the type from float to SFT
//     for(unsigned int i = 0; i < full_n; i++){
//       z[i] = target_vector[i]; 
//     }
//     M.from_canonical(yl,z);  //y is the coeefficient of target_vector wrt B*.
// }


// void Siever::initialize_target_vector(long* target_vector,  vector<SFT> &yl){
//     vector<SFT> z = vector<SFT>(full_n,0.); //use yl to store the <t,bi*>/<bi*,bi*>
//     //convert the type from float to SFT
//     for(unsigned int i = 0; i < full_n; i++){
//       z[i] = target_vector[i]; 
//     }
//     M.from_canonical(yl,z);  //y is the coeefficient of target_vector wrt B*.
// }


//yl: coeeficients of target vector wrt MatGSO basis
//pt: projected target vector
void Siever::initialize_projected_target_vector(){ 
  // Entry pt; 
  // for(unsigned int i = 0; i < n; i++){
    // pt.x[i] = 0;
    // pt.yr[i] = yl[i+full_n-n] * sqrt_rr[i];
  // }
 
  for(unsigned int i = 0; i < n; i++){
    pt.x[i] = 0;
    pt.yr[i] = yl[i+l] * sqrt_rr[i];  
    //cout<<l<<","<<pt.yr[i] <<"," <<yl[i+l] <<","<< sqrt_rr[i]<<endl;
  }
  update_entry(pt);
}


// void Siever::progressive_slicer_with_d4f(vector<VEC_ZT> target_vector, FT len_bound, int max_sample_times, fplll::MatGSO<SZT, SFT> M, unsigned int f, vector<VEC_ZT> &w, vector<VEC_ZT> &ee){  
//   /*
//   input: target_vector: the target vector t that |t-w|<γ·λ1(L)
//          w: the closest vector to t, the output vector
//          len_bound(γ): |t-w|<γ·λ1(L)
//          max_sample_times: the maximum sample times for randomized slicer
//          M: MatGSO class for lattice basis B
//          f: dimension-for-free value, π_[f:d](L) 
//   */ 
//     long max_db_size;
//     full_n = M.d;
//     vector<SFT> yl = vector<SFT>(full_n,0.);

//     //Compute full-dimensional GH
//     vector<double> gs = vector<double>(full_n,0.);
//     FP_NR<mpfr_t> tmp;
//     for(unsigned int i = 0; i < full_n; i++){
//       M.get_r(tmp,i,i);
//       gs[i] = tmp.get_d();
//     }
//     double gh = gaussian_heuristic(gs);


//     initialize_local(M);
//     preprocess_vector(target_vector, M, yl);

//     cout<<"--Progressive Slicer Process--"<<endl;
//     Entry pt;
//     while(n <= full_n ){
      
      
//       //5.Start Sieve.
//       if(n < params.gauss_crossover){
//         //Gauss Sieve
//         max_db_size = 500 + 10 * n + 2 * params.db_size_factor * pow(params.db_size_base, n);
//         // max_db_size = params.db_size_factor * pow(params.db_size_base, n);
//         // resize_db(max_db_size);
//         gauss_sieve(max_db_size);
//       }
//       else{
//         // max_db_size = params.db_size_factor * pow(params.db_size_base, n);
//         max_db_size = params.db_size_factor * pow(params.db_size_base, n+4);
//         LFT B = params.bgj1_bucket_size_factor * pow(max_db_size, params.bgj1_bucket_size_expo);
//         LFT y = B/(LFT) max_db_size;
//         // LFT x = invincompletebeta((n+1.)/2., 0.5, y);
//         LFT x = invincompletebeta((n+5.)/2., 0.5, y);
//         LFT alpha = pow(1. - x, 0.5);
//         resize_db(max_db_size);
//         hk3_sieve(alpha);
//       }

//       if(true){
//         //4. Compute projected t, i.e. pt = pi_l(t)
//         construct_projected_entry(yl, pt);

//         Entry pe = randomized_iterative_slicer(pt, len_bound, max_sample_times, true); 
//         // cout<<"pe.len = "<<pe.len<<endl;
//         // cout<<"n/(full_n-f) = "<<n<<"/"<<full_n-f<<endl;
//         if(pe.len <= len_bound){
//           cout<<n<<"/"<<full_n<<": "<< db.size() <<endl;
//           // cout<<"n/(full_n-f) = "<<n<<"/"<<full_n-f<<endl;
//           //6. recover the projected closest vector to a full-dimensional state.
//           recover_vector_from_yr(w, pt, pe, yl, M);

//           // Compute ee = t - closest_vector(t).
//           for(unsigned int i = 0; i < full_n; i++){
//             ee[i] = target_vector[i] - w[i];
//           }
//           print_vector(ee,0,full_n);
//           if(compute_norm(ee) <= len_bound * gh){
//             cout<<"n/(full_n-f) = "<<n<<"/"<<full_n-f<<endl;
//             cout<<"norm(ee) = "<<compute_norm(ee)<< "<= γ^2 gh = "<< len_bound * gh <<endl;
//             return;
//           }
//         }
//       }

//       //7.sample vectors and sieve to grow db_size.
//       //grow db through dimension, set a progressive sieve.
//       // initialize_local(0,ff,full_n);
//       if(n==full_n-f)//full-dimension
//         break;
//       extend_left(1);
//     }
//     cout<<"n/(full_n-f) = "<<n<<"/"<<full_n-f<<endl;
//     cout<<"Fail to find the solution"<<endl;
// }
