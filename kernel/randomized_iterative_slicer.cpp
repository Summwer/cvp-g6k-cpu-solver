#include <iostream>
#include <stack>
#include "siever.h"
// #include "utils.h"
#include <time.h>
#include "/usr/include/libalglib/stdafx.h"
#include "/usr/include/libalglib/specialfunctions.h"

using namespace alglib;



// void Siever::update_vector(Vec<LFT> &vec){
//   vec.c = sim_hashes.compress(vec.v);
//   vec.norm = std::inner_product(vec.v.begin(), vec.v.begin()+full_n, vec.v.begin(),  static_cast<LFT>(0.0));
//   if(vec.norm < 0) {
//     std::cout<<vec.norm;
//     throw "The data type is too small to present the data!";
//   }
// }


// void Siever::initialize_pt(long* vec){
//   for(unsigned int i = 0; i < full_n; i++){
//     pt.x[i] = vec[i];
//   }
// }


void Siever::update_entry(Entry &e){
  e.c = sim_hashes.compress(e.yr);
  
  e.len = std::inner_product(e.yr.begin(), e.yr.begin()+n, e.yr.begin(),  static_cast<FT>(0.0)); //normalization

  // double rho;
  // int min_prec = gso_min_prec(rho, n, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec     = 53;
  int old_prec = FP_NR<mpfr_t>::set_prec(prec);
  e.len_prec = std::inner_product(e.yr.begin(), e.yr.begin()+n, e.yr.begin(),  static_cast<FP_NR<mpfr_t>>(0.0));
  
  if(e.len < 0) {
    std::cout<<e.len;
    throw "The data type is too small to present the data!";
  }
}


// Siever::Vec<LFT> Siever::recover_vector(Entry e, fplll::MatGSO<SZT, SFT> M, unsigned int l){ //,int q
//   Vec<LFT> vec;

//   //v = sum_i=l^d-1(x_i-l * pi_l(b_{i}))
//   vector<SFT> y = vector<SFT>(full_n,0.), z = vector<SFT>(full_n,0.);
//   SFT mu;

//   for(unsigned int i = l; i < full_n; i++){
//     y[i] = e.yr[i-l]/sqrt_rr[i-l];
//   }

//   M.to_canonical(z,y); 
  
//   for(unsigned int i = 0; i < full_n; i++){
//     vec.v[i] = z[i].get_d();
//   }
  
//   vec.close_x = e.x;
//   update_vector(vec);
//   return vec;
// }
 

// Siever::Vec<LFT> Siever::recover_vector(std::array<ZT,MAX_SIEVING_DIM> x, fplll::MatGSO<SZT, SFT> M, unsigned int l){ //,int q
//   Vec<LFT> vec;

//   //v = sum_i=l^d-1(x_i-l * pi_l(b_{i}))
//   vector<SFT> y = vector<SFT>(full_n,0.), z = vector<SFT>(full_n,0.);
//   SFT mu;
  
//   for(unsigned int i = l; i < full_n; i++){
//     // y[i] = yr[i]/sqrt_rr[i-l];
//     for(unsigned int j = l; j < i; j++){
//       M.get_mu(mu, i, j);
//       y[j].add(y[j], (double)x[i-l]*mu.get_d());
//     }
//     y[i] += (double) x[i-l];
//   }

 
//   M.to_canonical(z,y); 
//   // cout<<l<<endl;
//   // print_vector(z,0.,full_n);
//   for(unsigned int i = 0; i < full_n; i++){
//     vec.v[i] = z[i].get_d();
//   }
  
//   vec.close_x = x;
//   update_vector(vec);
//   return vec;
// }
 

//the db after sieve is the coordinate of basis, thus we should recover them.
// void Siever::recover_db(fplll::MatGSO<SZT, SFT> M, unsigned int l){ //int q,
//   vecs.resize(db.size());
//   for(unsigned int j = 0; j < db.size(); j++){
//     Entry e = db[cdb[j].i];
//     vecs[j] = recover_vector(e, M, l);// q,
//     // printf("yr=");
//     // vector<LFT> y = vector<LFT>(n,0);
//     // for(unsigned int i = 0; i < n; i++)
//     //   y[i] = e.yr[i]/sqrt_rr[i];
//     // print_vector(y,0,n);
//     // print_vector(sqrt_rr);
//   } 
// }


/*
//the db after sieve is the coordinate of basis, thus we should recover them, recover an input db
void Siever::recover_db(std::vector<std::vector<ZT>> input_db, fplll::MatGSO<SZT, fplll::FP_NR<FT>> M, int q,unsigned int full_n){
  initialize_local_params(M,0,0, full_n,full_n);
  //std::cout<<M.b[1][0]<<std::endl;
  unsigned long N = input_db.size();
  reserve(N);
  cdb.resize(N);
  db.resize(N);
  vecs.resize(N);
  //std::array<ZT,MAX_SIEVING_DIM> x {};
  Entry e;
  for(unsigned int j = 0; j < input_db.size(); j++){
    for(unsigned int i = 0; i < full_n; i++){
      e.x[i] = input_db[j][i];
    }  

    // if(j==db.size()-1){
    //   std::cout<<"x = ";
    //   for(unsigned int k = 0; k < full_n; k++){
    //     for(unsigned int t = 0; t < full_n; t++){
    //       //std::cout<<e.x[t]<<",";
    //       std::cout<<t<<","<<k<<",";
    //       std::cout<<M.b[t][k].get_si()<<" ";
    //     //std::cout<<vecs[db.size()-1].v[j]<<" ";
    //     }
    //     std::cout<<"--------------------"<<std::endl;
    //     break;
    //   }
      
    //   std::cout<<std::endl;
    // }
    
   
    //std::cout<<vecs[0].len<<std::endl;
    //return;
    int large = 0;
    size_t np = full_n / 2;
    np = (large > 3) ? 0 : np;
    //recompute_data_for_entry_babai<Recompute::recompute_all_and_consider_otf_lift>(e,np);
    recompute_data_for_entry<Recompute::recompute_all_and_consider_otf_lift>(e);

    //if (!uid_hash_table.insert_uid(e.uid)) continue;
    histo[histo_index(e.len)] ++;
    db[j] = e;
            
    CompressedEntry ce;
    ce.len = e.len;
    ce.c = e.c;
    ce.i = j;
    cdb[j] = ce;
    // std::cout<<ce.c<<std::endl;

    vecs[j] = recover_vector(e.x, M.b, q, full_n);

    // if(j==db.size()-1){
    //   std::cout<<"x = ";
    //   for(unsigned int k = 0; k < full_n; k++)
    //     std::cout<<e.x[k]<<" ";
    //   std::cout<<std::endl;
    // }
    

    
  } 
}
*/



//sample  a vector t' in L(A+t), just sample a vector in L(A), and plus it to t.
Entry Siever::sample_t_(Entry pt){ // int q,s
  Entry t_;
  t_ = sample(0); //large = 0 

  for(unsigned int i = 0; i < n; i++){
    t_.yr[i] = pt.yr[i] - t_.yr[i];
  }

  //recompute t_ in simhash and its norm. represent t as an Entry & CompressedEntry
  update_entry(t_);
 
  return t_;
}



// void Siever::initialize_local(fplll::MatGSO<SZT, SFT> M){
//     std::cout<<"--Initialize Local Parameters--"<<std::endl;
//     full_n = M.d;
//     FP_NR<mpfr_t> tmp;
//     FP_NR<mpfr_t> rri;
//     LFT rr;
//     FT* mu = new FT[full_n*full_n]; //initialize the 1-full_n array for mu
//     for(unsigned int i=0; i < full_n; i++){
//         for(unsigned int j=0; j < full_n; j++){
//             //std::cout<<M.get_mu(mu,i,j)<<" ";
//             if(i==j){
//               rr = M.get_r(rri,i,i).get_d();
//               mu[i * full_n + i] = rr;
//             }else{
//               M.get_mu(tmp,i,j);
//               mu[i * full_n + j] = tmp.get_d();
//             }
//         }
//     }
//     load_gso(full_n, mu);
//     initialize_local(0,full_n-30,full_n);

//     //3.initialize db
//     reset_stats(); //clear all previous statistics

// }


// void Siever::progressive_sieve(fplll::MatGSO<SZT, SFT> M, unsigned int l){
//     // initialize_local(M);
//     long max_db_size;
//     // std::cout<<"--Initialize Local Parameters--"<<std::endl;
//     // full_n = M.d;
//     // FP_NR<mpfr_t> tmp;
//     // FP_NR<mpfr_t> rri;
//     // 
//     // LFT rr;
//     // FT* mu = new FT[full_n*full_n]; //initialize the 1-full_n array for mu

//     // for(unsigned int i=0; i < full_n; i++){
//     //     for(unsigned int j=0; j < full_n; j++){
//     //         //std::cout<<M.get_mu(mu,i,j)<<" ";
//     //         if(i==j){
//     //           rr = M.get_r(rri,i,i).get_d();
//     //           mu[i * full_n + i] = rr;
//     //         }else{
//     //           M.get_mu(tmp,i,j);
//     //           mu[i * full_n + j] = tmp.get_d();
//     //         }
//     //     }
//     // }
//     // load_gso(full_n, mu);
//     // initialize_local(0,full_n-30,full_n);
    
//     // //3.initialize db
//     // reset_stats(); //clear all previous statistics
  
//     //4.sample vectors to grow db_size.
    
//     //grow db through dimension, set a progressive sieve.
//     cout<<"--Sieve Process--"<<endl;
//     while(n < full_n - l){
//       extend_left(1);
      
//       //5.Start Gauss Sieve.
//       if(n < params.gauss_crossover){
//         max_db_size = 500 + 10* n + 2 * params.db_size_factor * pow(params.db_size_base, n);
//         resize_db(max_db_size);
//         gauss_sieve(max_db_size);
//       }
//       else{
//         max_db_size = params.db_size_factor * pow(params.db_size_base, n);
//         LFT B = params.bgj1_bucket_size_factor * pow(max_db_size, params.bgj1_bucket_size_expo);
//         LFT y = B/(LFT) max_db_size;
//         LFT x = invincompletebeta((n+1.)/2., 0.5, y);
//         LFT alpha = pow(1. - x, 0.5);
//         resize_db(max_db_size);
//         hk3_sieve(alpha);
//       }
//     }
// }


void Siever::randomized_iterative_slicer(double* y, long* x, FT len_bound, int max_sample_times){ 
  update_entry(pt);
  
  if(max_sample_times == 0)
    max_sample_times = 1000;
  Entry t_new1, t_new2, t_ = pt, mint_ = pt;


  // //Compute  GH
  // vector<double> gs = vector<double>(n,0.);
  // FP_NR<mpfr_t> tmp;
  // for(unsigned int i = 0; i < n; i++){
  //   M.get_r(tmp,i,i);
  //   gs[i] = tmp.get_d();
  // }
  // gh = gaussian_heuristic(gs);
  // cout<<"gh = "<<gh<<endl;


  // cout<<"initial(mint_.yr) =";
  // for(int i = 0; i< n; i++)
  //   cout<<mint_.yr[i]<<" ";
  // cout<<endl;

  int sample_times = 1;
start_sample_t_:
  //3. sample  a vector t' in L(A+t), just sample a vector in L(A), and plus it to t. 
  if(max_sample_times > 0)
    t_ = sample_t_(pt);
  //4. sort db in order from nearby with pt to farway with pt. 
  //reduce the sample_vector by slicer algorithm, and get a short enough vector.
start_over:
  // std::cout<<t_.len<<std::endl;
  for(unsigned int i = 0; i < db.size(); i++){
    if( t_.len <= len_bound){
        break;
    }
    if(UNLIKELY(is_reducible_maybe<XPC_THRESHOLD>(t_.c,db[i].c))){
      // cout<<"reducible maybe"<<endl;
      for(unsigned int j = 0; j < n; j++){
        t_new1.yr[j] = t_.yr[j] + db[i].yr[j];
        t_new1.x[j] = t_.x[j] - db[i].x[j];
      }
      update_entry(t_new1);
      for(unsigned int  j = 0; j < n; j++){
        t_new2.yr[j] = t_.yr[j] - db[i].yr[j];
        t_new2.x[j] = t_.x[j] + db[i].x[j];
      }
      update_entry(t_new2);


      // cout<<"t_.x:";
      // for(unsigned int j = 0; j < n; j++)
      //   cout<<t_.x[j]<<" ";
      // cout<<endl;

      // cout<<"t_new1.len:"<< (t_new1.len < t_.len) << "t_new2.len_prec:"<<(t_new2.len < t_.len) <<", t_.len:" << t_.len<<endl; 
      // cout<<"t_new1.len_prec:"<< (t_new1.len_prec < t_.len_prec) << "t_new2.len_prec:"<<(t_new2.len_prec < t_.len_prec) <<", t_.len_prec:" << t_.len_prec<<endl; 
      // cout<<t_new1.len - t_.len<<endl;
      if(t_new1.len < t_.len){
        t_ = t_new1;
        goto start_over;
      }
      else if(t_new2.len < t_.len){
        t_ = t_new2;
        goto start_over;
      }    
    }
  }

  if ( t_.len < mint_.len)
    mint_ = t_;



  if( t_.len > len_bound and sample_times < max_sample_times){
    sample_times += 1;
    goto start_sample_t_;
  }



  
  if( mint_.len > pt.len and sample_times >= max_sample_times){
    mint_ = pt;
  }

  // if(sample_times > 1 and verbose and sample_times < max_sample_times)
  //   std::cout<<"\n Sample Times = "<<sample_times<<std::endl;
  // cout<<"mint_.yr =";
  // for(int i = 0; i< n; i++)
  //   cout<<mint_.yr[i]<<" ";
  // cout<<endl;
  recover_vector_from_yr(y, x, mint_);
  // cout<<"y = ";
  // for(int i = 0; i< full_n; i++)
  //   cout<<y[i]<<", ";
  // cout<<endl;
}





// //Call randomlized slicer algorithm to find the closest vector.
// void Siever::run_randslicer(long* pt, FT len_bound, int max_sample_times, long* &w, long* &ee){
//   vector<SFT> yl = vector<SFT>(M.d,0.); 
//   Entry t, e;
//   initialize_local(M);
//   progressive_sieve(M);  //extend to full_n, function progressive_sieve must stay before function construct_projected_entry, otherwise n<full_n
//   preprocess_vector(pt, M, yl);
//   construct_projected_entry(yl, t);
//   e = randomized_iterative_slicer(t, len_bound, max_sample_times,true);
//   recover_vector_from_yr(w, t, e, yl, M);

//   // Compute ee = t - closest_vector(t).
//   for(int i = 0; i < M.d; i++){
//     ee[i] = pt[i] - w[i];
//   }
// }