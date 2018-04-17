// Copyright (C) 2018 by Fred Hutchinson Cancer Research Center
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <exception>
#include <algorithm>
#include <string>
#include <vector>
#include <functional>
#include <ctime>
#include <cstdio>
#include "./flowReMix.h"

#ifdef PROFILE
std::vector<double> timings(1);
#endif

/* some convenience functions for debugging */
void printDims(arma::vec a, std::string c) {
  Rcout << c << " " << a.n_rows << " x " << a.n_cols << "\n";
};
void printDims(arma::uvec a, std::string c) {
  Rcout << c <<" " << a.n_rows << " x " << a.n_cols << "\n";
};

void printDims(arma::mat a, std::string c) {
  Rcout << c << " " << a.n_rows << " x " << a.n_cols << "\n";
};

/* Function that does the stochastic EM per sample.
   Meant to be run by a single thread.             */
void thread_me(const int tid,
               const int subject,
               int nsamp_floor,
               const arma::mat  &proportions,
               const arma::mat  &preassign,
               const arma::vec  &doNotSample,
               int nsubsets,
               double iterAssignCoef,
               bool zeroPosteriorProbs,
               double prior,
               int intSampSize,
               bool betaDispersion,
               int keepEach,
               bool sampleRandom,
               int mat_size,
               // Progress  &prog,
               const arma::mat &subpopInd,
               const arma::mat &covariance,
               const arma::mat  &MHcoef,
               const arma::mat  &M,
               const arma::mat  &nullEta,
               const arma::vec  &index,
               const arma::vec  &condvar,
               const arma::mat  &altEta,
               const arma::mat  &Y,
               const arma::mat &N,
               const arma::mat &isingCoefs,
               const arma::mat &invcov,
               const arma::mat &rand,
               arma::mat  &clusterassignments,
               arma::cube  &assignmentMats,
               arma::cube  &randomEffectMats,
               arma::mat  &MHsuccessrates,
               std::mutex  &mut,
               std::vector<bool>& threadIsFinished) {
#ifdef PROFILE
  std::clock_t start;
  double duration;
#endif
  arma::vec this_success(nsubsets);
  arma::mat this_clusterDensity(intSampSize, 2);
  mut.lock();
  // copy shared data to local variables
  arma::mat thisclusterassignment(nsubsets, 1);
  std::copy(clusterassignments.begin_col(subject),
            clusterassignments.end_col(subject),
            thisclusterassignment.begin_col(0));

  arma::mat this_assignmentMats(mat_size, nsubsets);
  std::copy(assignmentMats.slice(subject).begin(),
            assignmentMats.slice(subject).end(),
            this_assignmentMats.begin());
  arma::mat this_randomEffects(mat_size, nsubsets);
  std::copy(randomEffectMats.slice(subject).begin(),
            randomEffectMats.slice(subject).end(),
            this_randomEffects.begin());
  mut.unlock();

  int abort = 0;
  int unifPosition = 0;
  double isingOffset = 0;
  int assignNum = 0;
  arma::uvec subject_indicator(1);
  subject_indicator(0) = subject;
  int sample = 0, subset = 0;
  arma::vec vsample(intSampSize);
#ifdef PROFILE
  start = std::clock();
#endif
  for (sample = 0; sample < nsamp_floor; sample++) {
    for (subset = 0; subset < nsubsets; subset++) {
      if (!abort) {
        if (doNotSample.at(subset) == 1.0 ||
            (preassign.at(subset, subject) != -1 &
             iterAssignCoef < 10e-4 &
             !zeroPosteriorProbs)) {
          thisclusterassignment.at(subset, 0) = 0;
          continue;
        } else if (preassign.at(subset, subject) != -1 &
                   iterAssignCoef < 10e-4 &
                   !zeroPosteriorProbs) {
          thisclusterassignment.at(subset, 0) = preassign.at(subset, subject);
          continue;
        } else if (preassign.at(subset, subject) == 0) {
          isingOffset = -prior;
        } else if (preassign.at(subset, subject) == 1) {
          isingOffset = prior;
        } else {
          isingOffset = 0;
        }
        arma::uvec subset_indicator = flowReMix::find(
            subpopInd.col(subject), subset + 1);
        if (subset_indicator.size() == 0) {
          continue;
        }
        arma::mat subsetNullEta, subsetAltEta,
            subsetN, subsetProp, subsetCount;
        subsetNullEta = nullEta(subset_indicator, subject_indicator);
        subsetAltEta = altEta(subset_indicator, subject_indicator);
        subsetN = N(subset_indicator, subject_indicator);
        subsetProp = proportions(subset_indicator, subject_indicator);
        subsetCount = Y(subset_indicator, subject_indicator);

        // Some cell populations for a specific subject may be missing
        double sigmaHat = sqrt(covariance.at(subset, subset));
        arma::vec empEta = logit_arma(
            flowReMix::pmax(subsetProp, 1.0 / subsetN));

        // integrating densities

        double mcoef;

        mcoef = MHcoef.at(subset);

        mcoef = std::max(1.0, mcoef);
        flowReMix::myrnorm3(vsample,
                            0,
                            sigmaHat * mcoef,
                            tid);
        double prevMuHat;
        for (int cluster = 0; cluster < 2; cluster++) {
          arma::mat eta;
          if (cluster == 0) {
            eta = subsetNullEta;
          } else {
            eta = subsetAltEta;
          }
          arma::vec importanceWeights;
          arma::vec sampNormDens = flowReMix::dnorm4(vsample,
                                                     0,
                                                     sigmaHat * mcoef,
                                                     true);
          arma::vec normDens = flowReMix::dnorm4(vsample,
                                                 0,
                                                 sigmaHat,
                                                 true);
          // presumably we're tuning mcoef.
          importanceWeights = normDens - sampNormDens;

          arma::mat randomEta;
          randomEta = computeRandomEta_arma(eta, vsample);
          double disp;
          disp = M.at(subset);


          int sampSize = randomEta.n_cols;
          int subsetSize = subsetCount.size();
          int count, n;
          double density, prob;

          arma::vec binomDensity(sampSize);
          int i = 0, j = 0;
          for (i = 0; i < sampSize; i++) {
            density = 0;
            for (j = 0; j < subsetSize; j++) {
              prob = randomEta.at(j, i);
              prob = expit(prob);
              count = subsetCount[j];
              n = subsetN[j];
              if (betaDispersion & (disp < 150000)) {
                density += betaBinomDens(count, n , prob, disp);
              } else {
                density += R::dbinom(count, n, prob, 1) / subsetSize;
              }
            }
            binomDensity.at(i) = density;
          }

          this_clusterDensity.col(cluster) =
              binomDensity +
              importanceWeights;
        }

        arma::vec integratedDensities;
        integratedDensities = computeIntegratedDensities_arma(
            this_clusterDensity);


        double priorProb = 0.5;
        if (sample >= 0) {
          thisclusterassignment.at(subset, 0) = 1;

          {
            arma::rowvec corow(nsubsets);
            arma::rowvec clust(nsubsets);
            clust = thisclusterassignment.col(0).t();
            corow = isingCoefs.row(subset);



            priorProb = expit(sum(corow % clust) + isingOffset);
          }
        } else {
          priorProb = 0.5;
        }

        double pResponder;
        double densityRatio;


        densityRatio = integratedDensities.at(0) / integratedDensities.at(1) *
            (1.0 - priorProb) / priorProb;
        pResponder = 1.0 / (1.0 + densityRatio);


        if (preassign.at(subset, subject) == 1) {
          pResponder = 1 - (1 - pResponder) * iterAssignCoef;
        } else if (preassign.at(subset, subject) == 0) {
          pResponder = pResponder * iterAssignCoef;
        }

        // if (unifVec.at(unifPosition++, subject) < pResponder) {
        if (R::runif(0, 1) < pResponder) {
          thisclusterassignment.at(subset, 0) = 1;
        } else {
          thisclusterassignment.at(subset, 0) = 0;
        }
      }
    }
    if ((sample % keepEach) == 0) {
      for (int i = 0; i < thisclusterassignment.col(0).n_elem; i ++) {
        this_assignmentMats.at(assignNum, i) =
            thisclusterassignment.at(i, 0);
      }
      assignNum++;
    }
  }

#ifdef PROFILE
    duration = (std::clock() - start) / static_cast<double>(CLOCKS_PER_SEC);
  mut.lock();
  /////////////////////////////////////////////////////////////////
  // number of loop iterations is nsubsets * nsubsets: inner and //
  // outer loop, times integral size times posterior samples,    //
  // times 2 for the two cluster components. This is really      //
  // inefficient. There must be something we can do              //
  // to improve this.                                            //
  /////////////////////////////////////////////////////////////////
  timings.push_back(static_cast<double>(nsamp_floor *
                                        nsubsets *
                                        nsubsets * 2 *
                                        intSampSize) /
                    duration);
  mut.unlock();
#endif
  // arma::vec unifVec2;
  // unifVec2 = flowReMix::myrunif(nsamp_floor * nsubsets, tid);


  arma::vec newEta(altEta.n_rows);
  newEta = nullEta.col(subject);   // should copy
  arma::rowvec assignment;
  assignment = this_assignmentMats.row(this_assignmentMats.n_rows-1);
  std::vector<int> responderSubset;
  responderSubset = flowReMix::match2(subpopInd.col(subject),
                                      flowReMix::which2(assignment));
  flowReMix::mapByVec(newEta, altEta.col(subject), responderSubset);
  arma::vec randomEst;
  randomEst = rand.col(subject);   // randomEst will  be modified?
  arma::vec MHattempts(nsubsets, arma::fill::zeros);
  arma::vec MHsuccess(nsubsets, arma::fill::zeros);


  if (sampleRandom) {
    double index_subject = index(subject);

    simRandomEffectCoordinateMH_mc(this_randomEffects,
                                   Y.col(subject),
                                   N.col(subject), index_subject,
                                   nsamp_floor, nsubsets, MHcoef,
                                   assignment.t(),
                                   subpopInd.col(subject), newEta,
                                   randomEst, condvar,
                                   covariance, invcov,
                                   MHattempts, MHsuccess,
                                   M, betaDispersion,
                                   keepEach, mat_size);
  }
  this_success = MHsuccess/MHattempts;
  mut.lock();
  std::copy(this_success.begin(),
            this_success.end(),
            MHsuccessrates.begin_col(subject));
  std::copy(this_randomEffects.begin(),
            this_randomEffects.end(),
            randomEffectMats.slice(subject).begin());
  std::copy(this_assignmentMats.begin(),
            this_assignmentMats.end(),
            assignmentMats.slice(subject).begin());
  std::copy(thisclusterassignment.begin_col(0),
            thisclusterassignment.end_col(0),
            clusterassignments.begin_col(subject));
  threadIsFinished[tid] = true;
  mut.unlock();
  return;
}

// [[Rcpp::export]]
List CppFlowSstepList_mc_vec(const int nsubjects, const arma::mat& Y,
                             const arma::mat& N, const arma::mat& subpopInd,
                             arma::mat& clusterassignments,
                             const arma::mat& nullEta, const arma::mat altEta,
                             const arma::mat& rand, arma::vec& index,
                             const arma::mat& preassign,
                             const int nsamp,
                             const int nsubsets,
                             const int intSampSize ,
                             const arma::mat& isingCoefs,
                             const arma::mat& covariance, const int keepEach,
                             const arma::vec& MHcoef, const bool betaDispersion,
                             const double randomAssignProb,
                             const double iterAssignCoef, const double prior,
                             const bool zeroPosteriorProbs,
                             const arma::vec& M, const arma::mat& invcov,
                             const bool mixed,
                             const bool sampleRandom,
                             const arma::vec& doNotSample,
                             const bool markovChainEM, int cpus, int seed) {
  int nsamp_floor = floor(nsamp / keepEach) * keepEach;
  int mat_size = static_cast<int> (nsamp_floor/keepEach);
  try {
    if (nsubjects != Y.n_cols) {
      throw std::domain_error("nsubjects and dimensions of Y don't match!");
    }
    if (nsubsets != rand.n_rows) {
      throw std::domain_error("nSubsets and dimensions of rand don't match!");
    }

    // allocate objects for return values; //
    arma::cube assignmentMats(mat_size, nsubsets, nsubjects);
    arma::cube randomEffectMats(mat_size, nsubsets, nsubjects);
    arma::mat MHsuccessrates(nsubsets, nsubjects);
    arma::mat proportions = Y / N;
    arma::vec condvar = (1.0) /(invcov.diag());

    if (mixed) {
      std::fill(assignmentMats.begin(), assignmentMats.end(), 1.0);
    }

    // Gibbs for each subject
    // arma::cube clusterDensities(intSampSize, 2, nsubjects);
    arma::mat iterPosteriors(nsubsets, nsubjects);

    auto max_threads = std::thread::hardware_concurrency();

    ////////////////////////////////////////////
    // reinitialize the parallel generator    //
    // only if the requested number           //
    //  of cpus is greater or less than       //
    // the size of the generator vector       //
    ////////////////////////////////////////////
    if (!ParallelNormalGenerator::isinit(cpus)) {
      ParallelNormalGenerator::initialize(cpus, seed);
    }
    if (!ParallelUnifGenerator::isinit(cpus)) {
      ParallelUnifGenerator::initialize(cpus, seed);
    }
    std::vector<std::thread> thread_vector;
    std::mutex mut;

    auto subject = 0;
    int tid = 0;
    if (cpus > nsubjects)
      cpus = nsubjects;
    Progress prog(nsubjects, true);
    std::vector<bool> threadIsFinished(cpus);
    std::fill(threadIsFinished.begin(),
              threadIsFinished.end(),
              true);   // no threads have started any work yet.

    // Start the initial set of workers.
    for (int tid = 0; tid < cpus; tid++) {
        thread_vector.push_back(std::thread(thread_me,
                                            tid,
                                            subject,
                                            nsamp_floor,
                                            std::ref(proportions),
                                            std::ref(preassign),
                                            std::ref(doNotSample),
                                            nsubsets,
                                            iterAssignCoef,
                                            zeroPosteriorProbs,
                                            prior,
                                            intSampSize,
                                            betaDispersion,
                                            keepEach,
                                            sampleRandom,
                                            mat_size,
                                            std::ref(subpopInd),
                                            std::ref(covariance),
                                            std::ref(MHcoef),
                                            std::ref(M),
                                            std::ref(nullEta),
                                            std::ref(index),
                                            std::ref(condvar),
                                            std::ref(altEta),
                                            std::ref(Y),
                                            std::ref(N),
                                            std::ref(isingCoefs),
                                            std::ref(invcov),
                                            std::ref(rand),
                                            std::ref(clusterassignments),
                                            std::ref(assignmentMats),
                                            std::ref(randomEffectMats),
                                            std::ref(MHsuccessrates),
                                            std::ref(mut),
                                            std::ref(threadIsFinished)));
        subject++;
#ifdef DEBUG
        std::cout <<
            "Initial launch of thread for subject "
                  << subject << " of "
                  << nsubjects << "\n";
#endif
    }

    // monitor the threads for completion.
    while (subject < nsubjects) {
      for (int tid = 0; tid < cpus; tid++) {
        if (threadIsFinished[tid] & thread_vector[tid].joinable()) {
            thread_vector[tid].join();
            prog.increment();
#ifdef DEBUG
            std::cout << "joined thread " << tid <<"\n";
#endif
            if (subject <  nsubjects) {
              mut.lock();
              threadIsFinished[tid] = false;
              mut.unlock();
              thread_vector[tid] = std::thread(thread_me,
                                               tid,
                                               subject,
                                               nsamp_floor,
                                               std::ref(proportions),
                                               std::ref(preassign),
                                               std::ref(doNotSample),
                                               nsubsets,
                                               iterAssignCoef,
                                               zeroPosteriorProbs,
                                               prior,
                                               intSampSize,
                                               betaDispersion,
                                               keepEach,
                                               sampleRandom,
                                               mat_size,
                                               std::ref(subpopInd),
                                               std::ref(covariance),
                                               std::ref(MHcoef),
                                               std::ref(M),
                                               std::ref(nullEta),
                                               std::ref(index),
                                               std::ref(condvar),
                                               std::ref(altEta),
                                               std::ref(Y),
                                               std::ref(N),
                                               std::ref(isingCoefs),
                                               std::ref(invcov),
                                               std::ref(rand),
                                               std::ref(clusterassignments),
                                               std::ref(assignmentMats),
                                               std::ref(randomEffectMats),
                                               std::ref(MHsuccessrates),
                                               std::ref(mut),
                                               std::ref(threadIsFinished));
#ifdef DEBUG
              std::cout << "Launched new thread " <<
                  tid << " for subject " <<
                  subject << "/" << nsubjects << "\n";
#endif
              subject++;
            } else {
              continue;
            }
        }  // end of if clause
      }  // end for loop over threads
    }  // end while loop over subjects
      // Join any dangling threads
    int alldone = 0;
    int ndone = 0;
    while (alldone != 2) {
      for (tid = 0; tid < thread_vector.size(); tid++) {
        if (thread_vector[tid].joinable()&threadIsFinished[tid]) {
          thread_vector[tid].join();
          prog.increment();
        }
      }  // end for loop
      ndone = std::accumulate(threadIsFinished.begin(),
                  threadIsFinished.end(), 0);
      if (ndone == thread_vector.size()) {
        alldone++;  // loop through twice to ensure we
                    // join all completed threads.
      }
    }
#ifdef PROFILE
    double sum = 0;
    double size =  timings.size();
    std::transform(timings.begin(), timings.end(),
                   timings.begin(),
                   [size] (double d) {
                     return d / size;
                   });
    sum = std::accumulate(timings.begin(), timings.end(), 0.0);
    std::cout << sum  <<
        " loop iterations per second per thread (avg) \n";
      timings.clear();
#endif
    List retval = List::create(
        Named("assign") = wrap(assignmentMats),
        Named("rand") = wrap(randomEffectMats),
        Named("rate") = wrap(MHsuccessrates));
    return retval;
  }  // end try block
  catch(std::exception const &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rf_error("C++ exception (unknown reason)");
  }
  return(List());
}
