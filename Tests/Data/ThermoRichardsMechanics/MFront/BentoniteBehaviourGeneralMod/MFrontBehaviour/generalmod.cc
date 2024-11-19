/*
    Copyright (c) since 2003, David Masin (masin@natur.cuni.cz)

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/****************** Hypoplasticity for unsaturated soils, expansive soils with *******************/
/****************** double-porosity structure, thermal effects ***********************************/

#include <string.h>
#include "generalmod.h"

using namespace std;

/**
  The function calculates increments of net stress, Sr and state variables dqstatev based on their
  initial values and increments of strains, pore fluid pressures and temperature. Sign convention: compression negative for stresses, strains and pressures.

  int c_nparms c_nstatev, ngstrain, c_ngstress jsou proměnné struktury General_model inicializované automaticky při jejím vytvoření.

  @param strain_gen[ngstrain] - (input): controlled variables [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]^T, updated to new values.
  @param stress_gen[c_ngstress] - (input): implied variables, updated to new values [signet11, signet22, signet33, signet12, signet13, signet23, Sr]^T. "signet" is net stress, that is total stress minus krondelta x ua.
  @param qstatev[c_nstatev] - (input): state variable vector, updated to new values.
  @param dstrain_gen[ngstrain] - (input): increment of controlled variables [deps11, deps22, deps33, dgamma12, dgamma13, dgamma23, ds, dT]^T
  @param dtime - (input): increment of time
  @param DDtan_gen[c_ngstress][ngstrain] - (output): generalised tangent stiffness d(dstress_gen)/d(dstrain_gen)
  @param parms[c_nparms] - (input): model parameters
  @param rkf_statev[nrkf_statev] - (input): array of state variables from SIFEL for Runge-Kutta, updated to new values.
  @param rkf_parms[nrkf_parms] - (input): array of parameters from SIFEL for Runge-Kutta
  @param flag - (input): 0 for state variable initialisation, 1 for dstress_gen & DDtan_gen, 2 for dstress_gen only, 3 for DDtan_gen only

  Order of statev variables in the arrays qstatev and dqstatev
  qstatev[0] = e (required on input)
  qstatev[1] = suction
  qstatev[2] = Sr
  qstatev[3] = T
  qstatev[4] = em
  qstatev[5] = eM
  qstatev[6] = SrM
  qstatev[7] = a_scan (required on input)
  qstatev[8] = re

  Order of parameters in the array rkf_parms
  rkf_parms[0] = err_sig
  rkf_parms[1] = h_min
  rkf_parms[2] = ni
  rkf_parms[3] = sv_rkf_zero
  rkf_parms[4] = rkt

  Order of parameters in the array rkf_statev
  rkf_statev[0] = neval
  rkf_statev[1] = nstep
  rkf_statev[2] = mindt
  rkf_statev[3] = dtsub

  @return The function does not return anything.

  Masin[1] : Double structure hydromechanical coupling formalism and a model for unsaturated expansive clays, Engineering Geology 165 (2013) 73–88
  Masin[2] : Coupled thermo-hydro-mechanical double structure model for expansive soils, ASCE Journal of Engineering Mechanics,
             special issue on "Mechanics of Unsaturated Porous Media", 2016


  Created by D. Masin 13.7.2015, modified by G. Scaringi, February 2019
*/

long Hypoplasti_unsat_expansive_thermal::soil_model(double strain_gen[], double stress_gen[],
  double qstatev[], double dstrain_gen[], double dtime, double *DDtan_gen, double parms[], double rkf_statev[], double rkf_parms[], int flag, int kinc/*=0*/, int ipp/*=0*/) {

  long errorrkf = 0;
  long errormath = 0;

  //Remember the original values for perturbation
  double stress_gen_init[max_ngstress];
  double strain_gen_init[max_ngstrain];
  double qstatev_init[max_nstatev];
  garray_set(stress_gen_init, 0.0, max_ngstress);
  garray_set(strain_gen_init, 0.0, max_ngstrain);
  garray_set(qstatev_init, 0.0, max_nstatev);
  garray_move(stress_gen, stress_gen_init, ngstress);
  garray_move(strain_gen, strain_gen_init, ngstrain);
  garray_move(qstatev, qstatev_init, nstatev);

  if (parms)
    initialise_parameters(parms);
  if (flag < 0)
    return 0;

  //Initialises state variables and parameters (flag 0)
  if (flag == 0) {
    initialise_variables(strain_gen, stress_gen, qstatev, dstrain_gen, 1.0, parms, kinc);
  }

  //direct stiffness matrix (approximate, but can be calculated with zero strain increment provided)
  else if (flag == 3) {
    double dstrain_gen_zero[max_ngstrain];
    garray_set(dstrain_gen_zero, 0, max_ngstrain);
    errormath += correct_statev_values(strain_gen, stress_gen, qstatev, dstrain_gen_zero, 1);

    direct_stiffness_matrix(strain_gen, stress_gen, qstatev, dstrain_gen, 1.0, DDtan_gen, parms, kinc, flag);
  }

  //Updates variables, if flag == 1 also stiffness matrix by perturbation
  else if (flag == 1 || flag == 2) {

    errormath += correct_statev_values(strain_gen, stress_gen, qstatev, dstrain_gen, 1);

    //Update variables
    errorrkf += calc_rkf(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
    errormath += correct_statev_values(strain_gen, stress_gen, qstatev, dstrain_gen, 2);

    //Calculate stiffness matrix by perturbation of current strain increment
    if(flag == 1) {
      direct_stiffness_matrix(strain_gen_init, stress_gen_init, qstatev_init, dstrain_gen, dtime, DDtan_gen, parms, kinc, flag);
    }
  }
  else {
    cout<<"Generalmod called with wrong flag"<<endl;
  }

  long matherror_here=check_math_error();
  if(matherror_here || errormath) {
    cout<<"Generalmodel ended up with math error, errormath = "<<errormath<<", matherror_here = "<<matherror_here<<", errorrkf = "<<errorrkf<<", flag = "<<flag<<endl;
    return(1);
  }

  if(errorrkf) {
    garray_move(stress_gen_init, stress_gen, ngstress);
    garray_move(strain_gen_init, strain_gen, ngstrain);
    garray_move(qstatev_init, qstatev, nstatev);
  }
  return (0);
}

/**
  The function calculates increment of net stress vector, of state variable vector and new Sr.

  @param signet[9] - (input): total stress: old value. This is not to be updated, increment of stress is in dsignet[9].
  @param suction - (input): suction: old value. This is not to be updated.
  @param Temper - (input): temperature: old value. This is not to be updated.
  @param qstatev[] - (input): state variable vector: old values. This is only to be updated if flag==0, otherwise only their increment is produced in dqstatev[].
  @param deps[9] - (input): increment of total strain
  @param dsuction - (input): increment of pore air pressure
  @param dTemper - (input): increment of temperature
  @param dsignet[9] - (output): increment of net stress
  @param Sr - (output): new degree of saturation
  @param dqstatev[] - (output): increment of state variables
  @param kinc       - (input):

  @return The function does not return anything.

  Created by D. Masin 13.7.2015
*/
bool Hypoplasti_unsat_expansive_thermal::calc_fsigq(double signet[9], double suction, double Temper,
  double *qstatev, double deps[9], double dsuction, double dTemper, double dtime,
  double dsignet[9], double &Sr_new, double dqstatev[], int kinc) {

  long errorfsigq = 0;

  //BEGIN convert increments to rates
  dsuction/=dtime;
  dTemper/=dtime;
  garray_multiply(deps, deps, 1/dtime, 9);
  //END convert increments to rates

  //calculation of tangent operators L, N, Hs, HT
  double kron_delta[3][3];
  garray_set(&kron_delta[0][0],0,9);
  kron_delta[0][0]=kron_delta[1][1]=kron_delta[2][2]=1;                             // initialises kroneker delta matrix
  double hypo_L[3][3][3][3];
  double hypo_N[3][3];
  double H_unsat[3][3];
  double HT_unsat[3][3];
  garray_set(&hypo_L[0][0][0][0],0,81);
  garray_set(&hypo_N[0][0],0,9);
  garray_set(&H_unsat[0][0],0,9);
  garray_set(&HT_unsat[0][0],0,9);                                                  // initialises tensors
  double fm=0.5;                                        // initially sets fm double-structure coupling factor to 0.5
  double fm_lambda_act=0.5;
  double fu=0;                                          // initially sets fu factor (to Hs and HT) to zero
  double SrM=qstatev[6];                                // takes the macrostructural Sr from the state variables vector
  double Sr=qstatev[2];
  bool swelling=false;

  double eM=qstatev[5];                                 // takes the macrostructural void ratio from the state variables vector
  //double SrM=qstatev[6];                              // does NOT take Sr again
  double a_scan=qstatev[7];                             // takes scanning curve
  double re=qstatev[8];                                 // takes re (relative void ratio)
  double sewM=aer*sairentry0*eM0/eM*(aTparam+bTparam*Temper)/(aTparam+bTparam*Trparam);

  give_LNH(signet, suction, Temper, qstatev, hypo_L, hypo_N, H_unsat, HT_unsat, swelling, fu); // calculates tangent operators

  if(re<0) {
    if(debug) cout<<"re<0 for some reason, correcting it, just for fm calculation. re="<<re<<" eM="<<eM<<" em="<<qstatev[4]<<" e="<<qstatev[0]<<" ascan="<<a_scan<<endl;
    re=0;
  }

  //Macrostructural strain increment for the given iteration
  double depM[9];
  garray_move(deps,depM,9);                                     // moves deps into depM (newly created)
  //Macrostructural stress increment for the given iteration
  double dsigMef[9];
  garray_set(dsigMef,0,9);
  //Total stress increment for the given iteration
  double dsignet_iter[9];
  garray_set(dsignet_iter,0,9);
  //Iteration counters
  int num_iter=0;                                               // initialises the iteration count to zero
  double error_measure=1.e10;                                   // sets the error measure to 1e10 (will be changed anyway)
  //Numerical parameter controlling step size for trial-and-error search for depm
  /*
  double step_depm=garray_size(&deps[0], 9);
  double minstep=1.e-5;
  double mexerror=step_depm/500;
  mexerror=round_to_digits(mexerror,1);
  if(mexerror<5.e-10) mexerror=5.e-10;
  bool tried=false;
  if(step_depm<minstep) step_depm=minstep;
  */

  double mexerror=garray_size(&deps[0], 9)/1000;
  mexerror=round_to_digits(mexerror,1);
  if(mexerror<1.e-10)
    mexerror=1.e-10;

  //Trial trace of dep_m for the given iteration
  double trdepm=0;                                          // here trdepm is initialised to zero, and this value is used in the first iteration

  double rlambda_for_ascan=0;
  double rlambda_for_se=0;
  double rlambda=give_rlambda(suction, a_scan, dsuction, sewM, SrM, rlambda_for_ascan, rlambda_for_se);

  double iterfactor=0;

  while(1) {

  //calculation of the double structure coupling factor
    swelling=true;

    bool notbytrdepm=false;
    if(fabs(dsuction)>1.e-4) notbytrdepm=true;

    if(!notbytrdepm && trdepm<0) swelling=false;
    else if (notbytrdepm && dsuction<0) swelling=false;


    if(num_iter == 0 || num_iter == 1 || num_iter == 5 || num_iter % 50 == 0)
        give_LNH(signet, suction, Temper, qstatev, hypo_L, hypo_N, H_unsat, HT_unsat, swelling, fu);

    give_fm(fm, fm_lambda_act, suction, sewM, re, swelling);
    double depm[9];
    garray_set(depm,0,9);

    for(int i=0; i<3; i++) depM[3*i+i]=deps[3*i+i]-fm*trdepm/3;

    double DM[3][3];
    garray_set(&DM[0][0], 0, 3*3);
    garray_move(depM, &DM[0][0], 3*3);
    double normDM=garray_size(&DM[0][0], 3*3);

    double linear[3][3];
    garray_set(&linear[0][0], 0, 3*3);
    gmatrix_a4b(hypo_L, *DM, *linear);

  //Calculate effective stress rate of Macrostructure
    for(int i=0; i<3; i++) {
      for(int j=0; j<3; j++) {
        dsigMef[3*i+j]= linear[i][j] + hypo_N[i][j] * normDM;
        if(dsuction>0)
          dsigMef[3*i+j]+=rlambda*H_unsat[i][j]*dsuction;
        if(dTemper>0)
          dsigMef[3*i+j]+=HT_unsat[i][j]*dTemper;
      }
    }

  //Calculate total stress rate from the effective stress rate of Macrostructure
    double gamma=lambdap0;
    double corr_s=(1-gamma*rlambda)*dsuction;
    double deM=(1+eM)*(gtrace(DM)+(fm-1)*trdepm);
    double corr_eM=-rlambda_for_se*gamma*suction*deM/eM;
    double corr_T=+rlambda_for_se*gamma*suction*bTparam*dTemper/(aTparam+bTparam*Temper);
    double chi=SrM;

    for(int i=0; i<3; i++) {
      for(int j=0; j<3; j++) {
         dsignet_iter[3*i+j]=dsigMef[3*i+j]-kron_delta[i][j]*chi*(corr_s+corr_eM+corr_T);
      }
    }
    double em=qstatev[4];

    if(num_iter>1000) {//for zero eM all the deformation is due to microstructure
      fu=0;
      garray_set(depM,0,9);

      trdepm=(deps[0]+deps[4]+deps[8]);
      double pefsat=(signet[0]+signet[4]+signet[8])/3+suction;
      double demdpm=give_demdpm(pefsat, em);
      double dpefsat=(trdepm-alpha_s*dTemper)*(1+em)/demdpm;
      dsignet_iter[0]=dsignet_iter[4]=dsignet_iter[8]=(dpefsat-dsuction);
    }

  //Calculate microstructural strain corresponding to the obtained dsignet in the current iteration.
    double trdepm_iter = calc_trdepm(signet, suction, Temper, em, dsuction, dTemper, dsignet_iter, fu);

  //Calculate the difference between microstructural strain obtained in the current iteration and the trial strain trdepm
    double old_error_measure=error_measure;
    error_measure=trdepm-trdepm_iter;

  //Solution found
    if(gscalar_dabs(error_measure)<mexerror)
      break;                                                // stops the iteration if the error is small

  //Solution not found, update trdepm
    /*
    if(gscalar_dabs(error_measure)>gscalar_dabs(old_error_measure)) {
        if((error_measure<0 && old_error_measure<0) || (error_measure>0 && old_error_measure>0))
            if(num_iter!=0) step_depm*=-0.5;
        else step_depm*=-1;
    }
    else if((old_error_measure>0 && error_measure<0) || (old_error_measure<0 && error_measure>0)) {
        if(num_iter!=0) step_depm*=-0.5;
    }
    double minsteprun=1.e-10;
    if(step_depm>0 && step_depm<minsteprun) step_depm=minsteprun;
    else if(step_depm<0 && step_depm>-1*minsteprun) step_depm=-1*minsteprun;
    if(num_iter>800 && !tried) {
         cout<<"Microstructural strain: critical number of iterations, err="<<error_measure<<endl;
         step_depm*=1.e6;
         tried=true;
    }
    trdepm+=step_depm;
    */
    if((old_error_measure>0 && error_measure<0) || (old_error_measure<0 && error_measure>0)) iterfactor+=1;

    trdepm=(trdepm_iter+iterfactor*trdepm)/(iterfactor+1);                       // after the first iteration trdepm is updated from zero to 1/4 of 3*depm_iter[0] (+3/4 of trdepm later when trdepm is not zero)

    if(num_iter>1000) {
       errorfsigq = true;
       if(debug) {
         cout<<endl<<"-----------------------------------------------------"<<endl;
         cout<<"Error in microstructural iterations."<<endl;
         for(int i=0; i<6; i++) cout<<"signet["<<i<<"] "<<signet[i]<<endl;
         cout<<"suction "<<suction<<endl;
         for(int i=0; i<nstatev; i++) cout<<"qstatev["<<i<<"] "<<qstatev[i]<<endl;
         cout<<"-----------------------------------------------------"<<endl;
       }
       exit(1);

    }
    num_iter++;
  }
  garray_move(dsignet_iter, dsignet, 9);

  //BEGIN rates back to increments
  dsuction*=dtime;
  dTemper*=dtime;
  garray_multiply(deps, deps, dtime, 9);
  garray_multiply(dsignet, dsignet, dtime, 9);
  //END rates back to increments

  //this calculates dqstatev, it is not updating the original values
  int flag=2;
  update_hisv(signet, suction, Temper, qstatev, deps, dsuction, dTemper, dsignet, Sr_new, dqstatev, flag, fu, kinc);

  return(errorfsigq);
}

void Hypoplasti_unsat_expansive_thermal_2017::give_fm(double& fm, double& fm_lambda_act, double suction, double sewM, double re, bool swelling) {

  double fm_swelling=1-gscalar_power(re,mmparam); //swelling
  fm=(csh*(suction/sewM));

  if(swelling) fm=fm_swelling;
  else if(fm>(1-fm_swelling)) fm=1-fm_swelling;

  if(fm<0) fm=0;
  else if (fm>1) fm=1;

  fm_lambda_act = -fm; //in the paper 2017, sign fm was wrong in lambda_act calculation
}

void Hypoplasti_unsat_expansive_thermal_2022::give_fm(double& fm, double& fm_lambda_act, double suction, double sewM, double re, bool swelling) {

  double fm_swelling=1-gscalar_power(re,mmparam); //swelling
  fm=(csh*(suction/sewM));

  if(swelling) fm=fm_swelling;
  else if(fm>(1-fm_swelling)) fm=1-fm_swelling;

  if(fm<0) fm=0;
  else if (fm>1) fm=1;

  fm_lambda_act = -fm; //in the paper 2017, sign fm was wrong in lambda_act calculation
}

void Hypoplasti_unsat_expansive_thermal::give_fm(double& fm, double& fm_lambda_act, double suction, double sewM, double re, bool swelling) {

  double fm_swelling=1-gscalar_power(re,mmparam);
  double fm_shrinkage=csh*(suction/sewM);
  if(suction>sewM) fm_shrinkage=csh;
  if (fm_shrinkage > (1-fm_swelling)) fm_shrinkage=1-fm_swelling;

  if(fm_swelling<0) fm_swelling=0;
  else if (fm_swelling>1) fm_swelling=1;
  if(fm_shrinkage<0) fm_shrinkage=0;
  else if (fm_shrinkage>1) fm_shrinkage=1;

  if(swelling) fm=fm_swelling;
  else fm=fm_shrinkage;

  fm_lambda_act=fm_shrinkage;
}


/*void Hypoplasti_unsat_expansive_thermal_2022::give_fm(double& fm, double& fm_lambda_act, double suction, double sewM, double re, double swelling) {

}*/

/**
  The function calculates state variable increment vector dqstatev[]. The old values of qstatev[0] (void ratio) and qstatev[6] (a_scan) must be provided by the FEM (even for the initialisation stage). The order of state variables is as follows:
  [0]e (required on input)
  [1]suction
  [2]Sr
  [3]T
  [4]em
  [5]eM
  [6]SrM
  [7]a_scan (required on input)
  [8]re
  [9] swelling indicator ???

  @param signet[9] - (input): total stress: old value. This is not to be updated, increment of stress is in dsignet[9].
  @param suction - (input): suction: old value. This is not to be updated.
  @param Temper - (input): temperature: old value. This is not to be updated.
  @param qstatev[] - (input): state variable vector: old values. This is only to be updated if flag==0, otherwise only their increment is produced in dqstatev[].
  @param deps[9] - (input): increment of total strain
  @param dsuction - (input): increment of pore air pressure
  @param dTemper - (input): increment of temperature
  @param dsignet[9] - (output): increment of total stress
  @param Sr - (output): new degree of saturation
  @param dqstatev[] - (output): increment of state variables
  @param flag -
  @param fu -

  @return The function does not return anything.

  Created by D. Masin 13.7.2015
*/
void Hypoplasti_unsat_expansive_thermal::update_hisv(double signet[9], double suction, double Temper,
  double qstatev[], double deps[9], double dsuction, double dTemper, double dsignet[9], double &Sr, double dqstatev[], int flag, double fu, int kinc) {

  double D[3][3];
  garray_set(&D[0][0], 0.0, 3*3);
  garray_move(deps, &D[0][0], 3*3);

  //old values of state variables needed for calculations
  double evoid=qstatev[0];
  double a_scan=qstatev[7];

  //calculating state variable increments/new values
  double delta_e=(1+evoid)*gtrace(D);
  double evoid_new=evoid+delta_e;
  if(evoid_new<0) evoid_new=0;
  double Temper_new=Temper+dTemper;

  double signet_new[9];
  garray_add(signet, dsignet, signet_new, 9);
  double suction_new=suction+dsuction;

  double em=qstatev[4];
  double em_new=em;

  if(flag==0) {
      double pefsat_new=(signet_new[0]+signet_new[4]+signet_new[8])/3+suction_new;
      double pefsat=(signet[0]+signet[4]+signet[8])/3.0+suction;
      double emstarT=exp(alpha_s*(Temper-Trparam))-1.0;
      double emstarT_new=exp(alpha_s*(Temper_new-Trparam))-1.0;
      double em_pefsat_new=init_em_pefsat(pefsat_new);
      double em_pefsat=init_em_pefsat(pefsat);
      em_new=exp(log(1+em_pefsat_new)+log(1+emstarT_new))-1.0;
      em=exp(log(1+em_pefsat)+log(1+emstarT))-1.0;
  }

  double trdepm = calc_trdepm(signet, suction, Temper, em, dsuction, dTemper, dsignet, fu);
  double d_em=(1+em)*trdepm;
  em_new=em+d_em;

  double min_void=0.0001;
  if(em_new<min_void) em_new=min_void;
  if(em<min_void) em=min_void;
  if(em_new > evoid_new-min_void) em_new=evoid_new-min_void;
  if(em>evoid-min_void) em=evoid-min_void;

  // old value of void ratio at macrolevel - equation (32) of Masin[2]
  double eM=(evoid-em)/(1.0+em);
  // new value of void ratio at macrolevel - equation (32) of Masin[2]
  double eM_new=(evoid_new-em_new)/(1.0+em_new);

  if(eM<0.0) {
    if(debug) cout<<"eM<0, this is not good, eM="<<eM<<", e="<<evoid<<", em="<<em<<", signet_new[0]="<<signet_new[0]<<", signet_new[4]="<<signet_new[4]<<", signet_new[8]="<<signet_new[8]<<", suction_new="<<suction_new<<", signet[0]="<<signet[0]<<", signet[4]="<<signet[4]<<", signet[8]="<<signet[8]<<", suction="<<suction<<endl;
    //return;
  }
  if(em<0.0) {
    if(debug) cout<<"em<0, this is not good, eM="<<eM<<", e="<<evoid<<", em="<<em<<endl;
    //return;
  }

  if (a_scan>1.001 || a_scan<-0.001) {
    if(debug) cout<<"something wrong with a_scan, a_scan=" << a_scan << endl;
    //return;
  }
  double eMuse=eM;
  double min_eM_use=0.1;
  if(eM<min_eM_use)
    eMuse=min_eM_use;
  double eM_newuse=eM_new;
  if(eM_new<min_eM_use)
    eM_newuse=min_eM_use;

  // suction at air entry value for drying process and original values -  equation (40) of Masin [2]
  double sedM=sairentry0*eM0/eMuse*(aTparam+bTparam*Temper)/(aTparam+bTparam*Trparam);
  // old value of suction at air expulsion value (Masin [1]) s_{exp}=a_e*s_{en} for which we can determine d_a_scan from equation (41) of Masin [2]
  double sewM=sedM*aer;
  // suction at air entry value for drying process and new values -  equation (40) of Masin [2]
  double sedM_new=sairentry0*eM0/eM_newuse*(aTparam+bTparam*Temper_new)/(aTparam+bTparam*Trparam);
  // new value of suction at air expulsion value (Masin [1]) s_{exp}=a_e*s_{en} for which we can determine d_a_scan from equation (41) of Masin [2]
  double sewM_new=sedM_new*aer;

  // old value of suction at air entry value - equation (37) of Masin[2]
  double se=sedM*(aer+a_scan-aer*a_scan);

  double a_scan_new=0.0;
  double SrM=qstatev[6];
  if (suction<sewM) {
    // suction at the main drying curve - equation (39) of Masin[2]
    double sD=sedM/se*suction;
    double sWM=sewM/(gscalar_power(SrM,(1/lambdap0)));
    double sDM=sWM/aer;
    double d_a_scan=0.0;
    double rlambda_for_ascan=0;
    double rlambda_for_se=0;
    double rlambda=give_rlambda(suction, a_scan, dsuction, sewM, SrM, rlambda_for_ascan, rlambda_for_se);

    if(aer<1 && suction<sWM && suction>sDM ) d_a_scan = (1.0-rlambda_for_ascan)/(sD*(1.0-aer))*(suction_new-suction);
    a_scan_new=a_scan+d_a_scan;
    if (a_scan_new>1.0 )
      a_scan_new=1.0;
    if (a_scan_new<0.0 )
      a_scan_new=0.0;
  }
  // new value of suction at air entry level - equation (37) of Masin[2]
  double se_new=sedM_new*(aer+a_scan_new-aer*a_scan_new);

  /* functions which distinguish original and updated versions */
  double SrM_new=give_SrM_new(suction, qstatev[6], eMuse, Temper, sewM, se, dsuction, eM_newuse-eMuse, dTemper, a_scan, flag);
  if(SrM_new>1) SrM_new=1.0;
  if(SrM_new<0) SrM_new=0;

  // new value of total degree of saturation - equation (33) of Masin[2]
  double Sr_new=1;
  if(evoid_new!=0) Sr_new=(SrM_new*(evoid_new-em_new)+em_new)/evoid_new;
  Sr=Sr_new;

  double tmpscal=0.0;
  double Nuse_new=0.0;
  double lambda_use_new=0.0;

  //fill-in state variable increments
  dqstatev[0]=evoid_new-qstatev[0];
  dqstatev[1]=suction_new-qstatev[1];
  dqstatev[2]=Sr_new-qstatev[2];
  dqstatev[3]=Temper_new-qstatev[3];
  dqstatev[4]=em_new-qstatev[4];
  dqstatev[5]=eM_newuse-qstatev[5];
  dqstatev[6]=SrM_new-SrM;
  dqstatev[7]=a_scan_new-qstatev[7];
  dqstatev[8]=0; //re: to be filled later

  //calculate new state variables, needed just for re
  double qstatev_new[c_nstatev];
  garray_add(qstatev, dqstatev, qstatev_new, c_nstatev);

  //calculate new re
  calc_N_lambda_dpeds(Nuse_new, lambda_use_new, tmpscal, tmpscal, signet_new, suction_new, Temper_new, qstatev_new);
  double sigef_new[9];
  make_sig_unsat(signet_new, suction_new, SrM_new, sigef_new, suction_new);
  double pmeanef_new=-(sigef_new[0]+sigef_new[4]+sigef_new[8])/3.0;

  double emaxisot_new=min_void;
  if(pmeanef_new<1.0) emaxisot_new=exp(Nuse_new)-1.0;
  else if(Nuse_new > lambda_use_new*log(pmeanef_new)) emaxisot_new=exp(Nuse_new-lambda_use_new*log(pmeanef_new))-1.0; // maximum void ratio for the state corresponding to NCL - equation (88) of Masin[2]

  if ( /*((evoid_new>emaxisot_new) && (evoid_new>evoid)) ||*/ ((evoid_new<min_void) && (evoid_new<evoid)) ) {
    dqstatev[0]=0.0; // rate of total void ratio
    dqstatev[4]=0.0; // rate of void ratio on microlevel
    dqstatev[5]=0.0; // rate of void ratio on macrolevel
  }
  double re_new=min_void;
  if(emaxisot_new>em_new) re_new=(evoid_new-em_new)/(emaxisot_new-em_new);

  if(re_new<min_void) re_new=min_void;
  else if(re_new>1) re_new=1;

  double re=qstatev[8];
  double dre=re_new-re;
  dqstatev[8]=dre;
}

double Hypoplasti_unsat_expansive_thermal::init_em_pefsat(double pefsat) {
  double em_pefsat = exp(kappam*log(smstar/pefsat)+log(1+emstar))-1;
  return em_pefsat;
}


/**
  The function calculates the macrostructural effective stress.

  @param sig[9] - (input): net stress.
  @param suction - (input): pore air pressure.
  @param SrM - (input): degree of saturation of macrostructure.
  @param tensor_unsat[9] - (output): effective stress of macrostructure
  @param scalar_unsat - (output): suction

  @return The function does not return anything.

  Created by D. Masin 13.7.2015
*/
void Hypoplasti_unsat_expansive_thermal::make_sig_unsat(double sig[3*3], double suction, double SrM,
  double tensor_unsat[3*3], double &scalar_unsat) {

  garray_set(tensor_unsat, 0, 3*3);
  garray_move(sig, tensor_unsat, 3*3);
  scalar_unsat=suction;

  double chi=SrM;

  tensor_unsat[0]+=chi*scalar_unsat-pt;
  tensor_unsat[4]+=chi*scalar_unsat-pt;
  tensor_unsat[8]+=chi*scalar_unsat-pt;
}

/**
  The function calculates the microstructural strain.

  @param signet[9] - (input): total stress.
  @param suction - (input): suction
  @param Temper - (input): temperature.
  @param dsuction - (input): increment of matric suction
  @param dTemper - (input): increment of temperature
  @param dsignet[9] - (input!!!): increment of total stress
  @param depm[9] - (output): microstructural strain

  @return The function does not return anything.

  Created by D. Masin 13.7.2015
*/
double Hypoplasti_unsat_expansive_thermal::calc_trdepm(double signet[9], double suction, double Temper, double em, double dsuction, double dTemper, double dsignet[9], double fu) {

  double pefsat=(signet[0]+signet[4]+signet[8])/3+suction;
  double dpefsat=(dsignet[0]+dsignet[4]+dsignet[8])/3+dsuction;
  double demdpm=give_demdpm(pefsat, em);
  double trdepm=demdpm*dpefsat/(1+em)+alpha_s*dTemper;

  if(fu>1) fu=1;
  if(trdepm>0) trdepm*=(1.0-fu);
  return trdepm;
}

double Hypoplasti_unsat_expansive_thermal::give_demdpm(double pefsat, double em) {
  if(pefsat>-1) pefsat=-1;
  double demdpm = -kappam*(1+em)/pefsat;
  return demdpm;
}


/**
  The function calculates current position and slope of normal compression line Nuse and lambda_use, and also the derivatives (dpe/ds)/pe and (dpe/dT)/pe

  @param Nuse - (output): Nuse
  @param lambda_use - (output): lambda_use
  @param dpeds_divpe - (output): (dpe/ds)/pe
  @param dpedT_divpe - (output): (dpe/dT)/pe
  @param signet[9] - (input): total stress.
  @param suction - (input): suction.
  @param Temper - (input): temperature.
  @param qstatev[] - (input): state variable vector.

  @return The function does not return anything.

  Created by D. Masin 13.7.2015
*/

void Hypoplasti_unsat_expansive_thermal::calc_N_lambda_dpeds(double &Nuse, double &lambda_use, double &dpeds_divpe,
  double &dpedT_divpe, double signet[9], double suction, double Temper, double qstatev[]) {

  double evoid=qstatev[0];
  double SrM=qstatev[6];
  if(SrM>1) SrM=1.0;
  if(SrM<0) SrM=0;

  // ln(s/s_e) where s/s_e has been expressed from equation (36) of Masin[2]
  double logsse=log(1/(gscalar_power(SrM,(1/lambdap0))));

  double sig_ef[9];
  make_sig_unsat(signet, suction, SrM, sig_ef, suction);

  // specific volume for reference pressure - equation (75a) of Masin[2]
  Nuse=N+nparam*logsse+nTparam*log(Temper/Trparam);
  // slope of the normal consolidation line - equation (75b) of Masin[2]
  lambda_use=lambda+lparam*logsse+lTparam*log(Temper/Trparam);

  // Hvorslev equivalent pressure - equation (74) of Masin[2]
  double pe=exp((Nuse-log(1+evoid))/lambda_use);
  dpeds_divpe=0;

  double eM=qstatev[5];                                 // takes the macrostructural void ratio from the state variables vector
  double sewM=aer*sairentry0*eM0/eM*(aTparam+bTparam*Temper)/(aTparam+bTparam*Trparam);

  if(suction<sewM) dpeds_divpe=(nparam-log(pe)*lparam)/(lambda_use*suction); // \frac{\rm d}p_e}/{{\rm d}s}\frac{1}{p_e} according to equation (74) of Masin[2]
  dpedT_divpe=(nTparam-log(pe)*lTparam)/(lambda_use*Temper);  // \frac{\rm d}p_e}/{{\rm d}T}\frac{1}{p_e} according to equation (74) of Masin[2]
};

/**
  The function calculates tangent operators hypo_L, hypo_N, H_unsat and HT_unsat of the model for macrostructure

  @param signet[9] - (input): total stress.
  @param suction - (input): suction.
  @param Temper - (input): temperature.
  @param qstatev[] - (input): state variable vector.
  @param hypo_L[3][3][3][3] - (output): tangent operator hypo_L
  @param hypo_N[3][3] - (output): tangent operator hypo_N
  @param H_unsat[3][3] - (output): tangent operator H_unsat
  @param HT_unsat[3][3] - (output): tangent operator HT_unsat

  @return The function does not return anything.

  Created by D. Masin 13.7.2015
*/
void Hypoplasti_unsat_expansive_thermal::give_LNH(double signet[9], double suction,
  double Temper, double qstatev[], double hypo_L[3][3][3][3], double hypo_N[3][3], double H_unsat[3][3],
  double HT_unsat[3][3], bool swelling, double &fu) {

  double sig_ef[9];
  double stresssuction=0;
  double SrM=qstatev[6];

  make_sig_unsat(signet, suction, SrM, sig_ef, stresssuction);
  double T[3][3];
  garray_set(&T[0][0], 0, 3*3);
  garray_move(sig_ef, &T[0][0], 3*3);
  double p_mean=-gtrace(T)/3;
  if(p_mean<1) {
    if(debug) cout<<"Correcting mean stress, p_mean="<<p_mean<<endl;
    p_mean=1;
  }

  double kron_delta[3][3];
  garray_set(&kron_delta[0][0],0,9);
  kron_delta[0][0]=kron_delta[1][1]=kron_delta[2][2]=1;

  double T_star[3][3];
  garray_move( &T[0][0], &T_star[0][0], 3*3 );
  for ( int idim=0; idim<3; idim++ ) T_star[idim][idim] += p_mean;
  double T_str_star[3][3];
  garray_multiply(&T_star[0][0], &T_str_star[0][0], (1/gtrace(T)), 3*3);
  double T_str[3][3];
  garray_multiply(&T[0][0], &T_str[0][0], (1/gtrace(T)), 3*3);
  double evoid=qstatev[0];

  double I1=gtrace(T);
  double I2=(garray_inproduct(&T[0][0], &T[0][0], 3*3)-I1*I1)/2;
  double I3=gmatrix_determinant(&T[0][0], 3);
  double I3plusI1I2=I3+I1*I2;
  if(I3plusI1I2==0) I3plusI1I2=1.e-10;

  double Nuse=N;
  double lambda_use=lambda;
  double dpeds_divpe=0;
  double dpedT_divpe=0;
  calc_N_lambda_dpeds(Nuse, lambda_use, dpeds_divpe, dpedT_divpe, signet, suction, Temper, qstatev);

  double peast=exp((Nuse-log(1+evoid))/lambda_use);
//        cout<<"OCR="<<peast/p_mean<<endl;
  double fd=(2*(p_mean-pt))/peast;
  if(fd<1.e-10) fd=1.e-10;
  double cos2phic=1-sin(phic)*sin(phic);
  double sin2phim=(9*I3+I1*I2)/I3plusI1I2;
  bool isotrop=false;
  if(sin2phim<1.e-10) {
    sin2phim=0;
    isotrop=true;
  }
  if(sin2phim>(1-1.e-10)) {
      //ret=2; //this was uncommented but the error is still there
      sin2phim=(1-1.e-10);
  }

  double ashape=0.3;
  double npow=-log(cos2phic)/log(2)+ashape*(sin2phim-sin(phic)*sin(phic));
  double phimfactor=-8*I3/I3plusI1I2;
  if(phimfactor<1.e-10)
    phimfactor=1.e-10;
  double fdsbs=2*gscalar_power(phimfactor,(1/npow));
  double a=sqrt(3.)*(3-sin(phic))/(2*sqrt(2.)*sin(phic));
  double alpha_power=log((lambda-kappa)*(3+a*a)/((lambda+kappa)*a*sqrt(3.)))/log(2.);
  fd=gscalar_power(fd,alpha_power);
  fdsbs=gscalar_power(fdsbs,alpha_power);
  if(fdsbs<fd && sin2phim>(1-1.e-9))
    fdsbs=fd;

  double fddivfdA=fd/fdsbs;
  double nuhh=nuparam;
  double nuvh=nuhh/alphanu;

  double Am=nuvh*nuvh*(4*alphaE*alphanu-2*alphaE*alphaE*alphanu*alphanu+2*alphaE*alphaE-alphanu*alphanu)+nuvh*(4*alphaE+2*alphaE*alphanu)+1+2*alphaE;
  double fs=9*p_mean/2*(1/kappa+1/lambda_use)/Am;

  double a1=alphaE*(1-alphanu*nuvh-2*alphaE*nuvh*nuvh);
  double a2=alphaE*nuvh*(alphanu+alphaE*nuvh);
  double a3=alphaE*nuvh*(1+alphanu*nuvh-alphanu-alphaE*nuvh);
  double a4=(1-alphanu*nuvh-2*alphaE*nuvh*nuvh)*(alphaE*(1-alphaG))/alphaG;
  double a5=alphaE*(1-alphaE*nuvh*nuvh)+1-alphanu*alphanu*nuvh*nuvh-2*alphaE*nuvh*(1+alphanu*nuvh)-2*alphaE*(1-alphanu*nuvh-2*alphaE*nuvh*nuvh)/alphaG;

  garray_set(&hypo_L[0][0][0][0], 0, 3*3*3*3);
  double nvect[3];
  garray_set(nvect,0,3);
  nvect[0]=1;//vertical direction is "0"
  double pmat[3][3];
  garray_set(&pmat[0][0],0,9);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      pmat[i][j]=nvect[i]*nvect[j];
    }
  }
  double kck[3][3][3][3];
  double kdk[3][3][3][3];
  double pdk[3][3][3][3];
  double kdp[3][3][3][3];
  double pck[3][3][3][3];
  double pdp[3][3][3][3];
  garray_set(&kck[0][0][0][0], 0, 81);
  garray_set(&kdk[0][0][0][0], 0, 81);
  garray_set(&pdk[0][0][0][0], 0, 81);
  garray_set(&kdp[0][0][0][0], 0, 81);
  garray_set(&pck[0][0][0][0], 0, 81);
  garray_set(&pdp[0][0][0][0], 0, 81);
  for ( int i=0; i<3; i++ ) {
    for ( int j=0; j<3; j++ ) {
      for ( int k=0; k<3; k++ ) {
        for ( int l=0; l<3; l++ ) {
          kck[i][j][k][l]=(kron_delta[i][k]*kron_delta[j][l]+kron_delta[i][l]*kron_delta[j][k]+kron_delta[j][l]*kron_delta[i][k]+kron_delta[j][k]*kron_delta[i][l])/2;
          kdk[i][j][k][l]=kron_delta[i][j]*kron_delta[k][l];
          pdk[i][j][k][l]=pmat[i][j]*kron_delta[k][l];
          kdp[i][j][k][l]=kron_delta[i][j]*pmat[k][l];
          pck[i][j][k][l]=(pmat[i][k]*kron_delta[j][l]+pmat[i][l]*kron_delta[j][k]+pmat[j][l]*kron_delta[i][k]+pmat[j][k]*kron_delta[i][l])/2;
          pdp[i][j][k][l]=pmat[i][j]*pmat[k][l];
        }
      }
    }
  }
  for ( int i=0; i<3; i++ ) {
    for ( int j=0; j<3; j++ ) {
       for ( int k=0; k<3; k++ ) {
          for ( int l=0; l<3; l++ ) {
            hypo_L[i][j][k][l]=a1*kck[i][j][k][l]/2+a2*kdk[i][j][k][l]+a3*(pdk[i][j][k][l]+kdp[i][j][k][l])+a4*pck[i][j][k][l]+a5*pdp[i][j][k][l];
          }
       }
    }
  }
  
  //BEGIN get lambda_act for constant suction;
  double eM=qstatev[5];
  double em=qstatev[4];
  double lambda_act=lambda_use;
  double fm=0, fm_lambda_act=0;
  double re=qstatev[8];                                 // takes re (relative void ratio)
  double a_scan=qstatev[7];
  double sewM=aer*sairentry0*eM0/eM*(aTparam+bTparam*Temper)/(aTparam+bTparam*Trparam);
  double ddum=0;
  double rlambda=give_rlambda(suction, a_scan, 1, sewM, SrM, ddum, ddum);

  give_fm(fm, fm_lambda_act, suction, sewM, re, swelling);
  double pefsat=(signet[0]+signet[4]+signet[8])/3+suction;
  double demdpm=give_demdpm(pefsat, em);

  double p_mean_positive=p_mean;
  if(p_mean_positive<1.0) p_mean_positive=1.0;
  lambda_act=(lambda_use*eM*(1+em)+demdpm*p_mean_positive*(1+eM)*(nparam-lparam*log(p_mean_positive))+fm_lambda_act*demdpm*p_mean_positive*eM)/(eM*(1+em)-rlambda*(1+evoid)*(nparam-lparam*log(p_mean_positive)));
  //END get lambda_act for constant suction;
  
  fs=9*p_mean/2*(1/kappa+1/lambda_act)/Am;
  garray_multiply(&hypo_L[0][0][0][0], &hypo_L[0][0][0][0], fs, 81);

  double hypo_Dsom[3][3];
  garray_set(&hypo_Dsom[0][0], 0, 3*3);
  double kpow=1.7+3.9*sin(phic)*sin(phic);
  double sinphickpow=gscalar_power(sin(phic),kpow);
  double sinphimkpow=gscalar_power(sqrt(sin2phim),kpow);

  double TT[3][3];
  garray_set(&TT[0][0], 0, 3*3);
  gmatrix_ab(&T_str_star[0][0], &T_str_star[0][0], &TT[0][0], 3, 3, 3);
  double TTT[3][3];
  garray_set(&TTT[0][0], 0, 3*3);
  gmatrix_ab( &TT[0][0], &T_str_star[0][0], &TTT[0][0],3, 3, 3);
  double cos3theta=-sqrt(6.)*gtrace(TTT)/gscalar_power(gtrace(TT),1.5);
  if(isotrop)
    cos3theta=-1;

  double Amult=2./3-sqrt(sqrt(sin2phim))*(cos3theta+1)/4.;

  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      hypo_Dsom[i][j]=-T_str_star[i][j]+kron_delta[i][j]*(sinphimkpow-sinphickpow)/(1-sinphickpow)*Amult;
    }
  }
  double Dsomnorm = garray_size(&hypo_Dsom[0][0], 3*3);
  garray_multiply(&hypo_Dsom[0][0], &hypo_Dsom[0][0], (1/Dsomnorm), 3*3);

  double A[3][3][3][3];
  garray_set(&A[0][0][0][0], 0, 81);

  for ( int i=0; i<3; i++ ) {
    for ( int j=0; j<3; j++ ) {
      for ( int k=0; k<3; k++ ) {
        for ( int l=0; l<3; l++ ) {
          A[i][j][k][l]=hypo_L[i][j][k][l]+T[i][j]*kron_delta[k][l]/lambda_act;
        }
      }
    }
  }

  gmatrix_a4b(A, &hypo_Dsom[0][0], &hypo_N[0][0]);
  garray_multiply(&hypo_N[0][0], &hypo_N[0][0],-fddivfdA, 3*3);

  double pdivpsbs=gscalar_power(fddivfdA, (1/alpha_power));
  fu=gscalar_power(pdivpsbs, mparam);
  if(pdivpsbs>1 && pdivpsbs<5)
    fu=gscalar_power(pdivpsbs, 2);
  else if (pdivpsbs>5)
    fu=gscalar_power(5, 2);
  double fu_Hunsat = nparam<0 ? 0 : fu;

  bool m0_OCRindep_collapse=true;
  double cic=1;
  if(m0_OCRindep_collapse && fd<fdsbs) {
    cic=(lambda_act+kappa)*(gscalar_power(2.0,alpha_power)-fd)+2*kappa*fd;
    cic=cic/((lambda_act+kappa)*(gscalar_power(2.0,alpha_power)-fdsbs)+2*kappa*fdsbs);
  }

  garray_set(&H_unsat[0][0], 0, 3*3);
  garray_multiply(&T[0][0], &H_unsat[0][0], dpeds_divpe*cic*lambda_use/lambda_act, 9);
  garray_set(&HT_unsat[0][0], 0, 3*3);
  garray_multiply(&T[0][0], &HT_unsat[0][0], dpedT_divpe*cic*lambda_use/lambda_act, 9);
  garray_multiply(&H_unsat[0][0], &H_unsat[0][0], fu_Hunsat, 9);  //H times fd/fdSOM
  garray_multiply(&HT_unsat[0][0], &HT_unsat[0][0], fu_Hunsat, 9);  //HT times fd/fdSOM
}

double Hypoplasti_unsat_expansive_thermal::give_SrM_new(double suction, double SrM, double eM, double Temper, double sewM, double se, double dsuction, double deM, double dTemper, double a_scan, int flag) {

  double SrM_new=1;
  if(flag==0) {
    if(suction<sewM && suction!=0 && (se/suction)>0) SrM_new =gscalar_power((se/suction), lambdap0);
  }
  else {
    double rlambda_for_ascan=0;
    double rlambda_for_se=0;
    double rlambda=give_rlambda(suction, a_scan, dsuction, sewM, SrM, rlambda_for_ascan, rlambda_for_se);
    double dSrM_ds=0;
    if(suction<-1.0) dSrM_ds=-rlambda*lambdap0/suction*SrM;
    else dSrM_ds=-rlambda*lambdap0/-1*SrM;
    double dSrM_dse=rlambda_for_se*lambdap0*SrM/se;
    double dse_deM=-sairentry0*eM0/(eM*eM)*(aTparam+bTparam*Temper)/(aTparam+bTparam*Trparam);
    double dse_dT=sairentry0*eM0*bTparam/(eM*(aTparam+bTparam*Trparam));

    double dSrM=dSrM_ds*dsuction+dSrM_dse*(dse_deM*deM+dse_dT*dTemper);

    SrM_new=SrM+dSrM;
    if(SrM_new > 1) SrM_new=1.0;
  }
  return SrM_new;

}

double Hypoplasti_unsat_expansive_thermal_2017::give_rlambda(double suction, double a_scan, double dsuction, double sewM, double SrM, double& rlambda_for_ascan, double& rlambda_for_se) {

  double rlambda=0;
  if(a_scan<1-1.0e-10 && a_scan>1.0e-10 && aer<1)
    rlambda=rlambdascan;
  else if(((a_scan>=1-1.0e-10 && dsuction>0) || (a_scan<=1.0e-10 && dsuction<0 && suction<=sewM)) && aer<1)
    rlambda=rlambdascan;
  else if(suction<=sewM)
    rlambda=1;
  rlambda_for_ascan=rlambda;
  rlambda_for_se=(rlambda>1.e-10&&rlambda<(1-1.e-10))?1:rlambda;

  return rlambda;
}

double Hypoplasti_unsat_expansive_thermal_2022::give_rlambda(double suction, double a_scan, double dsuction, double sewM, double SrM, double& rlambda_for_ascan, double& rlambda_for_se) {

  double scanpower=3;
  double SrMlimit=0.75;
  double SrMmax=1.0;
  double WRCpower=1.1;
  double Sepower=3.0;

  double fact=a_scan; //wetting within main curves
  if(dsuction>0) fact=1-a_scan; //drying within main curves

  double rlambda=1;
  rlambda_for_se=1;
  if(fact<1.e-10) fact=0;

  if(aer<1) rlambda=gscalar_power(fact, scanpower);

  if(SrM>1) SrM=1.0;
  if(suction>=sewM) {//drying from 0
      rlambda=0;
      rlambda_for_se=0;
  }
  rlambda_for_ascan=rlambda;

  //wetting to 0, rlambda_for_ascan remains 1 here
  if(dsuction>0 && SrM>SrMlimit) {
      rlambda=gscalar_power((SrMmax-SrM)/(SrMmax-SrMlimit), WRCpower);
      rlambda_for_se=gscalar_power((SrMmax-SrM)/(SrMmax-SrMlimit), Sepower);
  }

  return rlambda;
}

double Hypoplasti_unsat_expansive_thermal::give_rlambda(double suction, double a_scan, double dsuction, double sewM, double SrM, double& rlambda_for_ascan, double& rlambda_for_se) {

  double scanpower=3;
  double SrMlimit=0.2; //was SrMlimit=0.75;
  double SrMmax=1.0;
  double WRCpower=1.0; //was WRCpower=1.1;
  double Sepower=3.0;

  double fact=a_scan; //wetting within main curves
  if(dsuction>0) fact=1-a_scan; //drying within main curves

  double rlambda=1;
  rlambda_for_se=1;
  if(fact<1.e-10) fact=0;

  if(aer<1) rlambda=gscalar_power(fact, scanpower);

  if(SrM>1) SrM=1.0;
  if(suction>=sewM) {//drying from 0
      rlambda=0;
      rlambda_for_se=0;
  }
  rlambda_for_ascan=rlambda;

  //wetting to 0, rlambda_for_ascan remains 1 here
  if(dsuction>0 && SrM>SrMlimit) {
      rlambda=gscalar_power((SrMmax-SrM)/(SrMmax-SrMlimit), WRCpower);
      rlambda_for_se=gscalar_power((SrMmax-SrM)/(SrMmax-SrMlimit), Sepower);
  }

  return rlambda;
}

void Hypoplasti_unsat_expansive_thermal::direct_stiffness_matrix(double strain_gen[], double stress_gen[], double qstatev[],
    double dstrain_gen[], double dtime, double *DDtan_gen, double parms[], int kinc, int flag) {

    double signet[9];
    double Sr=0;
    double othervar[5]={0,0,0,0,0};
    
    convert_from_general(stress_gen, signet, othervar, ngstress-6, 1);
    Sr=othervar[0];

    double ddum9[9];
    convert_from_general(strain_gen, ddum9, othervar, ngstrain-6, 2);
    double suction=othervar[0];
    double Temper=othervar[1];

    double pertinc_all[c_ngstrain];
    garray_set(pertinc_all, 0, ngstrain);
    for(int i=0; i<6; i++) pertinc_all[i]=1.e-6;
    pertinc_all[6]=-sairentry0/100; //ds
    pertinc_all[7]=0.1; //dT
    
    if(flag!=3) { //segmentation fault without this
      for(int i=0; i<ngstrain; i++) {
        if(dstrain_gen[i]<0) pertinc_all[i]=-gscalar_dabs(pertinc_all[i]);
        else pertinc_all[i]=gscalar_dabs(pertinc_all[i]);
      }
    }

    garray_set(DDtan_gen, 0, c_ngstrain*c_ngstress);

    //stiffness matrix by perturbation from 0. Non-mechanical components always by perturbation
    for(int i=0; i<c_ngstrain; i++) {

      double strain_gen_pert[c_ngstrain];
      garray_set(strain_gen_pert, 0, c_ngstrain);
      strain_gen_pert[i]=pertinc_all[i];

      double deps_pert[9];
      double dsignet_pert[9];
      garray_set(dsignet_pert,0,9);
      garray_set(deps_pert,0,9);
      convert_from_general(strain_gen_pert, deps_pert, othervar, (c_ngstrain-6), 2);
      double dsuction_pert=othervar[0];
      double dTemper_pert=othervar[1];
      double Sr_pert=0;
      double dqstatev_pert[c_nstatev];

      calc_fsigq(signet, suction, Temper, qstatev, deps_pert, dsuction_pert, dTemper_pert, 1.0, dsignet_pert, Sr_pert, dqstatev_pert, kinc); //nonmechanical components are isotropic

      double dsig_gen_pert[c_ngstress];
      othervar[0]=Sr_pert-Sr;
      convert_to_general(dsig_gen_pert, dsignet_pert, othervar, (c_ngstress-6), 1);

      for(int j=0; j<c_ngstress; j++) {
        DDtan_gen[j*c_ngstrain+i]=dsig_gen_pert[j]/pertinc_all[i];
      }
    }

    bool perturbatestiffness=false;

    if(!perturbatestiffness) {
      double hypo_L[3][3][3][3];
      double DDtan[3][3][3][3];
      double DDtan_abq[6][6];
      double hypo_N[3][3];
      double H_unsat[3][3];
      double HT_unsat[3][3];
      garray_set(&hypo_L[0][0][0][0],0,81);
      garray_set(&DDtan[0][0][0][0],0,81);
      garray_set(&DDtan_abq[0][0],0,36);
      garray_set(&hypo_N[0][0],0,9);
      garray_set(&H_unsat[0][0],0,9);
      garray_set(&HT_unsat[0][0],0,9);
      double fu=0;
      give_LNH(signet, suction, Temper, qstatev, hypo_L, hypo_N, H_unsat, HT_unsat, false, fu);
      garray_move(&hypo_L[0][0][0][0], &DDtan[0][0][0][0], 81);
      convert4th_to_abaqus(DDtan, DDtan_abq);

      for(int i=0; i<6; i++) {
        for(int j=0; j<6; j++) {
          DDtan_gen[i*c_ngstrain+j]=DDtan_abq[i][j];
        }
      }
    }
}

void Hypoplasti_unsat_expansive_thermal::initialise_variables(double strain_gen[], double stress_gen[], double qstatev[],
    double dstrain_gen[], double dtime, double parms[], int kinc) {

    double signet[9];
    garray_set(signet, 0.0, 9);
    double othervar[5]={0,0,0,0,0};
    convert_from_general(stress_gen, signet, othervar, ngstress-6, 1);
    double Sr=othervar[0];

    double dqstatev[c_nstatev];
    garray_set(dqstatev, 0.0, nstatev);

    double zeroinctens[9]={0,0,0,0,0,0,0,0,0};
    double ddum9[9]={0,0,0,0,0,0,0,0,0};
    double zeroincscal=0;

    convert_from_general(strain_gen, ddum9, othervar, (c_ngstrain-6), 2);

    double suction=othervar[0];
    double Temper=othervar[1];
    qstatev[1]=suction;
    qstatev[3]=Temper;

    double fu=0;

    update_hisv(signet, suction, Temper, qstatev, zeroinctens, zeroincscal, zeroincscal, zeroinctens, Sr, dqstatev, 0, fu, kinc);
    garray_add(qstatev, dqstatev, qstatev, c_nstatev);

    if(qstatev[5]<0) {
      cout<<"Initial value of eM is lower than 0 and this is not allowed, please re-initialise, eM="<<qstatev[5]<<endl;
      exit(1);
    }

    double ddum3333[3][3][3][3];
    double ddum33[3][3];
    double ddum=0;
    garray_set(&ddum3333[0][0][0][0],0,81);
    garray_set(&ddum33[0][0],0,9);
    give_LNH(signet, suction, Temper, qstatev, ddum3333, ddum33, ddum33, ddum33, false, fu);
    if(fu>1) {
        cout<<"Warning, void ratio and stress state were initialised outside SBS, fu="<<fu<<endl;
    }
    stress_gen[6]=Sr;
}

void Hypoplasti_unsat_expansive_thermal::initialise_parameters(double parms[]) {
    //start define parameters
    phic=parms[0]*3.141592653589793238462/180;
    lambda=parms[1];
    kappa=parms[2];
    N=parms[3];
    nuparam=parms[4];
    Ocparam=2;

    nparam=parms[5];
    lparam=parms[6];
    nTparam=parms[7];
    lTparam=parms[8];
    mmparam=parms[9]; //double structure coupling m
    mparam=10;        //collapse m

    alpha_s=parms[10];
    kappam=parms[11];
    smstar=parms[12];
    emstar=parms[13];
    csh=parms[14];

    sairentry0=parms[15];
    eM0=parms[16];
    Trparam=parms[17];
    aTparam=parms[18];
    bTparam=parms[19];

    aer=parms[20];
    lambdap0=parms[21];
    pt=parms[22];

    alphaG=1;
    alphanu=alphaG;
    alphaE=gscalar_power(alphaG,1./0.8);

    rlambdascan=0.1;
    //end define parameters
}

long General_model::calc_rkf(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[], double dtime,
    double parms[], double rkf_statev[], double rkf_parms[], int kinc) {

    long errorrkf = false;
    //if not specified, initialise rkf_parms & rkf_statev to default values
    if(rkf_parms==NULL) {

      //Order of parameters in the array rkf_parms
      //rkf_parms[0] = err_sig, default 1.e-3
      //rkf_parms[1] = h_min, default  1.e-17
      //rkf_parms[2] = ni, default 1000
      //rkf_parms[3] = sv_rkf_zero, default  1.e-3
      //rkf_parms[4] = rkt, default  4
      double rkf_parms_default[nrkf_parms] = {1.0e-3,1.0e-17,1000,1.0e-3,4}; //{1.0e-5,1.0e-17,1,1.0e-3,7};
      rkf_parms=&rkf_parms_default[0];
    }
    if(rkf_statev==NULL) {
      double rkf_statev_default[nrkf_statev] = {0,0,1,0.1};
      rkf_statev=&rkf_statev_default[0];
    }

    int rkt = rkf_parms[4];
    switch (rkt){
      case 1: { //case rkt=1 forward Euler
        double signet[9];
        garray_set(signet, 0.0, 9);
        double othervar[5]={0,0,0,0,0};
        convert_from_general(stress_gen, signet, othervar, ngstress-6, 1);
        double Sr=othervar[0];

        double ddum9[9];
        convert_from_general(strain_gen, ddum9, othervar, ngstrain-6, 2);
        double suction=othervar[0];
        double Temper=othervar[1];

        double dstress_gen[max_ngstress];
        double dqstatev[max_nstatev];
        garray_set(dstress_gen, 0.0, ngstress);
        garray_set(dqstatev, 0.0, nstatev);

        double gmod_deps[9];
        garray_set(gmod_deps, 0, 9);
        convert_from_general(dstrain_gen, gmod_deps, othervar, ngstrain-6, 2);
        double dsuction=othervar[0];
        double dTemper=othervar[1];

        double dsignet[9];
        garray_set(dsignet, 0, 9);
        double Sr_new=0;

        errorrkf = calc_fsigq(signet, suction, Temper, qstatev, gmod_deps, dsuction, dTemper, dtime, dsignet, Sr_new, dqstatev, kinc);

        othervar[0]=Sr_new-Sr;
        garray_set(dstress_gen,0,ngstress);
        convert_to_general(dstress_gen, dsignet, othervar,(ngstress-6),1);

        //Update external variables
        for (int i=0; i<ngstress; i++) stress_gen[i] += dstress_gen[i];
        for (int i=0; i<nstatev; i++) qstatev[i] += dqstatev[i];
        rkf_statev[0] += 1;
        rkf_statev[1] += 1;
        rkf_statev[2] = 1.0;
        break;
      }
      case 2: //forward euler with adaptive substepping
        errorrkf = adfwdeuler(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
        break;
      case 3:
        errorrkf = rkf23(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
        break;
      case 4:
        errorrkf = rkf23bs(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
        break;
      case 5:
        errorrkf = rkf34(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
        break;
      case 6:
        errorrkf = rkf45(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
        break;
      case 7: //forward euler with fixed number of substeps
        errorrkf = fwdeulerfsub(strain_gen, stress_gen, qstatev, dstrain_gen, dtime, parms, rkf_statev, rkf_parms, kinc);
        break;
      default:
        cout << "soil model called with unknown Runge-Kutta scheme" << endl;
    }
    //update also strain, stress and state is updated in RKF call
    for (int i=0; i<ngstrain; i++) strain_gen[i] += dstrain_gen[i];
    return(errorrkf);
}

//Functions for abaqus-full convention conversion
void General_model::convert_from_general( double general[], double stressstrain[], double othervar[], int nothv, int flag ) {
//flag==1: stress-like; flag==2: strain-like
//4 -> 12 , 5 -> 13, 6 -> 23
  stressstrain[0]=general[0];
  stressstrain[4]=general[1];
  stressstrain[8]=general[2];
  if(flag == 1) {
    stressstrain[1]=general[3];
    stressstrain[2]=general[4];
    stressstrain[5]=general[5];
    stressstrain[3]=general[3];
    stressstrain[6]=general[4];
    stressstrain[7]=general[5];
  }
  if(flag == 2) {
    stressstrain[1]=general[3]/2;
    stressstrain[2]=general[4]/2;
    stressstrain[5]=general[5]/2;
    stressstrain[3]=general[3]/2;
    stressstrain[6]=general[4]/2;
    stressstrain[7]=general[5]/2;
  }
  for(int i=0; i<nothv; i++) othervar[i]=general[6+i];
}
void General_model::convert_to_general( double general[], double stressstrain[], double othervar[], int nothv, int flag ){
//flag==1: stress-like; flag==2: strain-like
//4 -> 12 , 5 -> 13, 6 -> 23
  general[0]=stressstrain[0];
  general[1]=stressstrain[4];
  general[2]=stressstrain[8];
  if(flag == 1) {
    general[3]=stressstrain[1];
    general[4]=stressstrain[2];
    general[5]=stressstrain[5];
  }
  if(flag == 2) {
    general[3]=2*stressstrain[1];
    general[4]=2*stressstrain[2];
    general[5]=2*stressstrain[5];
  }
  for(int i=0; i<nothv; i++) general[6+i]=othervar[i];
}
void General_model::convert4th_to_abaqus( double DDfull[3][3][3][3], double DDabq[6][6]){
  int fi=0;
  int fj=0;
  int fk=0;
  int fl=0;
  for ( int ai=0; ai<6; ai++ ) {
    for ( int aj=0; aj<6; aj++ ) {
        map_indices_abqtofull( fi, fj, ai);
        map_indices_abqtofull( fk, fl, aj);
        DDabq[ai][aj]=(DDfull[fi][fj][fk][fl]+DDfull[fi][fj][fl][fk])/2;
    }
  }
}
void General_model::map_indices_fulltoabq( int fulli, int fullj, int &abq){
   if(fulli == 0 && fullj == 0) abq=0;
   if(fulli == 1 && fullj == 1) abq=1;
   if(fulli == 2 && fullj == 2) abq=2;
   if(fulli == 0 && fullj == 1) abq=3;
   if(fulli == 0 && fullj == 2) abq=4;
   if(fulli == 1 && fullj == 2) abq=5;
   if(fulli == 1 && fullj == 0) abq=3;
   if(fulli == 2 && fullj == 0) abq=4;
   if(fulli == 2 && fullj == 1) abq=5;
}
void General_model::map_indices_abqtofull( int &fulli, int &fullj, int abq){
   if(abq==0) {
     fulli=0;
     fullj=0;
   };
   if(abq==1) {
     fulli=1;
     fullj=1;
   };
   if(abq==2) {
     fulli=2;
     fullj=2;
   };
   if(abq==3) {
     fulli=0;
     fullj=1;
   };
   if(abq==4) {
     fulli=0;
     fullj=2;
   };
   if(abq==5) {
     fulli=1;
     fullj=2;
   };
}

//Functions for mathematical operations
double General_model::gscalar_power( double a, double b ) {
  /*  if(b==0.5 && a<0) {
	cout<<"Attempting to calculate square root of a negative number";
	exit(1);
  }
  double result=0.;
  if ( b==0. )
    result = 1.;
  else
  result = pow( a, b );
  return result;*/
  long errno_input=errno;
  errno=0;
  double ret = pow(a,b);
  if (test_math_err() == 0) {
    errno=errno_input;
    return ret;
  }
  else if (a>=0 && !isnan(a) && !isnan(b)) return 0;
  else
  {
    if(debug) cout<<"Error in calculating power "<<a<<"^"<<b<<endl;
    return(0.0);
  }
}
double General_model::gscalar_dabs(double k) {
  double tmp;
  if(k>= 0)
    tmp=k;
  if(k<0)
    tmp=-1*k;
  return tmp;
}
void General_model::garray_add( double a[], double b[], double c[], long int n ) {
  long int i=0;
  for ( i=0; i<n; i++ ) c[i] = a[i] + b[i];
}
double General_model::garray_inproduct( double a[], double b[], long int n ) {
  long int i=0;
  double result=0.;
  for ( i=0; i<n; i++ ) result += a[i]*b[i];
  return result;
}
void General_model::garray_move( double from[], double to[], long int n ) {
  //  long int i=0;
  //  for ( i=0; i<n; i++ ) to[i] = from[i];
  memcpy(to, from, sizeof(*from)*n);
}
void General_model::garray_multiply( double a[], double b[], double c, long int n ) {
  long int i=0;
  for ( i=0; i<n; i++ ) b[i] = c * a[i];
}
void General_model::garray_set( double *ptr, double value, long int n ) {				//@GS initialises the array to "value" (e.g. to zero)
  // long int i=0;
  // for ( i=0; i<n; i++ ) *(ptr+i) = value;
  if (value == 0.0)
    memset(ptr, 0, sizeof(*ptr)*n);
  else
    for (long int i=0; i<n; i++ ) *(ptr+i) = value;
}
double General_model::garray_size( double a[], long int n ) {
  double size=0.;
  size = sqrt( gscalar_dabs( garray_inproduct(a,a,n) ) );
  return size;
}
double General_model::gtrace( double a[3][3]) {
	return(a[0][0]+a[1][1]+a[2][2]);
}
double General_model::gmatrix_determinant( double a[], long int n ) {

  double result=0.;

  if ( n==1 )
    result = a[0];
  else if ( n==2 )
    result = a[0]*a[3] - a[1]*a[2];
  else {
    if(n == 3) {
      result = a[0]*(a[4]*a[8]-a[7]*a[5]) -
      a[1]*(a[3]*a[8]-a[6]*a[5]) + a[2]*(a[3]*a[7]-a[6]*a[4]);
    }
    else {
      cout<<"n must be 1-3"<<endl;
      exit(1);
    }
  }
  return result;
}
void General_model::gmatrix_ab( double *a, double *b, double *c, long int n, long int m, long int k ) { // c[n][k] = a[n][m] * b[m][k]
    long int i=0, j=0, l=0;

  for ( i=0; i<n; i++ ) {
    for ( j=0; j<k; j++ ) {
      *( c + i*k + j ) = 0;
      for ( l=0; l<m; l++ ) {
        *( c + i*k + j ) += ( *( a + i*m + l ) ) *
                             ( *( b + l*k + j ) );
       }
    }
  }
}
void General_model::gmatrix_a4b( double a[3][3][3][3], double b[], double c[] ) {
  long int i=0, j=0, k=0, l=0;

  for ( i=0; i<3; i++ ) {
    for ( j=0; j<3; j++ ) {
      c[i*3+j] = 0.;
      for ( k=0; k<3; k++ ) {
        for ( l=0; l<3; l++ ) {
          c[i*3+j] += a[i][j][k][l] * b[k*3+l];
        }
      }
    }
  }

}
double General_model::round_to_digits(double value, int digits) {
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;
}
long General_model::test_math_err() {
  if (errno == EDOM) return 1;
  else if(errno == ERANGE) return 2;
  else return 0;
}

//@GS Functions adapted from SIFEL
/**
  The function integrates ODEs by adaptive forward Euler rule.

  @param strain_gen[ngstrain] - (input): controlled variables [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]^T
  @param stress_gen[c_ngstress] - (input): implied variables, will not be updated [signet11, signet22, signet33, signet12, signet13, signet23, Sr]^T. "signet" is net stress, that is total stress minus krondelta x ua.
  @param qstatev[c_nstatev] - (input): state variable vector, will be updated only for flag==0.
  @param dstrain_gen[ngstrain] - (input): increment of controlled variables [deps11, deps22, deps33, dgamma12, dgamma13, dgamma23, ds, dT]^T
  @param dtime - (input): increment of time
  @param parms[c_nparms] - array of model parameters (input)
  @param rkf_statev[nrkf_statev] - (input): array of state variables from SIFEL for Runge-Kutta
  @param rkf_parms[nrkf_parms] - (input): array of parameters from SIFEL for Runge-Kutta

  @retval 0 - on success
  @retval 1 - in the case of substepping failure => global step length must be decreased

  Created by Tomas Koudelka, 6.10.2015; modified by G. Scaringi, January 2019 //@GS
*/
long General_model::adfwdeuler(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[],
      double dtime, double parms[], double rkf_statev[], double rkf_parms[], int kinc) {
  //{ initialization of variables
  long errorfsigq = 0;
  long neval = rkf_statev[0];
  long nstep = rkf_statev[1];
  double mindt = rkf_statev[2];
  double dtsub = rkf_statev[3];
  double h, h_opt, h_opt_sig, h_opt_stv;
  double dt = 0.0;
  double err_sig = rkf_parms[0];
  double err_stv = rkf_parms[0];
  double err_stv_req = rkf_parms[0];
  double err_stv_max = 1.0e-3;
  double norm_sig, norm_stv, norm_stvo, aux, norm_r;
  double dtmin = rkf_parms[1];
  double idtsub = dtsub;
  double sv_rkf_zero = rkf_parms[3]; //needed
  long k, j;
  long ni = rkf_parms[2];
  long stopk = 0;
  long max_ngcomp = max_ngstress + max_nstatev;
  long ngcomp = ngstress + nstatev;
  double k1[max_ngstress + max_nstatev], k2[max_ngstress + max_nstatev], y1[max_ngstress + max_nstatev], y2[max_ngstress + max_nstatev], y_hat[max_ngstress + max_nstatev], y_til[max_ngstress + max_nstatev], y_err[max_ngstress + max_nstatev];
  double strain_gen_aux[max_ngstrain], dstrain_gen_aux[max_ngstrain];
  long herr;
  // additional variables (for calc_fsigq)
  double g_suction=0, g_dsuction=0, g_Temper=0, g_dTemper=0, g_Sr=0, g_Sr_new=0;
  double g_signet[9], g_deps[9], g_ddum9[9], g_othervar[5], g_dsignet[9], g_dSr[1];
  long error=0;
    //initialising variables
  for (j=0; j<max_ngcomp; j++) {
     k1[j]=0;
     k2[j]=0;
     y1[j]=0;
     y2[j]=0;
     y_hat[j]=0;
     y_til[j]=0;
     y_err[j]=0;
    if(j<max_ngstrain) {
      strain_gen_aux[j]=0;
      dstrain_gen_aux[j]=0;
    }
  }
  //}
  //{ preparation
  // set initial substep length
  if (!dtsub) h = 1.0;
  else h = dtsub;
  norm_stvo = 1.e30;
  // y1 has to be initialized to the values of stress vector sig and state variables
  garray_rcopy(stress_gen, 0, y1, 0, ngstress);
  garray_rcopy(qstatev, 0, y1, ngstress, nstatev);
  garray_copy(strain_gen, strain_gen_aux, ngstrain);
  //}
  //{ ITERATION LOOP
  for (k=0; k<ni && dt<1.0; k++) {
    //{ preparation
    // strain_gen_aux = strain_gen + dt*dstrain_gen
    garray_addmult(strain_gen, dstrain_gen, dt, strain_gen_aux, ngstrain);
    // dstrain_gen_aux = h*dstrain_gen;
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);
    if (h==0.0)
      return 2;
    //}
    //{ first call to calc_fsigq
    // calculate Forward Euler
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k1+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k1, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    garray_add(y1, k1, y_til, ngcomp);
    garray_addmult(NULL, dstrain_gen, h/2.0, dstrain_gen_aux, ngstrain);
    //}
    //{ second call to calc_fsigq
    // calculate Forward Euler for half step
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k2+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k2, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    garray_add(y1, k2, y2, ngcomp);
    garray_addmult(strain_gen, dstrain_gen, dt+0.5*h, strain_gen_aux, ngstrain);
    //}
    //{ third call to calc_fsigq
    convert_from_general (y2, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y2+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k2+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k2, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ finalization and error computation
    garray_add(y2, k2, y_hat, ngcomp);
    // increase the number of model evaluations
    neval += 3;
    rkf_statev[0] = neval;
    // solution error is difference between half step solution y_hat and solution with full step y_til
    garray_subtract(y_hat, y_til, y_err, ngcomp);
    // stress tensor norm stored in the Voigt notation, additional components in y_hat are ignored
    norm_sig = tensor_dbldot_prod(y_hat, y_hat, 2.0);
    norm_sig = sqrt(norm_sig);
    // calculate normalized errors of stress components
    aux = 1.0/norm_sig;
    for (j=0; j<6; j++) y_err[j] = fabs(y_err[j])*aux;
    // calculate normalized error of Sr and particular state variables
    for(j=6; j<ngcomp; j++) {
      if (fabs(y_hat[j]) < sv_rkf_zero) y_err[j] = fabs(y_err[j]);
      else y_err[j] = fabs(y_err[j])/y_hat[j];
    }
    // calculate norm of residual
    norm_r = garray_snorm(y_err, 0, ngcomp);
    // calculate norm of stress components
    norm_sig = garray_snorm(y_err, 0, 6);
    // zero residual of suction
    y_err[ngstress+1] = 0.0;
    // zero residual of temperature
    y_err[ngstress+3] = 0.0;
    // zero residual of swelling indicator
    y_err[ngstress+nstatev-1] = 0.0;
    // calculate norm of remaining components of generalized stress vector and state variables
    norm_stv = garray_snorm(y_err, 6, 1);
    // check divergence of state variable residual (norm_stv >=norm_stvo)
    if (norm_stv >= norm_stvo) err_stv = err_stv_max; // in the case of the divergence, use maximum tolerance allowed for the state variables rather then required one (err_stv_req)
    // calculate optimum substep length
    if (norm_sig != 0.0) h_opt_sig = 0.9*h*pow(err_sig/norm_sig, 1.0/3.0);
    else  h_opt_sig = 1.0;
    if (norm_stv != 0.0) h_opt_stv = 0.9*h*pow(err_stv/norm_stv, 1.0/3.0);
    else h_opt_stv = 1.0;
    h_opt = min2(h_opt_sig, h_opt_stv);
    //}
    //{ rejection or acceptance of substep
    if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
      // substep is accepted, !FE solution taken
      garray_swap(y_hat, y1, ngcomp);
      // set the cumulative substep length
      dt += h;
      if (stopk) break; // dt = 1 => whole interval has been integrated
      // select optimal RK substep length
      h = min2(4.0*h, h_opt);
      if (h > 1.0) h = 1.0;
      // store the last attained RK substep length to be the initial substep length in the next RKF method call
      dtsub = h;
      rkf_statev[3] = dtsub;
      // check the maximum RK substep length
      if (1.0-dt < h) {
        h = 1.0-dt;
        stopk = 1;
      }
      // set initial values for handling of state variable residual in the next RKF step
      norm_stvo = 1.e30;
      err_stv = err_stv_req;
    } else {
      // substep is not accepted => reduce substep length
      h = max2(1.0/4.0*h, h_opt);
      // check for minimum step size
      if (h < dtmin) {
        dtsub = idtsub;
        rkf_statev[3] = dtsub;
        return 1;
      }
      if (mindt > h) mindt = h;
      rkf_statev[2] = mindt;
      // actualize attained norm of state variable residual
      norm_stvo = norm_stv;
    }
    //}
  }
  //} end of the iteration loop
  //{ finalization (storing results, adapting time step)
  if (dt > 1.0) dt = 1.0; // cut-off rounding error in the cumulative RK substep length
  if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
    // store the attained stresses - they are stored in y1 because the content of y_hat and y1 is swapped
    // in the case of accepted substep
    garray_rcopy(y1, 0, stress_gen, 0, ngstress);
    // state variables must be also integrated by Runge-Kutta => they must be copied from the end of y1 vector
    garray_rcopy(y1, ngstress, qstatev, 0, nstatev);
    nstep += k;
    rkf_statev[1] = nstep;
  } else {
    // RKF method cannot integrate the function with required error was in the given maximum number of substeps
    // dt must be less then 1.0
    dtsub = idtsub;
    rkf_statev[3] = dtsub;
    return 1;
  }
  //}
  return 0;
}

//Forward Euler with fixed number of substeps
long General_model::fwdeulerfsub(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[],
      double dtime, double parms[], double rkf_statev[], double rkf_parms[], int kinc) {
  //{ initialization of variables
  long errorfsigq = 0;
  long neval = rkf_statev[0];
  long nstep = rkf_statev[1];
  double dtsub = rkf_statev[3];
  double h;
  double dt = 0.0;
  double dtmin = rkf_parms[1];
  double idtsub = dtsub;
  double sv_rkf_zero = rkf_parms[3]; //needed
  long j, k;
  long ni = rkf_parms[2];
  long stopk = 0;
  long max_ngcomp = max_ngstress + max_nstatev;
  long ngcomp = ngstress + nstatev;
  double k1[max_ngstress + max_nstatev], y1[max_ngstress + max_nstatev];
  double strain_gen_aux[max_ngstrain], dstrain_gen_aux[max_ngstrain];
  // additional variables (for calc_fsigq)
  double g_suction=0, g_dsuction=0, g_Temper=0, g_dTemper=0, g_Sr=0, g_Sr_new=0;
  double g_signet[9], g_deps[9], g_ddum9[9], g_othervar[5], g_dsignet[9], g_dSr[1];
  long error=0;
    //initialising variables
  for (j=0; j<max_ngcomp; j++) {
     k1[j]=0;
     y1[j]=0;
    if(j<max_ngstrain) {
      strain_gen_aux[j]=0;
      dstrain_gen_aux[j]=0;
    }
  }
  //}
  //{ preparation
  // set initial substep length
  h = 1/(double)ni;
  // y1 has to be initialized to the values of stress vector sig and state variables
  garray_rcopy(stress_gen, 0, y1, 0, ngstress);
  garray_rcopy(qstatev, 0, y1, ngstress, nstatev);
  garray_copy(strain_gen, strain_gen_aux, ngstrain);
  //}
  //{ ITERATION LOOP
  for (k=0; k<ni && dt<1.0; k++) {
    //{ preparation
    // strain_gen_aux = strain_gen + dt*dstrain_gen
    garray_addmult(strain_gen, dstrain_gen, dt, strain_gen_aux, ngstrain);
    // dstrain_gen_aux = h*dstrain_gen;
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k1+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k1, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq) return 1;
    garray_add(y1, k1, y1, ngcomp);

    neval += 3;
    rkf_statev[0] = neval;
    // substep is accepted, !FE solution taken
    dt += h;
    if (stopk) break; // dt = 1 => whole interval has been integrated
    // select optimal RK substep length
    // store the last attained RK substep length to be the initial substep length in the next RKF method call
    dtsub = h;
    rkf_statev[3] = dtsub;
    // check the maximum RK substep length
    if (1.0-dt < h) {
      h = 1.0-dt;
      stopk = 1;
    }
  }
  //} end of the iteration loop
  //{ finalization (storing results, adapting time step)
  if (dt > 1.0) dt = 1.0; // cut-off rounding error in the cumulative RK substep length
  // store the attained stresses - they are stored in y1 because the content of y_hat and y1 is swapped
  // in the case of accepted substep
  garray_rcopy(y1, 0, stress_gen, 0, ngstress);
  // state variables must be also integrated by Runge-Kutta => they must be copied from the end of y1 vector
  garray_rcopy(y1, ngstress, qstatev, 0, nstatev);
  nstep += k;
  rkf_statev[1] = nstep;
  //}
  return 0;
}

/**
  The function integrates ODEs by the kind of Runge-Kutta-Fehlberg 23 method.

  @param strain_gen[ngstrain] - (input): controlled variables [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]^T
  @param stress_gen[max_ngstress] - (input): implied variables, will not be updated [signet11, signet22, signet33, signet12, signet13, signet23, Sr]^T. "signet" is net stress, that is total stress minus krondelta x ua.
  @param qstatev[max_nstatev] - (input): state variable vector, will be updated only for flag==0.
  @param dstrain_gen[ngstrain] - (input): increment of controlled variables [deps11, deps22, deps33, dgamma12, dgamma13, dgamma23, ds, dT]^T
  @param dtime - (input): increment of time
  @param parms[c_nparms] - array of model parameters (input)
  @param rkf_statev[nrkf_statev] - (input): array of state variables from SIFEL for Runge-Kutta
  @param rkf_parms[nrkf_parms] - (input): array of parameters from SIFEL for Runge-Kutta

  @retval 0 - on success
  @retval 1 - in the case of substepping failure => global step length must be decreased

  Created by Tomas Koudelka, 6.10.2015; modified by G. Scaringi, January 2019 //@GS
*/
long General_model::rkf23(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[],
      double dtime, double parms[], double rkf_statev[], double rkf_parms[], int kinc) {
  //{ initialization of variables
  long errorfsigq = 0;
  long neval = rkf_statev[0];
  long nstep = rkf_statev[1];
  double mindt = rkf_statev[2];
  double dtsub = rkf_statev[3];
  double h, h_opt, h_opt_sig, h_opt_stv;
  double dt = 0.0;
  double err_sig = rkf_parms[0];
  double err_stv = rkf_parms[0];
  double norm_sig, norm_stv, /*norm_stvo,*/ aux, norm_r;
  double dtmin = rkf_parms[1];
  double idtsub = dtsub;
  double sv_rkf_zero = rkf_parms[3]; //needed
  long k, j;
  long ni = rkf_parms[2];
  long stopk = 0;
  long max_ngcomp = max_ngstress + max_nstatev;
  long ngcomp = ngstress + nstatev;
  long stvcomp;
  double k1[max_ngstress+ max_nstatev], k2[max_ngstress + max_nstatev], k3[max_ngstress + max_nstatev], y1[max_ngstress + max_nstatev], y2[max_ngstress + max_nstatev], y3[max_ngstress + max_nstatev], y_hat[max_ngstress + max_nstatev], y_til[max_ngstress + max_nstatev], y_err[max_ngstress + max_nstatev];
  double strain_gen_aux[max_ngstrain], dstrain_gen_aux[max_ngstrain];
  long herr;
  // additional variables (for calc_fsigq)
  double g_suction=0, g_dsuction=0, g_Temper=0, g_dTemper=0, g_Sr=0, g_Sr_new=0;
  double g_signet[9], g_deps[9], g_ddum9[9], g_othervar[5], g_dsignet[9], g_dSr[1];
  long error=0;
    //initialising variables
  for (j=0; j<max_ngcomp; j++) {
     k1[j]=0;
     k2[j]=0;
     k3[j]=0;
     y1[j]=0;
     y2[j]=0;
     y3[j]=0;
     y_hat[j]=0;
     y_til[j]=0;
     y_err[j]=0;
    if(j<max_ngstrain) {
      strain_gen_aux[j]=0;
      dstrain_gen_aux[j]=0;
    }
  }
  //}
  //{ preparation
  // set initial substep length
  if (!dtsub) h = 1.0;
  else h = dtsub;
  // y1 has to be initialized to the values of stress vector sig and state variables
  garray_rcopy(stress_gen, 0, y1, 0, ngstress);
  garray_rcopy(qstatev, 0, y1, ngstress, nstatev);
  garray_copy(strain_gen, strain_gen_aux, ngstrain);
  //}
  //{ ITERATION LOOP
  for (k=0; k<ni && dt<1.0; k++) {
    //{ preparation
    // the number state variables stored after stress components in the vector y_hat that should be checked for tolerance
    stvcomp = 1;
    // strain_gen_aux = strain_gen + dt*dstrain_gen
    garray_addmult(strain_gen, dstrain_gen, dt, strain_gen_aux, ngstrain);
    // dstrain_gen_aux = h*dstrain_gen;
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);
    if (h==0.0)
      return 2;
    //}
    //{ first call to calc_fsigq
    // calculate initial vector of coefficients for the first stage of Runge-Kutta method
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k1+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k1, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // correct ascan value after integration by RKF
    // ascan value was cut because it was out of prescribed range <0;1>
    // do not check tolerance for state variables, i.e. the number of checked state variables is zero
    aux = y1[ngstress+7]+k1[ngstress+7];
    if (((fabs(1.0-aux) < 1.0e-7) && (fabs(1.0-y1[ngstress+7]) > 1.0e-7)) || ((fabs(aux) < 1.0e-7) && (y1[ngstress+7] > 1.0e-7))) stvcomp = 0;
    // correct em value after integration by RKF
    aux = y1[ngstress+4]+k1[ngstress+4];
    // em value was cut on 0.01
    // do not check tolerance for state variables, i.e. the number of checked state variables is zero
    if (((fabs(0.01-aux) < 1.0e-7) && (fabs(y1[ngstress+7]-0.01) > 1.0e-7))) stvcomp = 0;

    // y2(i) = y1(i) + 1/2*h*k1(i)
    garray_addmult(y1, k1, 0.5, y2, ngcomp);
    // strain_gen_aux = strain_gen + (dt+0.5*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+0.5*h, strain_gen_aux, ngstrain);
    //}
    //{ second call to calc_fsigq
    // calculate coefficients for the second stage of Runge-Kutta method
    convert_from_general (y2, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y2+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k2+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k2, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y3(i) = y1 - h*k1(i) + 2.0*h*k2(i)
    for(j=0; j<ngcomp; j++) y3[j]=y1[j]+(-k1[j] + 2.0*k2[j]);
    if (stvcomp == 0) {
      y3[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y3[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+h)*dstrain_gen
    garray_addmult(strain_gen, dstrain_gen, dt+h, strain_gen_aux, ngstrain);
    //}
    //{ third call to calc_fsigq
    // calculate coefficients for the third stage of Runge-Kutta method
    convert_from_general (y3, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y3+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k3+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k3, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ finalization and error computation
    // calculate third stage of embedded Runge-Kutta-Fehlberg method
    for (j=0; j<ngcomp; j++) {
      //      y_til(j)=y1(j) + h*k2(j);
      //      y_hat(j)=y1(j) + h*((1.0/6.0)*k1(j) + (2.0/3.0)*k2(j) + (1.0/6.0)*k3(j));
      y_til[j]=y1[j] + k2[j];
      y_hat[j]=y1[j] + ((1.0/6.0)*k1[j] + (2.0/3.0)*k2[j] + (1.0/6.0)*k3[j]);
    }
    if (stvcomp == 0) {
      // use forward Euler rule for selected state variables
      y_hat[ngstress+7] = y_til[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y_hat[ngstress+4] = y_til[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    // increase the number of model evaluations
    neval += 3;
    rkf_statev[0] = neval;
    // solution error is difference between 3rd and 2nd stages
    garray_subtract(y_hat, y_til, y_err, ngcomp);
    // stress tensor norm stored in the Voigt notation, additional components in y_hat are ignored
    norm_sig = tensor_dbldot_prod(y_hat, y_hat, 2.0);
    norm_sig = sqrt(norm_sig);
    // calculate normalized errors of stress components
    aux = 1.0/norm_sig;
    for (j=0; j<6; j++) y_err[j] = fabs(y_err[j])*aux;
    // calculate normalized error of Sr and particular state variables
    for(j=6; j<ngcomp; j++) {
      if (fabs(y_hat[j]) < sv_rkf_zero) y_err[j] = fabs(y_err[j]);
      else y_err[j] = fabs(y_err[j])/y_hat[j];
    }
    // calculate norm of residual
    norm_r = garray_snorm(y_err, 0, ngcomp);
    // calculate norm of stress components
    norm_sig = garray_snorm(y_err, 0, 6);
    // zero residual of suction
    y_err[ngstress+1] = 0.0;
    // zero residual of temperature
    y_err[ngstress+3] = 0.0;
    // zero residual of swelling indicator
    y_err[ngstress+nstatev-1] = 0.0;
    // calculate norm of remaining components of generalized stress vector and state variables eventually
    norm_stv = garray_snorm(y_err, 6, stvcomp);
    // calculate optimum substep length
    if (norm_sig != 0.0) h_opt_sig = 0.9*h*pow(err_sig/norm_sig, 1.0/3.0);
    else h_opt_sig = 1.0;
    if (norm_stv != 0.0) h_opt_stv = 0.9*h*pow(err_stv/norm_stv, 1.0/3.0);
    else h_opt_stv = 1.0;
    h_opt = min2(h_opt_sig, h_opt_stv);
    //}
    //{ rejection or acceptance of substep
    if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
      // substep is accepted,
      // update y1 to the values of attained generalized stress vector and state variables
      garray_swap(y_hat, y1, ngcomp);
      // set the cumulative substep length
      dt += h;
      if (stopk) // dt = 1 => whole interval has been integrated
        break;
      // select optimal RK substep length
      h = min2(4.0*h, h_opt);
      if (h > 1.0)
        h = 1.0;
      // store the last attained RK substep length to be the initial substep length in the next RKF method call
      dtsub = h;
      rkf_statev[3] = dtsub;
      // check the maximum RK substep length
      if (1.0-dt < h) {
        h = 1.0-dt;
        stopk = 1;
      }
    }
    else {
      // substep is not accepted => reduce substep length
      h = max2(1.0/4.0*h, h_opt);
      // check for minimum step size
      if (h < dtmin) {
        dtsub = idtsub;
        rkf_statev[3] = dtsub;
        return 1;
      }
      if (mindt > h) mindt = h;
      rkf_statev[2] = mindt;
      // zero possibly set flag for the main loop breaking due to limit value dt
      stopk = 0;
    }
    //}
  }
  //} end of the iteration loop
  //{ finalization (storing results, adapting time step)
  if (dt > 1.0) dt = 1.0; // cut-off rounding error in the cumulative RK substep length
  if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
    // store the attained stresses - they are stored in y1 because the content of y_hat and y1 is swapped
    // in the case of accepted substep
    garray_rcopy(y1, 0, stress_gen, 0, ngstress);
    // state variables must be also integrated by Runge-Kutta => they must be copied from the end of y1 vector
    garray_rcopy(y1, ngstress, qstatev, 0, nstatev);
    nstep += k;
    rkf_statev[1] = nstep;
  } else {
    // RKF method cannot integrate the functvector &ion with required error was in the given maximum number of substeps
    // dt must be less then 1.0
    dtsub = idtsub;
    rkf_statev[3] = dtsub;
    return 1;
  }
  //}
  return 0;
}

/**
  The function integrates ODEs by the kind of Runge-Kutta-Fehlberg 23 method proposed by
  Bogacki & Shampine.

  @param strain_gen[ngstrain] - (input): controlled variables [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]^T
  @param stress_gen[max_ngstress] - (input): implied variables, will not be updated [signet11, signet22, signet33, signet12, signet13, signet23, Sr]^T. "signet" is net stress, that is total stress minus krondelta x ua.
  @param qstatev[max_nstatev] - (input): state variable vector, will be updated only for flag==0.
  @param dstrain_gen[ngstrain] - (input): increment of controlled variables [deps11, deps22, deps33, dgamma12, dgamma13, dgamma23, ds, dT]^T
  @param dtime - (input): increment of time
  @param parms[c_nparms] - array of model parameters (input)
  @param rkf_statev[nrkf_statev] - (input): array of state variables from SIFEL for Runge-Kutta
  @param rkf_parms[nrkf_parms] - (input): array of parameters from SIFEL for Runge-Kutta

  @retval 0 - on success
  @retval 1 - in the case of substepping failure => global step length must be decreased

  Created by Tomas Koudelka, 6.10.2015; modified by G. Scaringi, January 2019 //@GS
*/
long General_model::rkf23bs(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[],
  double dtime, double parms[], double rkf_statev[], double rkf_parms[], int kinc) {
  //{ initialization of variables
  long errorfsigq = 0;
  long neval = rkf_statev[0];
  long nstep = rkf_statev[1];
  double mindt = rkf_statev[2];
  double dtsub = rkf_statev[3];
  double h, h_opt, h_opt_sig, h_opt_stv, h_old;
  double dt = 0.0;
  double err_sig = rkf_parms[0];
  double err_stv = rkf_parms[0];
  double norm_sig, norm_stv, /*norm_stvo,*/ aux, norm_r;
  double dtmin = rkf_parms[1];
  double idtsub = dtsub;
  double sv_rkf_zero = rkf_parms[3]; //needed
  long k, j;
  long ni = rkf_parms[2];
  long stopk = 0;
  long max_ngcomp = max_ngstress + max_nstatev;
  long ngcomp = ngstress + nstatev;
  long stvcomp;
  double k1[max_ngstress + max_nstatev], k2[max_ngstress + max_nstatev], k3[max_ngstress + max_nstatev], k4[max_ngstress + max_nstatev];
  double y1[max_ngstress + max_nstatev], y2[max_ngstress + max_nstatev], y3[max_ngstress + max_nstatev], y_hat[max_ngstress + max_nstatev], y_til[max_ngstress + max_nstatev], y_err[max_ngstress + max_nstatev];
  double strain_gen_aux[max_ngstrain], dstrain_gen_aux[max_ngstrain];
  long herr;
  // additional variables (for calc_fsigq)
  double g_suction=0, g_dsuction=0, g_Temper=0, g_dTemper=0, g_Sr=0, g_Sr_new=0;
  double g_signet[9], g_deps[9], g_ddum9[9], g_othervar[5], g_dsignet[9], g_dSr[1];
  long error=0;
  //initialising variables
  for (j=0; j<max_ngcomp; j++) {
     k1[j]=0;
     k2[j]=0;
     k3[j]=0;
     k4[j]=0;
     y1[j]=0;
     y2[j]=0;
     y3[j]=0;
     y_hat[j]=0;
     y_til[j]=0;
     y_err[j]=0;
    if(j<max_ngstrain) {
      strain_gen_aux[j]=0;
      dstrain_gen_aux[j]=0;
    }
  }
  //}
  //{ preparation
  // set initial substep length
  if (!dtsub) h = 1.0;
  else h = dtsub;
  if(err_sig >= 1) h = 1/err_sig; //Constant number of substeps
  // y1 has to be initialized to the values of stress vector sig and state variables
  garray_rcopy(stress_gen, 0, y1, 0, ngstress);
  garray_rcopy(qstatev, 0, y1, ngstress, nstatev);
  garray_copy(strain_gen, strain_gen_aux, ngstrain);
  //}
  //{ first call to calc_fsigq in a do loop before the iteration loop
  // calculate initial values of increments k1
  do {
    // dstrain_gen_aux = h*dstrain_gen
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);

    // calculate initial vector of coefficients for the first stage of Runge-Kutta method
    // herr = calculate_flag_2(strain_gen_aux, y1, y1+ngstress, dstrain_gen_aux, k1, k1+ngstress, parms, kinc);
    errorfsigq += correct_statev_values(strain_gen_aux, y1, y1+ngstress, dstrain_gen_aux, 4);
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);
    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq += calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k1+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k1, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))
        return 1;
    } else {
      // increase the number of model evaluations
      neval += 1;
      rkf_statev[0]=neval;
      break;
    }
  } while (h >= dtmin);
  //}
  //{ ITERATION LOOP
  for (k=0; k<ni && dt<1.0; k++) {
    //{ preparation
    // the number state variables stored after stress components in the vector y_hat that should be checked for tolerance
    stvcomp = 1;
    // strain_gen_aux = strain_gen + dt*dstrain_gen;
    garray_addmult(strain_gen, dstrain_gen, dt, strain_gen_aux, ngstrain);
    // dstrain_gen_aux = h*dstrain_gen;
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);
    if (h==0.0) return 2;
    // correct ascan value after integration by RKF
    aux = y1[ngstress+7]+k1[ngstress+7];
    // ascan value was cut because it was out of prescribed range <0;1>
    // do not check tolerance for state variables, i.e. the number of checked state variables is zero
    if (((fabs(1.0-aux) < 1.0e-7) && (fabs(1.0-y1[ngstress+7]) > 1.0e-7)) || ((fabs(aux) < 1.0e-7) && (y1[ngstress+7] > 1.0e-7))) stvcomp = 0;
    // correct em value after integration by RKF
    aux = y1[ngstress+4]+k1[ngstress+4];
    if (((fabs(0.01-aux) < 1.0e-7) && (fabs(y1[ngstress+7]-0.01) > 1.0e-7))) {
      // em value was cut on 0.01
      // use forward Euler rule for state variables + non-stress components of generalized stress vector
      y_hat[ngstress+4] = y_til[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
      // do not check tolerance for state variables, i.e. the number of checked state variables is zero
      stvcomp = 0;
    }
    // y2(i) = y1(i) + 1/2*h*k1(i)
    //    addmultv(y1, k1, 0.5*h, y2);
    garray_addmult (y1, k1, 0.5, y2, ngcomp);
    // strain_gen_aux = strain_gen + (dt+0.5*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+0.5*h, strain_gen_aux, ngstrain);
    //}
    //{ second call to calc_fsigq
    // calculate coefficients for the second stage of Runge-Kutta method
    errorfsigq += correct_statev_values(strain_gen_aux, y2, y2+ngstress, dstrain_gen_aux, 4);
    convert_from_general (y2, g_signet, g_othervar, (ngstress-6), 1);
    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq += calc_fsigq(g_signet, g_suction, g_Temper, y2+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k2+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k2, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else if (err_sig >= 1)               return 1; //Constant number of substeps
      else                              continue;
    }
    //}
    //{ preparation
    // y3(i) = y1(i) + 3/4*h*k2
    garray_addmult (y1, k2, 0.75, y3, ngcomp);
    if (stvcomp == 0) {
      y3[ngstress+7] = y1[ngstress+7]+0.75*k1[ngstress+7];
      y3[ngstress+4] = y1[ngstress+4]+0.75*k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+0.75*h)*dstrain_gen(ngstrain);
    garray_addmult (strain_gen, dstrain_gen, dt+0.75*h, strain_gen_aux, ngstrain);
    //}
    //{ third call to calc_fsigq
    // calculate coefficients for the third stage of Runge-Kutta method
    errorfsigq += correct_statev_values(strain_gen_aux, y3, y3+ngstress, dstrain_gen_aux, 4);
    convert_from_general (y3, g_signet, g_othervar, (ngstress-6), 1);
    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq += calc_fsigq(g_signet, g_suction, g_Temper, y3+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k3+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k3, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else if (err_sig >= 1)               return 1; //Constant number of substeps
      else                              continue;
    }
    //}
    //{ preparation
    // calculate third stage of embedded Runge-Kutta-Fehlberg method
    for (j=0; j<ngcomp; j++)
      //      y_hat(j) = y1(j)+h*((2.0/9.0)*k1(j) + (1.0/3.0)*k2(j) + (4.0/9.0)*k3(j));
      y_hat[j] = y1[j]+((2.0/9.0)*k1[j] + (1.0/3.0)*k2[j] + (4.0/9.0)*k3[j]);

    if (stvcomp == 0) {
      // use forward Euler rule for selected state variables
      y_hat[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y_hat[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    //}
    //{ fourth call to calc_fsigq
    // calculate coefficients k4 for the second stage of embedded Runge-Kutta-Fehlberg method
    errorfsigq += correct_statev_values(strain_gen_aux, y_hat, y_hat+ngstress, dstrain_gen_aux, 4);
    convert_from_general (y_hat, g_signet, g_othervar, (ngstress-6), 1);
    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq += calc_fsigq(g_signet, g_suction, g_Temper, y_hat+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k4+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k4, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else if (err_sig >= 1)               return 1; //Constant number of substeps
      else                              continue;
    }
    //}
    //{ finalization and error computation
    for (j=0; j<ngcomp; j++)
      //      y_til(j) = y1(j)+h*((7.0/24.0)*k1(j) + (1.0/4.0)*k2(j) + (1.0/3.0)*k3(j) + (1.0/8.0)*k4(j));
      y_til[j] = y1[j]+((7.0/24.0)*k1[j] + (1.0/4.0)*k2[j] + (1.0/3.0)*k3[j] + (1.0/8.0)*k4[j]);
    if (stvcomp == 0) {
      // use forward Euler rule for selected state variables
      y_til[ngstress+7] = y_hat[ngstress+7];
      y_til[ngstress+4] = y_hat[ngstress+4];
    }
    // increase the number of model evaluations
    neval += 3;
    rkf_statev[0] = neval;
    // solution error is difference between 3rd and 2nd satges
    garray_subtract(y_hat, y_til, y_err, ngcomp);
    // stress tensor norm stored in the Voigt notation, additional components in y_hat are ignored
    norm_sig = tensor_dbldot_prod(y_hat, y_hat, 2.0);
    norm_sig = sqrt(norm_sig);
    // calculate normalized errors of stress components
    aux = 1.0/norm_sig;
    for (j=0; j<6; j++)
      y_err[j] = fabs(y_err[j])*aux;
    // calculate normalized error of Sr (and eventually other additional variables)
    for(j=6; j<ngcomp; j++) {
      if (fabs(y_hat[j]) < sv_rkf_zero) y_err[j] = fabs(y_err[j]);
      else y_err[j] = fabs(y_err[j])/y_hat[j];
    }
    // calculate norm of residual
    norm_r = garray_snorm(y_err, 0, ngcomp);
    // calculate norm of stress components
    norm_sig = garray_snorm(y_err, 0, 6);
    // zero residual of suction
    y_err[ngstress+1] = 0.0;
    // zero residual of temperature
    y_err[ngstress+3] = 0.0;
    // zero residual of swelling indicator
    y_err[ngstress+nstatev-1] = 0.0;
    // calculate norm of remaining components of generalized stress vector and state variables
    norm_stv = garray_snorm(y_err, 6, stvcomp);
    // calculate optimum substep length
    if (norm_sig != 0.0) h_opt_sig = 0.9*h*pow(err_sig/norm_sig, 1.0/3.0);
    else h_opt_sig = 1.0;
    if (norm_stv != 0.0) h_opt_stv = 0.9*h*pow(err_stv/norm_stv, 1.0/3.0);
    else h_opt_stv = 1.0;
    h_opt = min2(h_opt_sig, h_opt_stv);
    //}
    //{ rejection or acceptance of substep
    //if (((norm_sig < err_sig) && (norm_stv < err_stv))) {
    if (((norm_sig < err_sig) && (norm_stv < err_stv)) || norm_sig>1) {
      // substep is accepted,
      // update y1 to the values of attained generalized stress vector and state variables
      garray_swap(y_hat, y1, ngcomp);
      // use k4 coefficient evaluated at y_hat as k1 in the next step
      garray_swap(k4, k1, ngcomp);
      // set the cumulative RK substep length
      dt += h;
      if (stopk) // dt = 1 => whole interval has been integrated
        break;
      // store original step length due to FSAL concept
      h_old = h;
      // select optimal RK substep length
      h = min2(4.0*h, h_opt);
      if (h > 1.0) h = 1.0;
      if(err_sig >= 1) h = 1/err_sig; //Constant number of substeps
      // store the last attained RK substep length to be the initial substep length in the next RKF method call
      dtsub = h;
      rkf_statev[3] = dtsub;
      // check the maximum RK substep length
      if (1.0-dt < h) {
        h = 1.0-dt;
        stopk = 1;
      }
      // recalculate actual increments in k1 due to FSAL concept
      garray_addmult(NULL, k1, h/h_old, k1, ngcomp);
    } else {
      // store original step length due to FSAL concept
      h_old = h;
      // substep is not accepted => reduce substep length
      h = max2(1.0/4.0*h, h_opt);
      // recalculate actual increments in k1 due to FSAL concept
      garray_addmult(NULL, k1, h/h_old, k1, ngcomp);
      // check for minimum step size
      if (h < dtmin) {
        dtsub = idtsub;
        rkf_statev[3] = dtsub;
        return 1;
      }
      if (mindt > h) mindt = h;
      rkf_statev[2] = mindt;
      // zero possibly set flag for the main loop breaking due to limit value dt
      stopk = 0;
    }
    //}
  }
  //} end of the iteration loop
  //{ finalization (storing results, adapting time step)
  if (dt > 1.0) dt = 1.0; // cut-off rounding error in the cumulative RK substep length
  if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
    // store the attained stresses - they are stored in y1 because the content of y_hat and y1 is swapped
    // in the case of accepted substep
    garray_rcopy(y1, 0, stress_gen, 0, ngstress);
    // state variables must be also integrated by Runge-Kutta => they must be copied from the end of y1 vector
    garray_rcopy(y1, ngstress, qstatev, 0, nstatev);
    nstep += k;
    rkf_statev[1] = nstep;
  } else {
    // RKF method cannot integrate the function with required error was in the given maximum number of substeps
    // dt must be less then 1.0
    dtsub = idtsub;
    rkf_statev[3] = dtsub;
    return 1;
  }
  //}
  return 0;
}

/**
  The function integrates ODEs by the kind of Runge-Kutta-Fehlberg 34 method.

  @param strain_gen[ngstrain] - (input): controlled variables [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]^T
  @param stress_gen[max_ngstress] - (input): implied variables, will not be updated [signet11, signet22, signet33, signet12, signet13, signet23, Sr]^T. "signet" is net stress, that is total stress minus krondelta x ua.
  @param qstatev[max_nstatev] - (input): state variable vector, will be updated only for flag==0.
  @param dstrain_gen[ngstrain] - (input): increment of controlled variables [deps11, deps22, deps33, dgamma12, dgamma13, dgamma23, ds, dT]^T
  @param dtime - (input): increment of time
  @param parms[c_nparms] - array of model parameters (input)
  @param rkf_statev[nrkf_statev] - (input): array of state variables from SIFEL for Runge-Kutta
  @param rkf_parms[nrkf_parms] - (input): array of parameters from SIFEL for Runge-Kutta

  @retval 0 - on success
  @retval 1 - in the case of substepping failure => global step length must be decreased

  Created by Tomas Koudelka, 6.10.2015; modified by G. Scaringi, January 2019 //@GS
*/
long General_model::rkf34(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[],
  double dtime, double parms[], double rkf_statev[], double rkf_parms[], int kinc) {
  //{ initialization of variables
  long errorfsigq = 0;
  long neval = rkf_statev[0];
  long nstep = rkf_statev[1];
  double mindt = rkf_statev[2];
  double dtsub = rkf_statev[3];
  double h, h_opt, h_opt_sig, h_opt_stv;
  double dt = 0.0;
  double err_sig = rkf_parms[0];
  double err_stv = rkf_parms[0];
  double norm_sig, norm_stv, /*norm_stvo,*/ aux, norm_r;
  double dtmin = rkf_parms[1];
  double idtsub = dtsub;
  double sv_rkf_zero = rkf_parms[3]; //needed
  long k, j;
  long ni = rkf_parms[2];
  long stopk = 0;
  long max_ngcomp = max_ngstress + max_nstatev;
  long ngcomp = ngstress + nstatev;
  long stvcomp;
  double k1[max_ngstress + max_nstatev], k2[max_ngstress + max_nstatev], k3[max_ngstress + max_nstatev], k4[max_ngstress + max_nstatev], k5[max_ngstress + max_nstatev];
  double y1[max_ngstress + max_nstatev], y2[max_ngstress + max_nstatev], y3[max_ngstress + max_nstatev], y4[max_ngstress + max_nstatev], y_hat[max_ngstress + max_nstatev], y_til[max_ngstress + max_nstatev], y_err[max_ngstress + max_nstatev];
  double strain_gen_aux[max_ngstrain], dstrain_gen_aux[max_ngstrain];
  long herr;
  // additional variables (for calc_fsigq)
  double g_suction=0, g_dsuction=0, g_Temper=0, g_dTemper=0, g_Sr=0, g_Sr_new=0;
  double g_signet[9], g_deps[9], g_ddum9[9], g_othervar[5], g_dsignet[9], g_dSr[1];
  long error=0;
    //initialising variables
  for (j=0; j<max_ngcomp; j++) {
     k1[j]=0;
     k2[j]=0;
     k3[j]=0;
     k4[j]=0;
     y1[j]=0;
     y2[j]=0;
     y3[j]=0;
     y4[j]=0;
     y_hat[j]=0;
     y_til[j]=0;
     y_err[j]=0;
    if(j<max_ngstrain) {
      strain_gen_aux[j]=0;
      dstrain_gen_aux[j]=0;
    }
  }
  //}
  //{ preparation
  // set initial substep length
  if (!dtsub) h = 1.0;
  else h = dtsub;
  // y1 has to be initialized to the values of stress vector sig and state variables
  garray_rcopy(stress_gen, 0, y1, 0, ngstress);
  garray_rcopy(qstatev, 0, y1, ngstress, nstatev);
  garray_copy(strain_gen, strain_gen_aux, ngstrain);
  //}
  //{ ITERATION LOOP
  for (k=0; k<ni && dt<1.0; k++) {
    //{ preparation
    // the number state variables stored after stress components in the vector y_hat that should be checked for tolerance
    stvcomp = 1;
    // strain_gen_aux = strain_gen + dt*dstrain_gen;
    garray_addmult(strain_gen, dstrain_gen, dt, strain_gen_aux, ngstrain);
    // dstrain_gen_aux = h*dstrain_gen;
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);
    if (h==0.0) return 2;
    //}
    //{ first call to calc_fsigq
    // calculate initial vector of coefficients for the first stage of Runge-Kutta method
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k1+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k1, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // correct ascan value after integration by RKF
    aux = y1[ngstress+7]+k1[ngstress+7];
    // ascan value was cut because it was out of prescribed range <0;1>
    // do not check tolerance for state variables, i.e. the number of checked state variables is zero
    if (((fabs(1.0-aux) < 1.0e-7) && (fabs(1.0-y1[ngstress+7]) > 1.0e-7)) || ((fabs(aux) < 1.0e-7) && (y1[ngstress+7] > 1.0e-7))) stvcomp = 0;
    // correct em value after integration by RKF
    aux = y1[ngstress+4]+k1[ngstress+4];
    if (((fabs(0.01-aux) < 1.0e-7) && (fabs(y1[ngstress+7]-0.01) > 1.0e-7))) {
      // em value was cut on 0.01
      // use forward Euler rule for state variables + non-stress components of generalized stress vector
      y_hat[ngstress+4] = y_til[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
      // do not check tolerance for state variables, i.e. the number of checked state variables is zero
      stvcomp = 0;
    }
    // find y2(i) = y1(i) + 1/4*h*k1(i)
    garray_addmult (y1, k1, 0.25, y2, ngcomp);
    // strain_gen_aux = strain_gen + (dt+0.5*h)*dstrain_gen(ngstrain);
    garray_addmult (strain_gen, dstrain_gen, dt+0.5*h, strain_gen_aux, ngstrain);
    //}
    //{ second call to calc_fsigq
    // calculate coefficients for the second stage of Runge-Kutta method
    convert_from_general (y2, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y2+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k2+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k2, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y3(i)
    for(j=0; j<ngcomp; j++) y3[j]=y1[j] + ((4.0/81.0)*k1[j] + (32.0/81.0)*k2[j]);
    if (stvcomp == 0) {
      y3[ngstress+7] = y1[ngstress+7]+(36.0/81.0)*k1[ngstress+7];
      y3[ngstress+4] = y1[ngstress+4]+(36.0/81.0)*k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+36/81*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+36.0/81.0*h, strain_gen_aux, ngstrain);
    //}
    //{ third call to calc_fsigq
    // calculate coefficients for the third stage of Runge-Kutta method
    convert_from_general (y3, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y3+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k3+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k3, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y4(i)
    for(j=0; j<ngcomp; j++) y4[j]=y1[j] + ((57.0/98.0)*k1[j] - (432.0/343.0)*k2[j] + (1053.0/686.0)*k3[j]);
    if (stvcomp == 0) {
      y4[ngstress+7] = y1[ngstress+7]+(6.0/7.0)*k1[ngstress+7];
      y4[ngstress+4] = y1[ngstress+4]+(6.0/7.0)*k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+6/7*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+6.0/7.0*h, strain_gen_aux, ngstrain);
    //}
    //{ fourth call to calc_fsigq
    // calculate coefficients for the fourth stage of Runge-Kutta method
    convert_from_general (y4, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y4+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k4+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k4, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // calculate fourth stage of embedded Runge-Kutta-Fehlberg method
    for (j=0; j<ngcomp; j++) y_til[j]=y1[j] + ((1.0/6.0)*k1[j] + (27.0/52.0)*k3[j] + (49.0/156.0)*k4[j]);
    if (stvcomp == 0) {
      // use forward Euler rule for selected state variables
      y_til[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y_til[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + h*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+h, strain_gen_aux, ngstrain);
    //}
    //{ fifth call to calc_fsigq
    // calculate coefficients for the fifth stage of Runge-Kutta method
    convert_from_general (y_til, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y_til+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k5+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k5, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ finalization and error computation
    for (j=0; j<ngcomp; j++) y_hat[j]=y1[j] + ((43.0/288.0)*k1[j] + (243.0/416.0)*k3[j] + (343.0/1872.0)*k4[j] + (1.0/12.0)*k5[j]);
    if (stvcomp == 0) {
      // use forward Euler rule for selected state variables
      y_hat[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y_hat[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    // increase the number of model evaluations
    neval += 5;
    rkf_statev[0] = neval;
    // solution error is difference between 3rd and 2nd satges
    garray_subtract(y_hat, y_til, y_err, ngcomp);
    // stress tensor norm stored in the Voigt notation, additional components in y_hat are ignored
    norm_sig = tensor_dbldot_prod(y_hat, y_hat, 2.0);
    norm_sig = sqrt(norm_sig);
    // calculate normalized errors of stress components
    aux = 1.0/norm_sig;
    for (j=0; j<6; j++)
      y_err[j] = fabs(y_err[j])*aux;
    // calculate normalized error of Sr and particular state variables
    for(j=6; j<ngcomp; j++) {
      if (fabs(y_hat[j]) < sv_rkf_zero) y_err[j] = fabs(y_err[j]);
      else y_err[j] = fabs(y_err[j])/y_hat[j];
    }
    // calculate norm of residual
    norm_r = garray_snorm(y_err, 0, ngcomp);
    // calculate norm of stress components
    norm_sig = garray_snorm(y_err, 0, 6);
    // zero residual of suction
    y_err[ngstress+1] = 0.0;
    // zero residual of temperature
    y_err[ngstress+3] = 0.0;
    // zero residual of swelling indicator
    y_err[ngstress+nstatev-1] = 0.0;
    // calculate norm of remaining components of generalized stress vector and state variables
    norm_stv = garray_snorm(y_err, 6, stvcomp);
    // calculate optimum substep length
    if (norm_sig != 0.0) h_opt_sig = 0.9*h*pow(err_sig/norm_sig, 1.0/4.0);
    else h_opt_sig = 1.0;
    if (norm_stv != 0.0) h_opt_stv = 0.9*h*pow(err_stv/norm_stv, 1.0/4.0);
    else h_opt_stv = 1.0;
    h_opt = min2(h_opt_sig, h_opt_stv);
    //}
    //{ acceptance or rejection of substep
    if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
      // substep is accepted,
      // update y1 to the values of attained generalized stress vector and state variables
      garray_swap(y_hat, y1, ngcomp);
      // set the cumulative substep length
      dt += h;
      if (stopk) // dt = 1 => whole interval has been integrated
        break;
      // select optimal RK substep length
      h = min2(4.0*h, h_opt);
      if (h > 1.0) h = 1.0;
      // store the last attained RK substep length to be the initial substep length in the next RKF method call
      dtsub = h;
      rkf_statev[3] = dtsub;
      // check the maximum RK substep length
      if (1.0-dt < h) {
        h = 1.0-dt;
        stopk = 1;
      }
    } else {
      // substep is not accepted => reduce substep length
      h = max2(1.0/4.0*h, h_opt);
      // check for minimum step size
      if (h < dtmin) {
        dtsub = idtsub;
        rkf_statev[3] = dtsub;
        return 1;
      }
      if (mindt > h) mindt = h;
      rkf_statev[2] = mindt;
      // zero possibly set flag for the main loop breaking due to limit value dt
      stopk = 0;
    }
    //}
  }
  //} end of the iteration loop
  //{ finalization (storing results, adapting time step)
  if (dt > 1.0) dt = 1.0; // cut-off rounding error in the cumulative RK substep length
  if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
    // store the attained stresses - they are stored in y1 because the content of y_hat and y1 is swapped
    // in the case of accepted substep
    garray_rcopy(y1, 0, stress_gen, 0, ngstress);
    // state variables must be also integrated by Runge-Kutta => they must be copied from the end of y1 vector
    garray_rcopy(y1, ngstress, qstatev, 0, nstatev);
    nstep += k;
    rkf_statev[1] = nstep;
  } else {
    // RKF method cannot integrate the function with required error was in the given maximum number of substeps
    // dt must be less then 1.0
    dtsub = idtsub;
    rkf_statev[3] = dtsub;
    return 1;
  }
  //}
  return 0;
}

/**
  The function integrates ODEs by the kind of Runge-Kutta-Fehlberg 45 method.

  @param strain_gen[ngstrain] - (input): controlled variables [eps11, eps22, eps33, gamma12, gamma13, gamma23, s, T]^T
  @param stress_gen[max_ngstress] - (input): implied variables, will not be updated [signet11, signet22, signet33, signet12, signet13, signet23, Sr]^T. "signet" is net stress, that is total stress minus krondelta x ua.
  @param qstatev[max_nstatev] - (input): state variable vector, will be updated only for flag==0.
  @param dstrain_gen[ngstrain] - (input): increment of controlled variables [deps11, deps22, deps33, dgamma12, dgamma13, dgamma23, ds, dT]^T
  @param dtime - (input): increment of time
  @param parms[c_nparms] - array of model parameters (input)
  @param rkf_statev[nrkf_statev] - (input): array of state variables from SIFEL for Runge-Kutta
  @param rkf_parms[nrkf_parms] - (input): array of parameters from SIFEL for Runge-Kutta

  @retval 0 - on success
  @retval 1 - in the case of substepping failure => global step length must be decreased

  Created by Tomas Koudelka, 6.10.2015; modified by G. Scaringi, January 2019 //@GS
*/
long General_model::rkf45(double strain_gen[], double stress_gen[], double qstatev[], double dstrain_gen[],
  double dtime, double parms[], double rkf_statev[], double rkf_parms[], int kinc) {
  //{ initialization of variables
  long errorfsigq = 0;
  long neval = rkf_statev[0];
  long nstep = rkf_statev[1];
  double mindt = rkf_statev[2];
  double dtsub = rkf_statev[3];
  double h, h_opt, h_opt_sig, h_opt_stv;
  double dt = 0.0;
  double err_sig = rkf_parms[0];
  double err_stv = rkf_parms[0];
  double norm_sig, norm_stv, aux, norm_r;
  double dtmin = rkf_parms[1];
  double idtsub = dtsub;
  double sv_rkf_zero = rkf_parms[3]; //needed
  long k, j;
  long ni = rkf_parms[2];
  long stopk = 0;
  long max_ngcomp = max_ngstress + max_nstatev;
  long ngcomp = ngstress + nstatev;
  long stvcomp;
  double k1[max_ngstress + max_nstatev], k2[max_ngstress + max_nstatev], k3[max_ngstress + max_nstatev], k4[max_ngstress + max_nstatev], k5[max_ngstress + max_nstatev], k6[max_ngstress + max_nstatev];
  double y1[max_ngstress + max_nstatev], y2[max_ngstress + max_nstatev], y3[max_ngstress + max_nstatev], y4[max_ngstress + max_nstatev], y5[max_ngstress + max_nstatev], y6[max_ngstress + max_nstatev];
  double y_hat[max_ngstress + max_nstatev], y_til[max_ngstress + max_nstatev], y_err[max_ngstress + max_nstatev];
  double strain_gen_aux[max_ngstrain], dstrain_gen_aux[max_ngstrain];
  long herr;
  // additional variables (for calc_fsigq)
  double g_suction=0, g_dsuction=0, g_Temper=0, g_dTemper=0, g_Sr=0, g_Sr_new=0;
  double g_signet[9], g_deps[9], g_ddum9[9], g_othervar[5], g_dsignet[9], g_dSr[1];
  long error=0;
      //initialising variables
  for (j=0; j<max_ngcomp; j++) {
     k1[j]=0;
     k2[j]=0;
     k3[j]=0;
     k4[j]=0;
     k5[j]=0;
     k6[j]=0;
     y1[j]=0;
     y2[j]=0;
     y3[j]=0;
     y4[j]=0;
     y5[j]=0;
     y6[j]=0;
     y_hat[j]=0;
     y_til[j]=0;
     y_err[j]=0;
    if(j<max_ngstrain) {
      strain_gen_aux[j]=0;
      dstrain_gen_aux[j]=0;
    }
  }

  //}
  //{ preparation
  // set initial substep length
  if (!dtsub) h = 1.0;
  else h = dtsub;
  // y1 has to be initialized to the values of stress vector sig and state variables
  garray_rcopy(stress_gen, 0, y1, 0, ngstress);
  garray_rcopy(qstatev, 0, y1, ngstress, nstatev);
  garray_copy(strain_gen, strain_gen_aux, ngstrain);
  //}
  //{ ITERATION LOOP
  for (k=0; k<ni && dt<1.0; k++) {
    //{ preparation
    // the number state variables stored after stress components in the vector y_hat that should be checked for tolerance
    stvcomp = 1;
    // strain_gen_aux = strain_gen + dt*dstrain_gen;
    garray_addmult(strain_gen, dstrain_gen, dt, strain_gen_aux, ngstrain);
    // dstrain_gen_aux = h*dstrain_gen;
    garray_addmult(NULL, dstrain_gen, h, dstrain_gen_aux, ngstrain);
    if (h==0.0)
      return 2;
    //}
    //{ first call to calc_fsigq
    // calculate initial vector of coefficients for the first stage of Runge-Kutta method
    convert_from_general (y1, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y1+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k1+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k1, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // correct ascan value after integration by RKF
    aux = y1[ngstress+7]+k1[ngstress+7];
    // ascan value was cut because it was out of prescribed range <0;1>
    // do not check tolerance for state variables, i.e. the number of checked state variables is zero
    if (((fabs(1.0-aux) < 1.0e-7) && (fabs(1.0-y1[ngstress+7]) > 1.0e-7)) || ((fabs(aux) < 1.0e-7) && (y1[ngstress+7] > 1.0e-7))) stvcomp = 0;
    // correct em value after integration by RKF
    aux = y1[ngstress+4]+k1[ngstress+4];
    if (((fabs(0.01-aux) < 1.0e-7) && (fabs(y1[ngstress+7]-0.01) > 1.0e-7))) {
      // em value was cut on 0.01
      // use forward Euler rule for state variables + non-stress components of generalized stress vector
      y_hat[ngstress+4] = y_til[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
      // do not check tolerance for state variables, i.e. the number of checked state variables is zero
      stvcomp = 0;
    }
    // y2(i) = y1(i) + 1/4*k1(i)
    garray_addmult (y1, k1, 0.25, y2, ngcomp);
    // strain_gen_aux = strain_gen + (dt+0.25*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+0.25*h, strain_gen_aux, ngstrain);
    //}
    //{ second call to calc_fsigq
    // calculate coefficients for the second stage of Runge-Kutta method
    convert_from_general (y2, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y2+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k2+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k2, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y3(i)
    for(j=0; j<ngcomp; j++) y3[j]=y1[j] + h*((3.0/32.0)*k1[j] + (9.0/32.0)*h*k2[j]);
    if (stvcomp == 0) {
      y3[ngstress+7] = y1[ngstress+7]+(12.0/32.0)*k1[ngstress+7];
      y3[ngstress+4] = y1[ngstress+4]+(12.0/32.0)*k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+12/32*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+12.0/32.0*h, strain_gen_aux, ngstrain);
    //}
    //{ third call to calc_fsigq
    // calculate coefficients for the third stage of Runge-Kutta method
    convert_from_general (y3, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y3+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k3+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k3, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y4(i)
    for(j=0; j<ngcomp; j++) y4[j]=y1[j] + h*((1932.0/2197.0)*k1[j] - (7200.0/2197.0)*k2[j] + (7296.0/2197.0)*k3[j]);
    if (stvcomp == 0) {
      y4[ngstress+7] = y1[ngstress+7]+(2028.0/2197.0)*k1[ngstress+7];
      y4[ngstress+4] = y1[ngstress+4]+(2028.0/2197.0)*k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+2028/2197*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+2028.0/2197.0*h, strain_gen_aux, ngstrain);
    //}
    //{ fourth call to calc_fsigq
    // calculate coefficients for the fourth stage of Runge-Kutta method
    convert_from_general (y4, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y4+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k4+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k4, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y5(i)
    for(j=0; j<ngcomp; j++) y5[j]=y1[j] + h*((439.0/216.0)*k1[j] - 8.0*k2[j] + (3680.0/513.0)*k3[j] - (845.0/4104.0)*k4[j]);
    if (stvcomp == 0) {
      y5[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y5[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+h, strain_gen_aux, ngstrain);
    //}
    //{ fifth call to calc_fsigq
    // calculate coefficients for the fifth stage of Runge-Kutta method
    convert_from_general (y5, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y5+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k5+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k5, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ preparation
    // find y6(i)
    for(j=0; j<ngcomp; j++) y6[j]=y1[j] + h*((-8.0/27.0)*k1[j] + 2.0*k2[j] - (3544.0/2565.0)*k3[j] + (1859.0/4104.0)*k4[j] - (11.0/40.0)*k5[j]);
    if (stvcomp == 0) {
      y6[ngstress+7] = y1[ngstress+7]+0.5*k1[ngstress+7];
      y6[ngstress+4] = y1[ngstress+4]+0.5*k1[ngstress+4];
    }
    // strain_gen_aux = strain_gen + (dt+0.5*h)*dstrain_gen(ngstrain);
    garray_addmult(strain_gen, dstrain_gen, dt+0.5*h, strain_gen_aux, ngstrain);
    //}
    //{ sixth call to calc_fsigq
    // calculate coefficients for the sixth stage of Runge-Kutta method
    convert_from_general (y6, g_signet, g_othervar, (ngstress-6), 1);

    g_Sr = g_othervar[0];
    convert_from_general (strain_gen_aux, g_ddum9, g_othervar, (ngstrain-6), 2);
    g_suction=g_othervar[0];
    g_Temper=g_othervar[1];
    for (int i=0; i<6; i++) if (gscalar_dabs(dstrain_gen_aux[i])>1) dstrain_gen_aux[i]=0;
    convert_from_general(dstrain_gen_aux, g_deps, g_othervar, (ngstrain-6), 2);
    g_dsuction=g_othervar[0];
    g_dTemper=g_othervar[1];

    errorfsigq = calc_fsigq(g_signet, g_suction, g_Temper, y6+ngstress, g_deps, g_dsuction, g_dTemper, dtime*h, g_dsignet, g_Sr_new, k6+ngstress, kinc);

    g_dSr[0] = g_Sr_new - g_Sr;
    convert_to_general(k6, g_dsignet, g_dSr,(ngstress-6),1);

    if (check_math_error() || errorfsigq){
      if (rkf_redstep(0.5, h, dtmin))   return 1;
      else                              continue;
    }
    //}
    //{ finalization and error computation
    // calculate sixth stage of embedded Runge-Kutta-Fehlberg method
    for (j=0; j<ngcomp; j++) {
      y_til[j]=y1[j] + h*((25.0/216.0)*k1[j] + (1408.0/2565.0)*k3[j] + (2197.0/4104.0)*k4[j] - 0.2*k5[j]);
      y_hat[j]=y1[j] + h*((16.0/135.0)*k1[j] + (6656.0/12825.0)*k3[j] + (28561.0/56430.0)*k4[j] - (9.0/50.0)*k5[j] + (2.0/55.0)*k6[j]);
    }
    if (stvcomp == 0) {
      // use forward Euler rule for selected state variables
      y_hat[ngstress+7] = y_til[ngstress+7] = y1[ngstress+7]+k1[ngstress+7];
      y_hat[ngstress+4] = y_til[ngstress+4] = y1[ngstress+4]+k1[ngstress+4];
    }
    // increase the number of model evaluations
    neval += 6;
    rkf_statev[0] = neval;
    // solution error is difference between 3rd and 2nd satges
    garray_subtract(y_hat, y_til, y_err, ngcomp);
    // stress tensor norm stored in the Voigt notation, additional components in y_hat are ignored
    norm_sig = tensor_dbldot_prod(y_hat, y_hat, 2.0);
    norm_sig = sqrt(norm_sig);
    // calculate normalized errors of stress components
    aux = 1.0/norm_sig;
    for (j=0; j<6; j++)  y_err[j] = fabs(y_err[j])*aux;
    // calculate normalized error of Sr and particular state variables
    for(j=6; j<ngcomp; j++) {
      if (fabs(y_hat[j]) < sv_rkf_zero) y_err[j] = fabs(y_err[j]);
      else y_err[j] = fabs(y_err[j])/y_hat[j];
    }
    // calculate norm of residual
    norm_r = garray_snorm(y_err, 0, ngcomp);
    // calculate norm of stress components
    norm_sig = garray_snorm(y_err, 0, 6);
    // zero residual of suction
    y_err[ngstress+1] = 0.0;
    // zero residual of temperature
    y_err[ngstress+3] = 0.0;
    // zero residual of swelling indicator
    y_err[ngstress+nstatev-1] = 0.0;
    // calculate norm of remaining components of generalized stress vector and state variables
    norm_stv = garray_snorm(y_err, 6, stvcomp);
    // calculate optimum substep length
    if (norm_sig != 0.0) h_opt_sig = 0.9*h*pow(err_sig/norm_sig, 0.2);
    else h_opt_sig = 1.0;
    if (norm_stv != 0.0) h_opt_stv = 0.9*h*pow(err_stv/norm_stv, 0.2);
    else h_opt_stv = 1.0;
    h_opt = min2(h_opt_sig, h_opt_stv);
    //}
    //{ acceptance or rejection of substep
    if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
      // substep is accepted,
      // update y1 to the values of attained generalized stress vector and state variables
      garray_swap(y_hat, y1, ngcomp);
      // set the cumulative substep length
      dt += h;
      if (stopk) // dt = 1 => whole interval has been integrated
        break;
      // select optimal RK substep length
      h = min2(4.0*h, h_opt);
      if (h > 1.0) h = 1.0;
      // store the last attained RK substep length to be the initial substep length in the next RKF method call
      dtsub = h;
      rkf_statev[3] = dtsub;
      if (1.0-dt < h) {
        h = 1.0-dt;
        stopk = 1;
      }
    } else {
      // substep is not accepted => reduce substep length
      h = max2(1.0/4.0*h, h_opt);
      // check for minimum step size
      if (h < dtmin) {
        dtsub = idtsub;
        rkf_statev[3] = dtsub;
        return 1;
      }
      if (mindt > h) mindt = h;
      rkf_statev[2] = mindt;
      // zero possibly set flag for the main loop breaking due to limit value dt
      stopk = 0;
    }
    //}
  }
  //} end of the iteration loop
  //{ finalization (storing results, adapting time step)
  if (dt > 1.0) dt = 1.0; // cut-off rounding error in the cumulative RK substep length
  if ((norm_sig < err_sig) && (norm_stv < err_stv)) {
    // store the attained stresses - they are stored in y1 because the content of y_hat and y1 is swapped
    // in the case of accepted substep
    garray_rcopy(y1, 0, stress_gen, 0, ngstress);
    // state variables must be also integrated by Runge-Kutta => they must be copied from the end of y1 vector
    garray_rcopy(y1, ngstress, qstatev, 0, nstatev);
    nstep += k;
    rkf_statev[1] = nstep;
  } else {
    // RKF method cannot integrate the function with required error was in the given maximum number of substeps
    // dt must be less than 1.0
    dtsub = idtsub;
    rkf_statev[3] = dtsub;
    return 1;
  }
  //}
  return 0;
}

//@GS more routines
/**
  The function reduces step length in Runge-Kutta-Fehlberg methods according to given
  coefficient and minimum step length.

  @param rc   - reduction coefficient of step length (input)
  @param h    - actual/modified step length (input/output)
  @param hmin - minimum step length (input)

  @retval 0 - On sucessfull step reduction.
  @retval 1 - The minimum step length has been attained and the step cannot be reduced further.

  Created by Tomas Koudelka, 24.5.2016; modified by G. Scaringi, January 2019 //@GS
*/
long General_model::rkf_redstep(double rc, double &h, double hmin) {
  if (h > hmin)
  {
    h *= rc;
    if (h < hmin)
      h = hmin;
    return 0;
  }
  else
    return 1;
}

double General_model::tensor_dbldot_prod (double a[], double b[], double k) {
  long i;
  double nor;

  nor=0.0;
  for (i=0;i<3;i++){
    nor+=a[i]*b[i];
  }
  for (i=3;i<6;i++){
    nor+=k*a[i]*b[i];
  }
  return nor;
}

int General_model::garray_copy (double source[], double destination[], int length) {
  for (int i=0; i<length; i++) destination[i] = source[i];
  return(0);
}

int General_model::garray_rcopy (double source[], int source_index, double destination[], int destination_index, int number_of_elements) {
  for (int i=0; i<number_of_elements; i++) destination[destination_index+i] = source[source_index+i];
  return(0);
}

void General_model::garray_subtract (double a[], double b[], double c[], long int n) {
  long int i=0;
  for (i=0; i<n; i++) c[i] = a[i] - b[i];
}

int General_model::garray_addmult (double source1[], double source2[], double multiplier, double destination[], int length) {
  if (source1==NULL) for (int i=0; i<length; i++) destination[i] = source2[i]*multiplier;
  else for (int i=0; i<length; i++) destination[i] = source1[i] + source2[i]*multiplier;
  return(0);
}

int General_model::garray_swap (double a[], double b[], int length) {
  for (int i=0; i<length; i++) {
    double tmp = a[i];
    a[i] = b[i];
    b[i] = tmp;
  }
  return(0);
}

double General_model::garray_snorm(double a[], long fi, long nc) {
  double s = 0.0;
  long li = fi+nc;
  for (long i=fi; i < li; i++)
    s += a[i] * a[i];
  return sqrt(s);
}

long Hypoplasti_unsat_expansive_thermal::correct_statev_values(double strain_gen[], double stress_gen[],
  double qstatev[], double dstrain_gen[],  int call) {

  //BEGIN tension cutoff
  double othervar[5]={0,0,0,0,0};
  double ddum9[9]={0,0,0,0,0,0,0,0,0};
  convert_from_general (strain_gen, ddum9, othervar, ngstrain-6, 2);
  double suction=othervar[0];
  double Temper=othervar[1];

  double signet[9]={0,0,0,0,0,0,0,0,0};
  convert_from_general (stress_gen, signet, othervar, ngstress-6, 1);

  double sigef[9]={0,0,0,0,0,0,0,0,0};
  double SrM=qstatev[6];
  make_sig_unsat(signet, suction, SrM, sigef, suction);
  double chixsuction=sigef[0]-signet[0];

  double signet_corrected[9]={0,0,0,0,0,0,0,0,0};
  garray_move(signet, signet_corrected, 9);

  bool printstat=false;
  if(sigef[0]>0) {
      printstat=true;
      signet_corrected[0]=chixsuction;
  }
  if(sigef[4]>0) {
      printstat=true;
      signet_corrected[4]=chixsuction;
  }
  if(sigef[8]>0) {
      printstat=true;
      signet_corrected[8]=chixsuction;
  }
  if((sigef[0]+sigef[4]+sigef[8])/3>0) {
      printstat=true;
  }
  long errormath = 0;
  if(debug && printstat && call != 4) {

    if (call == 3) errormath++;
    cout<<"Correcting signet, call "<<call<<", signet[0]="<<signet[0]<<", signet[4]="<<signet[4]<<", signet[8]="<<signet[8]<<", sigef[0]="<<sigef[0]<<", sigef[4]="<<sigef[4]<<", sigef[8]="<<sigef[8]<<", suction="<<suction<<", dsuction="<<dstrain_gen[6]<<", SrM="<<SrM<<endl;
    garray_move(signet_corrected, signet, 9);
  }

  convert_to_general(stress_gen, signet, othervar, ngstress-6, 1);
  //END tension cutoff

  // SrM cutoff
  if (SrM < 0.0)
    SrM = 0.0;
  if (SrM > 1.0)
    SrM = 1.0;
  qstatev[6] = SrM;
  // Sr cutoff
  double Sr = stress_gen[ngstress-1];
  if (Sr < 0.0)
    Sr = 0.0;
  if (Sr > 1.0)
    Sr = 1.0;
  stress_gen[ngstress-1] = qstatev[2] = Sr;
  

  for(int i=0; i<ngstrain; i++) {if(isnan( strain_gen[i] )) errormath++;}
  for(int i=0; i<ngstrain; i++) {if(isnan( dstrain_gen[i] )) errormath++;}
  for(int i=0; i<ngstress; i++) {if(isnan( stress_gen[i] )) errormath++;}
  for(int i=0; i<nstatev; i++) {if(isnan( qstatev[i] )) errormath++;}

  return(errormath);
}

long Hypoplasti_unsat_expansive_thermal::correct_DDtan_gen(double strain_gen[], double stress_gen[], double *DDtan_gen) {
  //Check that dSr/ds is not negative, which is not physical
  if(strain_gen[6]>sairentry0 && DDtan_gen[54]<0) DDtan_gen[54]=0;
  //check NAN error of stiffness matrix
  long errormath = 0;
  for(int j=0; j<ngstress*ngstrain; j++) if(isnan(DDtan_gen[j])) errormath += 1;
  return (errormath);
};

long General_model::check_math_error() {
  long matherror=test_math_err();
  errno = 0;
  if(matherror) return(matherror);
  else return(0);
}
