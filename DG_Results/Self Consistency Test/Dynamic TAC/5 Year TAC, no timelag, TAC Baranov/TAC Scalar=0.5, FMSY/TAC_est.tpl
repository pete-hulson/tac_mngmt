

GLOBALS_SECTION
  #include "admodel.h"
  #define EOUT(var) cout <<#var<<" "<<var<<endl
TOP_OF_MAIN_SECTION
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);
DATA_SECTION
  init_int nages
  init_int nyrs_SIM
  init_int nyrs
  init_int nyrs_landings
  init_int nyrs_survey  
  init_int nyrs_TAC_SIM

  init_int phase_dummy
  init_int phase_beta1
  init_int phase_beta2
  init_int phase_beta1_index
  init_int phase_beta2_index
  init_int phase_ln_sel_age
  init_int phase_ln_sel_age_index
  init_int phase_F
  init_int phase_F_devs
  init_int phase_rec_dev
  init_int phase_ln_R_ave
  init_int phase_h
  init_int phase_ln_R0_offset
  init_int phase_ln_init_abund
  init_int phase_q

  init_number beta1_initial
  init_number beta2_initial
  init_number beta1_index_initial
  init_number beta2_index_initial
  init_number ln_sel_age_initial
  init_number ln_sel_age_index_initial
  init_number ln_F_initial
  init_number ln_F_devs_initial
  init_number ln_rec_devs_initial
  init_number ln_R_ave_initial
  init_number ln_h_initial
  init_number ln_R0_offset_initial
  init_number ln_init_abund_initial
  init_number ln_q_initial

  init_int SIM_REC_FXN_type
  init_int ASS_REC_FXN_type
  init_int est_steep_switch //==0 fix at sim_h, ==1 estimate
  init_int init_abund_switch  //how to estimate initial abundance, ==0 assume equilibrium age strcuture from an offset from R0, ==1 estimate initial abundance at each age
  init_int SIM_survey_sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int survey_sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int SIM_sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int TAC_CNST_switch //TAC switch ==1 TAC is set constant to level in yr_TAC_start, ==0 TAC is recalculated each year
  init_int TAC_timelag   //==0 no timelage, otherwise a 1 year lag was implemented in SIM to determine TAC (i.e., based on previous year abundance)
  init_int TAC_type  //==0 simple TAC=FMSY*exploitable Biomass, ==1 TAC=Baranov' catch equation fishing at FMSY
  init_int Fit_OBS  //==0 yes fit observed data, ==1 no fit true values (i.e., perform self-consistency run)  
  init_int F_devs_pen_switch //should include penalty on F_devs to avoid extreme values?  ==1 TRUE, ==0 FALSE
  init_int R_devs_pen_switch //should include penalty on R_devs to avoid extreme values?  ==1 TRUE, ==0 FALSE
  init_int F_devs_pen_mult
  init_int R_devs_pen_mult
  
  init_vector wt_at_age(1,nages)
  init_vector maturity(1,nages)
  init_number sex_scalar
  init_vector M(1,nages)
  init_number tspawn

  init_number sigma_recruit
  init_number ASS_sigma_recruit
  init_number sigma_F
  init_number ncatch
  init_number nindex
  init_vector SIM_sigma_index(1,nyrs_survey)
  init_vector sigma_index(1,nyrs_survey)
  init_vector sigma_landings(1,nyrs_landings)

  init_number wt_index
  init_number wt_recruit
  init_number wt_index_prop
  init_number wt_landings
  init_number wt_catch_prop
  
  init_number SIM_h
  init_number SIM_R_ave
  init_number SIM_SSB_zero
  init_number Fmsy
  init_number Fmsy_scalar
  init_number TAC_scalar
  init_number SIM_beta1_index
  init_number SIM_beta2_index
  init_number SIM_beta1
  init_number SIM_beta2
  init_number SIM_q
  init_vector SIM_F(1,nyrs_landings)
  init_matrix SIM_TAC_age(1,nyrs_TAC_SIM+1,1,nages) //TAC only enacted for nyrs_TAC, but calculated for nyrs_TAC+1 because of assumed one year timelag in implementation
  init_vector SIM_TAC(1,nyrs_TAC_SIM+1)
  init_matrix SIM_F_at_age(1,nyrs_landings,1,nages)
  init_matrix SIM_selectivity(1,nyrs_landings,1,nages)
  init_matrix SIM_selectivity_index(1,nyrs_survey,1,nages)
  init_vector SIM_ssb(1,nyrs_SIM)
  init_vector SIM_SPR_yr(1,nyrs_SIM)
  init_vector SIM_BH_REC(1,nyrs_SIM-1)
  init_vector SIM_recruits(1,nyrs_SIM) 
  init_matrix SIM_abundance_at_age(1,nyrs_SIM,1,nages) //before movement and mortality
  init_matrix SIM_abundance_spawn(1,nyrs_SIM,1,nages)
  init_vector SIM_init_abund(1,nages)
  init_vector SIM_biomass(1,nyrs_SIM)
  init_vector SIM_rec_devs(1,nyrs_SIM-1)  //don't have deviation in first year because assume at equilibrium

  init_matrix TRUE_catch_prop(1,nyrs_landings,1,nages)
  init_matrix OBS_catch_prop(1,nyrs_landings,1,nages)
  init_vector TRUE_landings(1,nyrs_landings)
  init_vector OBS_landings(1,nyrs_landings)
  init_vector TRUE_Index(1,nyrs_survey)
  init_vector OBS_Index(1,nyrs_survey)
  init_matrix TRUE_Index_prop(1,nyrs_survey,1,nages)
  init_matrix OBS_index_prop(1,nyrs_survey,1,nages)

  init_number MNconst
  init_int debug
  number t_survey
  
  int a
  int y
  int z
  int k
  int j
  int i
  int s
  int r
  int n
  int w
  int p
  int v
  int x
  int u
  int d

 !! cout << "input read" << endl;
 !! cout << "If Debug=-999 then input read correctly" << endl;
 !! cout << "Debug=" << endl;
 !! cout << debug << endl;


PARAMETER_SECTION
  !! cout << "reading parameters"<<endl;
  init_number dummy(phase_dummy)
  init_bounded_number beta1(0.00000001,10.,phase_beta1)  //logisitic fishery selectivity parameters
  init_bounded_number beta2(0.00000001,10.,phase_beta2)
  init_bounded_number beta1_index(0.00000001,10.,phase_beta1_index)  //logisitic fishery selectivity parameters
  init_bounded_number beta2_index(0.00000001,10.,phase_beta2_index)  
  init_bounded_vector ln_sel_age(1,nages,-100,50.,phase_ln_sel_age)  
  init_bounded_vector ln_sel_age_index(1,nages,-100,50.,phase_ln_sel_age_index)
  matrix selectivity(1,nyrs_landings,1,nages)
  matrix selectivity_index(1,nyrs_survey,1,nages)
  
  init_bounded_number ln_F(-50,50,phase_F)  //F
  number F
  init_bounded_dev_vector ln_F_devs(1,nyrs_landings,-20,20,phase_F_devs) //yearly F devs
  sdreport_vector fmult(1,nyrs_landings)  //total F 
  matrix F_at_age(1,nyrs,1,nages)
  sdreport_vector F_devs(1,nyrs_landings)

  init_bounded_dev_vector ln_rec_dev(1,nyrs-1,-20,20,phase_rec_dev)
  sdreport_vector rec_devs(1,nyrs-1)
  init_bounded_number ln_R_ave(1,30,phase_ln_R_ave)  //beverton-holt s-r parameters
  init_bounded_number ln_h(-30,30,phase_h)
  sdreport_number h
  sdreport_number R_ave
  init_bounded_number ln_R0_offset(-10,0,phase_ln_R0_offset)  //estimated offset of recruitment in first year from R0 (if assuming not starting from virgin conditions); only used if init_abund_switch==0
  sdreport_number R0_offset

  init_bounded_vector ln_init_abund(1,nages,1,30,phase_ln_init_abund)  //vector of initial abundance
  sdreport_vector init_abund(1,nages)
  
  init_bounded_number ln_q(-40,100.,phase_q) //catchability coeff.
  sdreport_number q

  vector SPR_N(1,nages)
  vector SPR_SSB(1,nages)
  number SPR
  number SSB_zero
  
  sdreport_vector recruits(1,nyrs) 
  vector BH_REC(1,nyrs-1)
  sdreport_vector ssb(1,nyrs)
  vector SPR_yr(1,nyrs)
  matrix ssb_temp(1,nyrs,1,nages)
  matrix abundance_at_age(1,nyrs,1,nages) //before movement and mortality
  matrix abundance_spawn(1,nyrs,1,nages)
  vector wt_mat_mult(1,nages)
  vector biomass(1,nyrs)
  matrix biomass_temp(1,nyrs,1,nages)

  matrix pred_catch_at_age(1,nyrs_landings,1,nages)
  vector catchsum(1,nyrs_landings)
  matrix pred_catch_prop(1,nyrs_landings,1,nages)
  vector pred_landings(1,nyrs_landings)
  matrix pred_landings_temp(1,nyrs,1,nages)

  matrix Index_age(1,nyrs_survey,1,nages)
  vector Indexsum(1,nyrs_survey)
  matrix pred_Index_prop(1,nyrs_survey,1,nages)
  vector pred_Index(1,nyrs_survey)

  vector res_recruit(1,nyrs-1)
  number rss_recruit
  vector res_landings(1,nyrs_landings)
  matrix res_catch_prop(1,nyrs_landings,1,nages)
  vector res_index(1,nyrs_survey)
  matrix res_index_prop(1,nyrs_survey,1,nages)
  number rss_landings
  vector rss_landings_temp(1,nyrs_landings)
  number rss_catch_prop
  vector rss_catch_prop_temp(1,nyrs_landings)
  number rss_index
  vector rss_index_temp(1,nyrs_survey)
  number rss_index_prop
  vector rss_index_prop_temp(1,nyrs_survey)
  number F_devs_pen
  number R_devs_pen
  
  number weighted_catch_prop
  number weighted_index
  number weighted_index_prop
  number weighted_landings
  number weighted_recruit

  objective_function_value f
  
  !! cout << "parameters set" << endl;
  
INITIALIZATION_SECTION  //set initial values
  beta1 beta1_initial;
  beta2 beta2_initial;
  beta1_index beta1_index_initial;
  beta2_index beta2_index_initial;
  ln_sel_age ln_sel_age_initial;
  ln_sel_age_index ln_sel_age_index_initial;
  ln_F ln_F_initial;
  ln_F_devs ln_F_devs_initial;
  ln_rec_dev ln_rec_devs_initial;
  ln_R_ave ln_R_ave_initial;
  ln_h ln_h_initial;
  ln_R0_offset ln_R0_offset_initial;
  ln_init_abund ln_init_abund_initial;
  ln_q ln_q_initial;
  
PROCEDURE_SECTION

  get_survey_selectivity();

  get_selectivity();

  get_abundance();

  get_biomass();

  get_catch();

  get_landings();

  //get_fish_dep_indices();

  get_fish_ind_indices();
  
  evaluate_the_objective_function();

FUNCTION get_survey_selectivity
 for (int a=1;a<=nages;a++)
  {
    for (int y=1;y<=nyrs_survey;y++)
        {
      if(survey_sel_switch==1)
       {
         selectivity_index(y,a)=mfexp(ln_sel_age_index(a))/(mfexp(ln_sel_age_index(a))+1);
        }
      if(survey_sel_switch==2)
       {     
         selectivity_index(y,a)=1/(1+mfexp(-log(19)*(a-(beta1_index))/(beta2_index)));
        }        
     }
   }

FUNCTION get_selectivity
  F=mfexp(ln_F);
  F_devs=mfexp(ln_F_devs);

  for (int a=1;a<=nages;a++)
   {
    for (int y=1;y<=nyrs_landings;y++)
      {
               if(sel_switch==1)
                 {
                 selectivity(y,a)=mfexp(ln_sel_age(a))/(mfexp(ln_sel_age(a))+1);
                 fmult(y)=(F)*mfexp(ln_F_devs(y));
                 }
               if(sel_switch==2)
                 {
                 selectivity(y,a)=1/(1+mfexp(-log(19)*(a-(beta1))/(beta2)));
                 fmult(y)=(F)*mfexp(ln_F_devs(y));                 
                 }
           F_at_age(y,a)=fmult(y)*selectivity(y,a);                
        }
     }                         
FUNCTION get_abundance

  if(est_steep_switch==0)
   {
    h=SIM_h;
   }
  if(est_steep_switch==1)
   {
    h=(mfexp(ln_h))/(mfexp(ln_h)+1);
   }
  R_ave=mfexp(ln_R_ave);
  R0_offset=mfexp(ln_R0_offset);
  
  for (int n=1;n<=nyrs-1;n++)
   {
    rec_devs(n)=mfexp(ln_rec_dev(n)-0.5*square(ASS_sigma_recruit));
   }

    for (int a=1;a<=nages;a++)
    {
      if(a==1)
        {
         SPR_N(a)=1000;
         wt_mat_mult(a)=sex_scalar*maturity(a)*wt_at_age(a);
         SPR_SSB(a)=wt_mat_mult(a)*SPR_N(a);
        }
      if(a>1)
        {
         SPR_N(a)=SPR_N(a-1)*mfexp(-M(a-1));
         wt_mat_mult(a)=sex_scalar*maturity(a)*wt_at_age(a);
         SPR_SSB(a)=wt_mat_mult(a)*SPR_N(a);
        }
      if(a==nages)
        {
         SPR_N(a)=SPR_N(a-1)*mfexp(-M(a-1))*(1/(1-mfexp(-M(a))));
         wt_mat_mult(a)=sex_scalar*maturity(a)*wt_at_age(a);         
         SPR_SSB(a)=wt_mat_mult(a)*SPR_N(a);
        }
    }
  SPR=sum(SPR_SSB)/1000;
  SSB_zero=R_ave*SPR;
  
 for (int y=1;y<=nyrs;y++)
  {
    for (int a=1;a<=nages;a++)
    {
      if (y==1)
      {
      if(init_abund_switch==0)  //assume init abundance in equilibrium based on an offset from R0
       {
    //equilibrium abundance in first year
        if(a==1) //age-0 recruits
        {
          recruits(y)=R_ave*mfexp(ln_R0_offset); 
          abundance_at_age(y,a)=recruits(y); 
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
          init_abund(a)=abundance_at_age(y,a);
        }
        if(a==2)
        {
          abundance_at_age(y,a)=abundance_at_age(y,a-1)*mfexp(-M(a-1)*(1-tspawn));
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
          init_abund(a)=abundance_at_age(y,a);
        }        
        if(a>2 && a!=nages)
        {
          abundance_at_age(y,a)=abundance_at_age(y,a-1)*mfexp(-M(a-1));
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
          init_abund(a)=abundance_at_age(y,a);
        }
        if(a==nages)
        {
         abundance_at_age(y,a)=abundance_at_age(y,a-1)*mfexp(-M(a-1))/(1-mfexp(-M(a)));
         abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
         init_abund(a)=abundance_at_age(y,a);
        }
        ssb_temp(y,a)=abundance_spawn(y,a)*wt_mat_mult(a);
        ssb(y)=sum(ssb_temp(y));
        SPR_yr(y)=ssb(y)/SSB_zero;
       }
      if(init_abund_switch==1)  //estimate initial abundance
       {
    //equilibrium abundance in first year
        if(a==1) //age-0 recruits
        {
          recruits(y)=mfexp(ln_init_abund(a)); 
          abundance_at_age(y,a)=recruits(y); 
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
          init_abund(a)=abundance_at_age(y,a);
        }
        if(a>1)
        {
          abundance_at_age(y,a)=mfexp(ln_init_abund(a));
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
          init_abund(a)=abundance_at_age(y,a);
        }
       ssb_temp(y,a)=abundance_spawn(y,a)*wt_mat_mult(a);
       ssb(y)=sum(ssb_temp(y));
       SPR_yr(y)=ssb(y)/SSB_zero;      
      }
     }
      if(y>1) 
      {
       if(a==1) //age-0 recruits
       {
         if(ASS_REC_FXN_type==0)
          {
           BH_REC(y-1)=R_ave;
           recruits(y)=BH_REC(y-1)*rec_devs(y-1);
          }
         if(ASS_REC_FXN_type==1)
          {
           BH_REC(y-1)=((4*h*R_ave*ssb(y-1))/(SSB_zero*(1-h)+ssb(y-1)*(5*h-1)));
           recruits(y)=BH_REC(y-1)*rec_devs(y-1);
          }  
          abundance_at_age(y,a)=recruits(y); 
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn); 
        }
        if (a==2)
        {
          abundance_at_age(y,a)=abundance_at_age(y-1,a-1)*exp(-(F_at_age(y-1,a-1)+M(a-1))*(1-tspawn));
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn);
        }        
        if (a>2 && a<nages)
        {
          abundance_at_age(y,a)=abundance_at_age(y-1,a-1)*exp(-(F_at_age(y-1,a-1)+M(a-1)));
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn);
        }
        if (a==nages)
        {
          abundance_at_age(y,a)=abundance_at_age(y-1,a-1)*exp(-(F_at_age(y-1,a-1)+M(a-1)))+abundance_at_age(y-1,a)*exp(-(F_at_age(y-1,a)+M(a)));
          abundance_spawn(y,a)=abundance_at_age(y,a)*exp(-(F_at_age(y,a)+M(a))*tspawn);
         }
       ssb_temp(y,a)=abundance_spawn(y,a)*wt_mat_mult(a);
       ssb(y)=sum(ssb_temp(y));
       SPR_yr(y)=ssb(y)/SSB_zero;       
      }
     }
    }
FUNCTION get_biomass
  for (int y=1;y<=nyrs;y++)
    {
      for (int a=1;a<=nages;a++)
      {
        biomass_temp(y,a)=wt_at_age(a)*abundance_at_age(y,a);
        biomass(y)=sum(biomass_temp(y));
      }
     }
FUNCTION get_catch
  for (int y=1;y<=nyrs_landings;y++)
    {
      for (int a=1;a<=nages;a++)
      {
       if(a==1)
        {
          //baranov's catch equation
           pred_catch_at_age(y,a)=abundance_at_age(y,a)*(1.-exp(-(F_at_age(y,a)+M(a))*(1-tspawn)))*(F_at_age(y,a))/(F_at_age(y,a)+M(a));
        }
       if(a>1)
        {
          //baranov's catch equation
           pred_catch_at_age(y,a)=abundance_at_age(y,a)*(1.-exp(-(F_at_age(y,a)+M(a))))*(F_at_age(y,a))/(F_at_age(y,a)+M(a));
        }
       }
        catchsum(y)=sum(pred_catch_at_age(y));
      for (int a=1;a<=nages;a++)
      {        
        pred_catch_prop(y,a)=pred_catch_at_age(y,a)/catchsum(y);
      }
    }
FUNCTION get_landings
  for (int y=1;y<=nyrs_landings;y++)
    {
      for (int a=1;a<=nages;a++)
      {
        pred_landings_temp(y,a)=wt_at_age(a)*pred_catch_at_age(y,a);
        pred_landings(y)=sum(pred_landings_temp(y));
      }
    }
FUNCTION get_fish_ind_indices  //assume survey occurs Jan 1
  q=mfexp(ln_q)/(mfexp(ln_q)+1);
  t_survey=nyrs-nyrs_survey;
      for (int y=1;y<=nyrs_survey;y++)
        {
          for (int a=1;a<=nages;a++)
           {
             Index_age(y,a)=selectivity_index(y,a)*abundance_at_age(y+t_survey,a);
           }
         pred_Index(y)=q*sum(Index_age(y));
         Indexsum(y)=sum(Index_age(y));
          for (int a=1;a<=nages;a++)
           {
             pred_Index_prop(y,a)=Index_age(y,a)/Indexsum(y);
           }
         }

FUNCTION evaluate_the_objective_function
   f=0.0;

  for (int w=1;w<=nyrs-1;w++)
  {
    res_recruit(w)=log(ASS_sigma_recruit)+(.5*(square(log(recruits(w+1))-log(BH_REC(w))))/(square(ASS_sigma_recruit)));
  }
  rss_recruit=sum(res_recruit);
  weighted_recruit=wt_recruit*rss_recruit;
  f+=weighted_recruit;

 if(Fit_OBS==0)
 {
  res_index=log(OBS_Index+MNconst)-log(pred_Index+MNconst);
  for (int y=1;y<=nyrs_survey;y++)
  {  
    rss_index_temp(y)=log(sigma_index(y))+(.5*(square(log(OBS_Index(y)+MNconst)-log(pred_Index(y)+MNconst)))/(square(sigma_index(y))));
  }
  rss_index=sum(rss_index_temp);
  weighted_index=wt_index*rss_index;
  f+=weighted_index;

  res_index_prop=log(OBS_index_prop+MNconst)-log(pred_Index_prop+MNconst);
  rss_index_prop_temp=rowsum(elem_prod(OBS_index_prop,res_index_prop));
  rss_index_prop=sum(nindex*rss_index_prop_temp);

  weighted_index_prop=wt_index_prop*rss_index_prop;
  f+=weighted_index_prop;
  
  res_catch_prop=log(OBS_catch_prop+MNconst)-log(pred_catch_prop+MNconst);
  rss_catch_prop_temp=rowsum(elem_prod(OBS_catch_prop,res_catch_prop));
  rss_catch_prop=sum(ncatch*rss_catch_prop_temp);
  
  res_landings=log(OBS_landings+MNconst)-log(pred_landings+MNconst);
  for (int y=1;y<=nyrs_landings;y++)
  {    
    rss_landings_temp(y)=log(sigma_landings(y))+(.5*(square(log(OBS_landings(y)+MNconst)-log(pred_landings(y)+MNconst)))/(square(sigma_landings(y))));
  }
  rss_landings=sum(rss_landings_temp);

  weighted_landings=wt_landings*rss_landings;
  weighted_catch_prop=wt_catch_prop*rss_catch_prop;

  f+=weighted_catch_prop+weighted_landings; 
 }
 if(Fit_OBS==1)  //self-consistency runs
 {
  res_index=log(TRUE_Index+MNconst)-log(pred_Index+MNconst);
  for (int y=1;y<=nyrs_survey;y++)
  {  
    rss_index_temp(y)=log(sigma_index(y))+(.5*(square(log(TRUE_Index(y)+MNconst)-log(pred_Index(y)+MNconst)))/(square(sigma_index(y))));
  }
  rss_index=sum(rss_index_temp);
  weighted_index=wt_index*rss_index;
  f+=weighted_index;

  res_index_prop=log(TRUE_Index_prop+MNconst)-log(pred_Index_prop+MNconst);
  rss_index_prop_temp=rowsum(elem_prod(TRUE_Index_prop,res_index_prop));
  rss_index_prop=sum(nindex*rss_index_prop_temp);

  weighted_index_prop=wt_index_prop*rss_index_prop;
  f+=weighted_index_prop;
  
  res_catch_prop=log(TRUE_catch_prop+MNconst)-log(pred_catch_prop+MNconst);
  rss_catch_prop_temp=rowsum(elem_prod(TRUE_catch_prop,res_catch_prop));
  rss_catch_prop=sum(ncatch*rss_catch_prop_temp);
  
  res_landings=log(TRUE_landings+MNconst)-log(pred_landings+MNconst);
  for (int y=1;y<=nyrs_landings;y++)
  {    
    rss_landings_temp(y)=log(sigma_landings(y))+(.5*(square(log(TRUE_landings(y)+MNconst)-log(pred_landings(y)+MNconst)))/(square(sigma_landings(y))));
  }
  rss_landings=sum(rss_landings_temp);

  weighted_landings=wt_landings*rss_landings;
  weighted_catch_prop=wt_catch_prop*rss_catch_prop;

  f+=weighted_catch_prop+weighted_landings; 
 }
  if(F_devs_pen_switch==1)
  {
   if(active(ln_F_devs))
   {
    F_devs_pen=F_devs_pen_mult*norm2(ln_F_devs);
    f+=F_devs_pen;
   }
  }
  if(R_devs_pen_switch==1)
  {
   if(active(ln_rec_dev))
   {
    R_devs_pen=R_devs_pen_mult*norm2(ln_rec_dev);
    f+=R_devs_pen;
   }
  }  
REPORT_SECTION
  report<<"$max_gradient"<<endl;
  report<<objective_function_value::gmax<<endl;
  report<<"$Obj_Func"<<endl;
  report<<f<<endl;
  report<<"$SIM_REC_FXN_type"<<endl;
  report<<SIM_REC_FXN_type<<endl;
  report<<"$ASS_REC_FXN_type"<<endl;
  report<<ASS_REC_FXN_type<<endl;    
  report<<"$est_steep_switch"<<endl;
  report<<est_steep_switch<<endl;
  report<<"$init_abund_switch"<<endl;
  report<<init_abund_switch<<endl;
  report<<"$SIM_survey_sel_switch"<<endl;
  report<<SIM_survey_sel_switch<<endl;
  report<<"$survey_sel_switch"<<endl;
  report<<survey_sel_switch<<endl;
  report<<"$SIM_sel_switch"<<endl;
  report<<SIM_sel_switch<<endl;
  report<<"$sel_switch"<<endl;
  report<<sel_switch<<endl;
  report<<"$TAC_CNST_switch"<<endl;
  report<<TAC_CNST_switch<<endl;
  report<<"$TAC_timelag"<<endl;
  report<<TAC_timelag<<endl;
  report<<"$TAC_type"<<endl;
  report<<TAC_type<<endl;
  report<<"$Fit_OBS"<<endl;
  report<<Fit_OBS<<endl; 
  
  report<<"$sigma_recruit"<<endl;
  report<<sigma_recruit<<endl;
  report<<"$ASS_sigma_recruit"<<endl;
  report<<ASS_sigma_recruit<<endl;
  report<<"$sigma_F"<<endl;
  report<<sigma_F<<endl;
  report<<"$ncatch"<<endl;
  report<<ncatch<<endl;
  report<<"$nindex"<<endl;
  report<<nindex<<endl;
  report<<"$sigma_index"<<endl;
  report<<sigma_index<<endl;
  report<<"$SIM_sigma_index"<<endl;
  report<<SIM_sigma_index<<endl;  
  report<<"$sigma_landings"<<endl;
  report<<sigma_landings<<endl;

  report<<"$SIM_h"<<endl;
  report<<SIM_h<<endl;
  report<<"$h"<<endl;
  report<<h<<endl;
  report<<"$SIM_R_ave"<<endl;
  report<<SIM_R_ave<<endl;
  report<<"$R_ave"<<endl;
  report<<R_ave<<endl;
  report<<"$R0_offset"<<endl;
  report<<R0_offset<<endl;
  report<<"$SIM_SSB_zero"<<endl;
  report<<SIM_SSB_zero<<endl;
  report<<"$SSB_zero"<<endl;
  report<<SSB_zero<<endl;

  report<<"$SIM_F"<<endl;
  report<<SIM_F<<endl;
  report<<"$fmult"<<endl;
  report<<fmult<<endl;
  report<<"$SIM_ssb"<<endl;
  report<<SIM_ssb<<endl;
  report<<"$ssb"<<endl;
  report<<ssb<<endl;
  report<<"$SIM_SPR_yr"<<endl;
  report<<SIM_SPR_yr<<endl;
  report<<"$SPR_yr"<<endl;
  report<<SPR_yr<<endl;    
  report<<"$SIM_BH_REC"<<endl;
  report<<SIM_BH_REC<<endl;
  report<<"$BH_REC"<<endl;
  report<<BH_REC<<endl;
  report<<"$SIM_recruits"<<endl;
  report<<SIM_recruits<<endl;
  report<<"$recruits"<<endl;
  report<<recruits<<endl;
  report<<"$SIM_biomass"<<endl;
  report<<SIM_biomass<<endl;
  report<<"$biomass"<<endl;
  report<<biomass<<endl;
  
  report<<"$TAC_age"<<endl;
  report<<SIM_TAC_age<<endl;
  report<<"$TAC"<<endl;
  report<<SIM_TAC<<endl;
  report<<"$SIM_F_at_age"<<endl;
  report<<SIM_F_at_age<<endl;
  report<<"$F_at_age"<<endl;
  report<<F_at_age<<endl;
  report<<"$SIM_selectivity"<<endl;
  report<<SIM_selectivity<<endl;
  report<<"$selectivity"<<endl;
  report<<selectivity<<endl;
  report<<"$SIM_selectivity_index"<<endl;
  report<<SIM_selectivity_index<<endl;
  report<<"$selectivity_index"<<endl;
  report<<selectivity_index<<endl;
  report<<"$Fmsy"<<endl;
  report<<Fmsy<<endl;
  report<<"$Fmsy_scalar"<<endl;
  report<<Fmsy_scalar<<endl;
  report<<"$TAC_scalar"<<endl;
  report<<TAC_scalar<<endl;
  report<<"$SIM_beta1_index"<<endl;
  report<<SIM_beta1_index<<endl;
  report<<"$beta1_index"<<endl;
  report<<beta1_index<<endl;
  report<<"$SIM_beta2_index"<<endl;
  report<<SIM_beta2_index<<endl;
  report<<"$beta2_index"<<endl;
  report<<beta2_index<<endl;
  report<<"$SIM_beta1"<<endl;
  report<<SIM_beta1<<endl;
  report<<"$beta1"<<endl;
  report<<beta1<<endl;
  report<<"$SIM_beta2"<<endl;
  report<<SIM_beta2<<endl;
  report<<"$beta2"<<endl;
  report<<beta2<<endl;
  report<<"$SIM_q"<<endl;
  report<<SIM_q<<endl;
  report<<"$q"<<endl;
  report<<q<<endl;

  report<<"$SIM_abundance_at_age"<<endl;
  report<<SIM_abundance_at_age<<endl;
  report<<"$abundance_at_age"<<endl;
  report<<abundance_at_age<<endl;
  report<<"$SIM_abundance_spawn"<<endl;
  report<<SIM_abundance_spawn<<endl;
  report<<"$abundance_spawn"<<endl;
  report<<abundance_spawn<<endl;
  report<<"$SIM_init_abund"<<endl;
  report<<SIM_init_abund<<endl;
  report<<"$init_abund"<<endl;
  report<<init_abund<<endl;

  report<<"$SIM_rec_devs"<<endl;
  report<<SIM_rec_devs<<endl;
  report<<"$rec_devs"<<endl;
  report<<rec_devs<<endl;
  
  report<<"$TRUE_catch_prop"<<endl;
  report<<TRUE_catch_prop<<endl;
  report<<"$OBS_catch_prop"<<endl;
  report<<OBS_catch_prop<<endl;
  report<<"$pred_catch_prop"<<endl;
  report<<pred_catch_prop<<endl;
  report<<"$TRUE_landings"<<endl;
  report<<TRUE_landings<<endl;
  report<<"$OBS_landings"<<endl;
  report<<OBS_landings<<endl;
  report<<"$pred_landings"<<endl;
  report<<pred_landings<<endl;
  report<<"$TRUE_Index"<<endl;
  report<<TRUE_Index<<endl;
  report<<"$OBS_Index"<<endl;
  report<<OBS_Index<<endl;
  report<<"$pred_Index"<<endl;
  report<<pred_Index<<endl;
  report<<"$TRUE_Index_prop"<<endl;
  report<<TRUE_Index_prop<<endl;
  report<<"$OBS_index_prop"<<endl;
  report<<OBS_index_prop<<endl;
  report<<"$pred_Index_prop"<<endl;
  report<<pred_Index_prop<<endl;

  report<<"$pred_catch_at_age"<<endl;
  report<<pred_catch_at_age<<endl;

  report<<"$res_recruit"<<endl;
  report<<res_recruit<<endl;
  report<<"$rss_recruit"<<endl;
  report<<rss_recruit<<endl;
  report<<"$res_landings"<<endl;
  report<<res_landings<<endl;
  report<<"$rss_landings"<<endl;
  report<<rss_landings<<endl;
  report<<"$res_catch_prop"<<endl;
  report<<res_catch_prop<<endl;
  report<<"$rss_catch_prop"<<endl;
  report<<rss_catch_prop<<endl;
  report<<"$res_index"<<endl;
  report<<res_index<<endl;
  report<<"$rss_index"<<endl;
  report<<rss_index<<endl;
  report<<"$res_index_prop"<<endl;
  report<<res_index_prop<<endl;
  report<<"$rss_index_prop"<<endl;
  report<<rss_index_prop<<endl;
  report<<"$F_devs_pen_switch"<<endl;
  report<<F_devs_pen_switch<<endl;
  report<<"$F_devs_pen"<<endl;
  report<<F_devs_pen<<endl;
  report<<"$F_devs_pen_mult"<<endl;
  report<<F_devs_pen_mult<<endl;
  report<<"$R_devs_pen_switch"<<endl;
  report<<R_devs_pen_switch<<endl;
  report<<"$R_devs_pen"<<endl;
  report<<R_devs_pen<<endl;  
  report<<"$R_devs_pen_mult"<<endl;
  report<<R_devs_pen_mult<<endl;
  
  report<<"$nages"<<endl;
  report<<nages<<endl;
  report<<"$nyrs_SIM"<<endl;
  report<<nyrs_SIM<<endl;
  report<<"$nyrs_assessment"<<endl;
  report<<nyrs<<endl;  
  report<<"$nyrs_landings"<<endl;
  report<<nyrs_landings<<endl;  
  report<<"$nyrs_survey"<<endl;
  report<<nyrs_survey<<endl;  
  report<<"$nyrs_TAC_SIM"<<endl;
  report<<nyrs_TAC_SIM<<endl;  

  report<<"$weight"<<endl;
  report<<wt_at_age<<endl;
  report<<"$maturity"<<endl;
  report<<maturity<<endl;
  report<<"$sex_scalar"<<endl;
  report<<sex_scalar<<endl;
  report<<"$M"<<endl;
  report<<M<<endl;
  report<<"$tspawn"<<endl;
  report<<tspawn<<endl;
  
  report<<"$phase_dummy"<<endl;
  report<<phase_dummy<<endl;
  report<<"$phase_beta1"<<endl;
  report<<phase_beta1<<endl;
  report<<"$phase_beta2"<<endl;
  report<<phase_beta2<<endl;
  report<<"$phase_beta1_index"<<endl;
  report<<phase_beta1_index<<endl;
  report<<"$phase_beta2_index"<<endl;
  report<<phase_beta2_index<<endl;
  report<<"$phase_ln_sel_age"<<endl;
  report<<phase_ln_sel_age<<endl;
  report<<"$phase_ln_sel_age_index"<<endl;
  report<<phase_ln_sel_age_index<<endl;
  report<<"$phase_F"<<endl;
  report<<phase_F<<endl;
  report<<"$phase_F_devs"<<endl;
  report<<phase_F_devs<<endl;
  report<<"$phase_rec_dev"<<endl;
  report<<phase_rec_dev<<endl;
  report<<"$phase_ln_R_ave"<<endl;
  report<<phase_ln_R_ave<<endl;
  report<<"$phase_h"<<endl;
  report<<phase_h<<endl;
  report<<"$phase_ln_R0_offset"<<endl;
  report<<phase_ln_R0_offset<<endl;
  report<<"$phase_ln_init_abund"<<endl;
  report<<phase_ln_init_abund<<endl;
  report<<"$phase_q"<<endl;
  report<<phase_q<<endl;

  report<<"$wt_index"<<endl;
  report<<wt_index<<endl;
  report<<"$wt_recruit"<<endl;
  report<<wt_recruit<<endl;
  report<<"$wt_index_prop"<<endl;
  report<<wt_index_prop<<endl;
  report<<"$wt_landings"<<endl;
  report<<wt_landings<<endl;
  report<<"$wt_catch_prop"<<endl;
  report<<wt_catch_prop<<endl;
RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4, 1.0e-7
  maximum_function_evaluations 100000
  


