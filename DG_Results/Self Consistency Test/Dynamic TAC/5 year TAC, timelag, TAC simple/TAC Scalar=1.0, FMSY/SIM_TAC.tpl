

GLOBALS_SECTION
  #include "admodel.h"
  #define EOUT(var) cout <<#var<<" "<<var<<endl
TOP_OF_MAIN_SECTION
  arrmblsize=500000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000000);
DATA_SECTION
  init_int nages
  init_int nyrs
  init_int nfleets //not yet implemented for more than 1 fleet
  init_int yr_fishing_start
  init_int yr_TAC_start
  init_int yr_survey_start
  
  int nyrs_landings
  int nyrs_survey
  int nyrs_TAC
  int nyrs_F_pre_TAC

  int t_F
  int t_F_post
  int t_TAC
  int t_survey
  
!! nyrs_landings=nyrs-yr_fishing_start+1;
!! nyrs_survey=nyrs-yr_survey_start+1;
!! nyrs_TAC=nyrs-yr_TAC_start+1;
!! nyrs_F_pre_TAC=nyrs_landings-nyrs_TAC+1;

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
  init_number F_initial
  init_number F_devs_initial
  init_number rec_dev_initial
  init_number ln_R_ave_initial
  init_number h_initial
  init_number ln_R0_offset_initial
  init_number ln_init_abund_initial
  init_number q_initial

  init_int phase_dummy
  init_int phase_dummy_assessment
  init_number dummy_initial

  init_int SIM_REC_FXN_type //==0 Ave REC, ==1 BH
  init_int ASS_REC_FXN_type  //==0 Ave REC, ==1 BH
  init_int est_steep_switch //==0 fix steep at SIM_h, ==1 estimate
  init_int init_abund_switch // how to estimate initial abundance for the assessment model, ==0 assume equilibrium age strcuture from an offset from R0, ==1 estimate initial abundance at each age
  init_int SIM_survey_sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int survey_sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int SIM_sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int sel_switch //survey selectivity switch, ==1 input, ==2 logistic
  init_int TAC_CNST_switch //TAC switch ==1 TAC is set constant to level in yr_TAC_start, ==0 TAC is recalculated each year
  init_int TAC_timelag  // ==0 there is no timelag and TAC is based on abundance at start of year, otehrwise TAC includes a 1 year timelage and is based on abundance at start of previous year
  init_int TAC_type  //==0 simple TAC=FMSY*exploitable Biomass, ==1 TAC=Baranov' catch equation fishing at FMSY
  init_int Fit_OBS  //==0 yes fit observed data, ==1 no fit true values (i.e., perform self-consistency run)
  init_int SIM_Fmsy_switch //==1 do an Fmsy simulation with no TAC, ==0 use Fmsy to simulate population prior to TAC
  init_int SIM_deterministic_recruit //==1 perform deterministic recruitment simulations (rec_devs=1), ==0 perform stochastic recruit simulations
  init_int F_devs_pen_switch //should include penalty on F_devs to avoid extreme values?  ==1 TRUE, ==0 FALSE
  init_int R_devs_pen_switch //should include penalty on R_devs to avoid extreme values?  ==1 TRUE, ==0 FALSE

  init_number F_devs_pen_mult
  init_number R_devs_pen_mult
  init_number SIM_h
  init_number SIM_R_ave
  init_number sigma_recruit
  init_number ASS_sigma_recruit
  init_number sigma_F
  init_number Fmsy
  init_number Fmsy_scalar
  init_number TAC_scalar
  init_number SIM_beta1_index
  init_number SIM_beta2_index
  init_vector SIM_input_sel_index(1,nages)
  init_number SIM_beta1
  init_number SIM_beta2
  init_vector SIM_input_sel(1,nages)
  init_number SIM_q_survey
  init_vector wt_at_age(1,nages)
  init_vector maturity(1,nages)
  init_number sex_scalar
  init_vector M(1,nages)
  init_number tspawn

  init_number Fnew_start
  init_number NR_dev
  init_number NR_iterations
  init_number max_Fnew
  
  init_number SIM_ncatch
  init_number SIM_nindex
  init_vector sigma_index(1,nyrs_survey)
  init_vector ASS_sigma_index(1,nyrs_survey)
  init_vector sigma_landings(1,nyrs_landings)

  init_number wt_index
  init_number wt_recruit
  init_number wt_index_prop
  init_number wt_landings
  init_number wt_catch_prop
  init_number MNconst

  init_int debug
  init_int myseed
  
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

  vector rand_rec(1,nyrs-1)
  vector rand_F(1,nyrs_F_pre_TAC)
  
 !! cout << "input read" << endl;
 !! cout << "If Debug=-999 then input read correctly" << endl;
 !! cout << "Debug=" << endl;
 !! cout << debug << endl;


PARAMETER_SECTION
  !! cout << "reading parameters"<<endl;

  init_number dummy(phase_dummy)
  
  matrix SIM_selectivity(1,nyrs_landings,1,nages)
  matrix SIM_selectivity_index(1,nyrs_survey,1,nages)
  vector SIM_rec_devs(1,nyrs-1)  //don't have deviation in first year because assume at equilibrium
  vector SIM_F_devs(1,nyrs_F_pre_TAC)
  vector SIM_F_pre_TAC(1,nyrs_F_pre_TAC)
  vector wt_mat_mult(1,nages)
  vector SPR_N(1,nages)
  vector SPR_SSB(1,nages)
  number SPR
  number SIM_SSB_zero

  vector SIM_recruits(1,nyrs) 
  matrix SIM_abundance_at_age(1,nyrs,1,nages) //before movement and mortality
  matrix SIM_abundance_spawn(1,nyrs,1,nages)
  vector SIM_init_abund(1,nages)
  vector SIM_ssb(1,nyrs)
  vector SIM_SPR_yr(1,nyrs)
  matrix ssb_temp(1,nyrs,1,nages)
  vector SIM_BH_REC(1,nyrs-1)
  vector SIM_F(1,nyrs_landings)
  matrix SIM_F_at_age(1,nyrs_landings,1,nages)
  matrix TAC_age(1,nyrs_TAC+1,1,nages)
  vector TAC(1,nyrs_TAC+1)
  
  number Fnew
  number delt
  vector fofFvect(1,nages)
  vector fprimeFhigh(1,nages)
  vector fprimeFlow(1,nages)
  number fofF
  number fprimeF

  vector SIM_F_post_TAC(1,nyrs_TAC)
  matrix biomass_temp(1,nyrs,1,nages)
  vector SIM_biomass(1,nyrs)
  matrix TRUE_catch_at_age(1,nyrs_landings,1,nages)
  vector catchsum(1,nyrs_landings)
  matrix TRUE_catch_prop(1,nyrs_landings,1,nages)
  vector TRUE_landings(1,nyrs_landings)
  matrix landings_temp(1,nyrs_landings,1,nages)
  
  vector rand_SIM_catch_prop_temp(1,SIM_ncatch)
  vector rand_SIM_catch_prop_temp2(1,nages)
  matrix SIM_catch_prop(1,nyrs_landings,1,nages)
  matrix OBS_catch_prop(1,nyrs_landings,1,nages)
  vector OBS_landings(1,nyrs_landings)

  matrix TRUE_Index_age(1,nyrs_survey,1,nages)
  vector TRUE_Index(1,nyrs_survey)
  vector Indexsum(1,nyrs_survey)
  matrix TRUE_Index_prop(1,nyrs_survey,1,nages)
  vector OBS_Index(1,nyrs_survey)
  vector rand_SIM_index_prop_temp(1,SIM_nindex)
  vector rand_SIM_index_prop_temp2(1,nages)
  matrix SIM_index_prop(1,nyrs_survey,1,nages)
  matrix OBS_index_prop(1,nyrs_survey,1,nages)

  objective_function_value f
  !! cout << "input read" << endl;

INITIALIZATION_SECTION  //set initial values
  dummy dummy_initial;
PROCEDURE_SECTION

  get_survey_selectivity();

  get_selectivity();

  get_abundance();

  get_biomass();

  get_catch();

  get_rand_CAA();

  get_landings();

  get_rand_landings();

  //get_fish_dep_indices();

  //get_rand_fish_dep_indices();

  get_fish_ind_indices();

  get_rand_fish_ind_indices();
  
  evaluate_the_objective_function();

FUNCTION get_survey_selectivity
 for (int a=1;a<=nages;a++)
  {
    for (int y=1;y<=nyrs_survey;y++)
        {
      if(SIM_survey_sel_switch==1)
       {
        SIM_selectivity_index(y,a)=SIM_input_sel_index(a);
        }
      if(SIM_survey_sel_switch==2)
       {    
        SIM_selectivity_index(y,a)=1/(1+mfexp(-log(19)*(a-(SIM_beta1_index))/(SIM_beta2_index)));
       }   
     }
   }

FUNCTION get_selectivity

  for (int a=1;a<=nages;a++)
   {
    for (int y=1;y<=nyrs_landings;y++)
      {
               if(SIM_sel_switch==1)
                 {
                  SIM_selectivity(y,a)=SIM_input_sel(a);                 
                 }
               if(SIM_sel_switch==2)
                 {
                  SIM_selectivity(y,a)=1/(1+mfexp(-log(19)*(a-(SIM_beta1))/(SIM_beta2)));                 
                 }
        }
     }
FUNCTION get_abundance
   t_F=nyrs-nyrs_landings;  //diff between total years and number of years of landings+1 to match start year for F adjustment where prior to TAC
   t_F_post=nyrs-nyrs_TAC;  //F time adjustment for years where TAC used
   t_TAC=t_F_post-1; //TAc time adjustment for years where TAC calculated (TAC is calculated year prior to implementation so contains an extra year compared to F_post_TAC

  random_number_generator myrand(15334); //make seed constant because want the same set of devs for all simulations, unlike seeds for OBS data, which need to change each run
  random_number_generator myrand1(35484);
  rand_rec.fill_randn(myrand);
  rand_F.fill_randn(myrand1);

  for (int n=1;n<=nyrs-1;n++)
   {
    if(SIM_deterministic_recruit==1) //if performing deterministic simulations then turn off recruit devs (ie set to 1.0)
    {   
     SIM_rec_devs(n)=1;
    }
    if(SIM_deterministic_recruit==0) //if performing stochastic simulations then turn on recruit devs
    {
     SIM_rec_devs(n)=mfexp(rand_rec(n)*sigma_recruit-0.5*square(sigma_recruit));
    }  
   }

  for (int n=1;n<=nyrs_F_pre_TAC;n++) //prior to TAC fishery is assumed to fluctuate around scalar*FMSY if scalar=1.0 then fluctuates around FMSY otherwise can fluctuate around some hi or low F value to increase/decrease stock size relative to BMSY
   {
    SIM_F_devs(n)=mfexp(rand_F(n)*sigma_F-0.5*square(sigma_F));
    if(SIM_Fmsy_switch==1) //set F=input F for all years
    {
     SIM_F_pre_TAC(n)=Fmsy;
    }
    if(SIM_Fmsy_switch==0) //set F=scalar*Fmsy
    {    
     SIM_F_pre_TAC(n)=Fmsy_scalar*Fmsy*SIM_F_devs(n);
    }
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
  SIM_SSB_zero=SIM_R_ave*SPR;

 for (int y=1;y<=nyrs;y++)
  {
    for (int a=1;a<=nages;a++)
    {
      if (y==1)
      {
    //equilibrium abundance in first year
        if(a==1) //age-0 recruits
        {
          SIM_recruits(y)=SIM_R_ave; 
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 

          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn); 
          SIM_init_abund(a)=SIM_abundance_at_age(y,a);
        }
        if(a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y,a-1)*mfexp(-M(a-1)*(1-tspawn));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn); 
          SIM_init_abund(a)=SIM_abundance_at_age(y,a);
        }        
        if(a>2 && a!=nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y,a-1)*mfexp(-M(a-1));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn); 
          SIM_init_abund(a)=SIM_abundance_at_age(y,a);
        }
        if(a==nages)
        {
         SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y,a-1)*mfexp(-M(a-1))/(1-mfexp(-M(a)));
         SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn); 
         SIM_init_abund(a)=SIM_abundance_at_age(y,a);
        }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;
      }

      if(y>1 && y<yr_fishing_start) 
      {
       if(a==1) //age-0 recruits
       {
         if(SIM_REC_FXN_type==0)
          {
           SIM_BH_REC(y-1)=SIM_R_ave;
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }
         if(SIM_REC_FXN_type==1)
          {
           SIM_BH_REC(y-1)=((4*SIM_h*SIM_R_ave*SIM_ssb(y-1))/(SIM_SSB_zero*(1-SIM_h)+SIM_ssb(y-1)*(5*SIM_h-1)));
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }          
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn); 
        }
        if (a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-M(a-1)*(1-tspawn));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn);
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-M(a-1));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn);
        }
        if (a==nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-M(a-1))+SIM_abundance_at_age(y-1,a)*exp(-M(a));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-M(a)*tspawn);
         }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;
       }	
      if(y==yr_fishing_start) 
      {
       if(a==1) //age-0 recruits
       {
         if(SIM_REC_FXN_type==0)
          {
           SIM_BH_REC(y-1)=SIM_R_ave;
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }
         if(SIM_REC_FXN_type==1)
          {
           SIM_BH_REC(y-1)=((4*SIM_h*SIM_R_ave*SIM_ssb(y-1))/(SIM_SSB_zero*(1-SIM_h)+SIM_ssb(y-1)*(5*SIM_h-1)));
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }    
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn); 
        }
        if (a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-M(a-1)*(1-tspawn));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-M(a-1));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }
        if (a==nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-M(a-1))+SIM_abundance_at_age(y-1,a)*exp(-M(a));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
         }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;       
       SIM_F(y-t_F)=SIM_F_pre_TAC(y-t_F);
      }

      if(y>yr_fishing_start && y<(yr_TAC_start-1)) 
      {
       if(a==1) //age-0 recruits
       {
         if(SIM_REC_FXN_type==0)
          {
           SIM_BH_REC(y-1)=SIM_R_ave;
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }
         if(SIM_REC_FXN_type==1)
          {
           SIM_BH_REC(y-1)=((4*SIM_h*SIM_R_ave*SIM_ssb(y-1))/(SIM_SSB_zero*(1-SIM_h)+SIM_ssb(y-1)*(5*SIM_h-1)));
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }    
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn); 
        }
        if (a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1))*(1-tspawn));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }
        if (a==nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)))+SIM_abundance_at_age(y-1,a)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a)+M(a)));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
         }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;       
       SIM_F(y-t_F)=SIM_F_pre_TAC(y-t_F);
      }	

      if(y==(yr_TAC_start-1)) //determine TAC in year prior to implementation (1 year lag in implementation)
      {
       if(a==1) //age-0 recruits
       {
         if(SIM_REC_FXN_type==0)
          {
           SIM_BH_REC(y-1)=SIM_R_ave;
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }
         if(SIM_REC_FXN_type==1)
          {
           SIM_BH_REC(y-1)=((4*SIM_h*SIM_R_ave*SIM_ssb(y-1))/(SIM_SSB_zero*(1-SIM_h)+SIM_ssb(y-1)*(5*SIM_h-1)));
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }    
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn); 
          TAC_age(y-t_TAC,a)=Fmsy*(1-tspawn)*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a); //TAC assumes FMSY HCR...ie FMSY*ExploitableBiomass
          if(TAC_type==1)
           {
            TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))*(1-tspawn)))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
           }
        }
        if (a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1))*(1-tspawn));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
          TAC_age(y-t_TAC,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a); //TAC assumes FMSY HCR...ie FMSY*ExploitableBiomass
          if(TAC_type==1)
           {
            TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
           }
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
          TAC_age(y-t_TAC,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a); //TAC assumes FMSY HCR...ie FMSY*ExploitableBiomass
          if(TAC_type==1)
           {
            TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
           }
        }
        if (a==nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)))+SIM_abundance_at_age(y-1,a)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a)+M(a)));
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_pre_TAC(y-t_F)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
          TAC_age(y-t_TAC,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a); //TAC assumes FMSY HCR...ie FMSY*ExploitableBiomass
          if(TAC_type==1)
           {
            TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
           }
         }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;       
       SIM_F(y-t_F)=SIM_F_pre_TAC(y-t_F);
       TAC(y-t_TAC)=TAC_scalar*sum(TAC_age(y-t_TAC));
      }	
    }

      if(y==yr_TAC_start) 
      {
       for (int a=1;a<=nages;a++)
        {
       if(a==1) //age-0 recruits
       {
         if(SIM_REC_FXN_type==0)
          {
           SIM_BH_REC(y-1)=SIM_R_ave;
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }
         if(SIM_REC_FXN_type==1)
          {
           SIM_BH_REC(y-1)=((4*SIM_h*SIM_R_ave*SIM_ssb(y-1))/(SIM_SSB_zero*(1-SIM_h)+SIM_ssb(y-1)*(5*SIM_h-1)));
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }    
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 
        }
        if (a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1))*(1-tspawn));
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)));
        }
        if (a==nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)))+SIM_abundance_at_age(y-1,a)*exp(-(SIM_F_pre_TAC(y-t_F-1)*SIM_selectivity(y-t_F-1,a)+M(a)));
        }
          if(TAC_timelag==0)
          {
           if(a==1)
            {
              TAC_age(y-t_TAC-1,a)=Fmsy*(1-tspawn)*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC-1,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))*(1-tspawn)))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }
           if(a>1)
            {
              TAC_age(y-t_TAC-1,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC-1,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }
           }
          }
          if(TAC_timelag==0)
           {
             TAC(y-t_TAC-1)=TAC_scalar*sum(TAC_age(y-t_TAC-1));
           }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///// Newton Raphson to fit the TAC //////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  if(SIM_Fmsy_switch==0)
   {
            if(TAC(y-t_TAC-1)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
             {
               Fnew=0;
             }
            if(TAC(y-t_TAC-1)>0)
             {
              Fnew=Fnew_start;
              for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
               {
                delt=Fnew*NR_dev;  // NR_dev~0.001
                 for(int s=1;s<=nages;s++)
                  {
                    if(s==1)
                     {
                      fofFvect(s)=wt_at_age(s)*((Fnew*SIM_selectivity(y-t_F,s))/(Fnew*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*(Fnew*SIM_selectivity(y-t_F,s)+M(s))*(1-tspawn)));
                      fprimeFhigh(s)=wt_at_age(s)*(((Fnew+delt)*SIM_selectivity(y-t_F,s))/((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s))*(1-tspawn)));
                      fprimeFlow(s)=wt_at_age(s)*(((Fnew-delt)*SIM_selectivity(y-t_F,s))/((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s))*(1-tspawn)));                                              
                     }
                    if(s>1)
                     {
                      fofFvect(s)=wt_at_age(s)*((Fnew*SIM_selectivity(y-t_F,s))/(Fnew*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*(Fnew*SIM_selectivity(y-t_F,s)+M(s))));
                      fprimeFhigh(s)=wt_at_age(s)*(((Fnew+delt)*SIM_selectivity(y-t_F,s))/((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s))));
                      fprimeFlow(s)=wt_at_age(s)*(((Fnew-delt)*SIM_selectivity(y-t_F,s))/((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s))));
                     }
                   } 
                  fofF=sum(fofFvect)-TAC(y-t_TAC-1);
                  fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                  Fnew=Fnew-(fofF/fprimeF);
                 if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                 if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                }
               } 
           SIM_F_post_TAC(y-t_F_post)=Fnew;
           SIM_F(y-t_F)=SIM_F_post_TAC(y-t_F_post);
       }
  if(SIM_Fmsy_switch==1)
   {
     SIM_F_post_TAC(y-t_F_post)=Fmsy;
     SIM_F(y-t_F)=SIM_F_post_TAC(y-t_F_post);
   }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     for (int a=1;a<=nages;a++)
      {
       if(a==1) //age-0 recruits
       {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn); 
        }
        if (a==2)
        {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }
        if (a==nages)
        {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
         }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
        if(TAC_CNST_switch==0)
          {
           if(a==1)
            {
             TAC_age(y-t_TAC,a)=Fmsy*(1-tspawn)*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))*(1-tspawn)))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }
           if(a>1)
            {
             TAC_age(y-t_TAC,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }  
          }
        }
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;       

      if(TAC_CNST_switch==1)
      {
        TAC_age(y-t_TAC)=TAC_age(y-t_TAC-1);
      }
        TAC(y-t_TAC)=TAC_scalar*sum(TAC_age(y-t_TAC));
     }

      if(y>yr_TAC_start)
      {
       for (int a=1;a<=nages;a++)
       {
       if(a==1) //age-0 recruits
       {
         if(SIM_REC_FXN_type==0)
          {
           SIM_BH_REC(y-1)=SIM_R_ave;
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }
         if(SIM_REC_FXN_type==1)
          {
           SIM_BH_REC(y-1)=((4*SIM_h*SIM_R_ave*SIM_ssb(y-1))/(SIM_SSB_zero*(1-SIM_h)+SIM_ssb(y-1)*(5*SIM_h-1)));
           SIM_recruits(y)=SIM_BH_REC(y-1)*SIM_rec_devs(y-1);
          }    
          SIM_abundance_at_age(y,a)=SIM_recruits(y); 
        }
        if (a==2)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_post_TAC(y-t_F_post-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1))*(1-tspawn));
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_post_TAC(y-t_F_post-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)));
        }
        if (a==nages)
        {
          SIM_abundance_at_age(y,a)=SIM_abundance_at_age(y-1,a-1)*exp(-(SIM_F_post_TAC(y-t_F_post-1)*SIM_selectivity(y-t_F-1,a-1)+M(a-1)))+SIM_abundance_at_age(y-1,a)*exp(-(SIM_F_post_TAC(y-t_F_post-1)*SIM_selectivity(y-t_F-1,a)+M(a)));
        }
          if(TAC_timelag==0 && TAC_CNST_switch==0)
          {
           if(a==1)
            {
              TAC_age(y-t_TAC-1,a)=Fmsy*(1-tspawn)*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC-1,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))*(1-tspawn)))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }
           if(a>1)
            {
              TAC_age(y-t_TAC-1,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC-1,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }
           }
         }
          if(TAC_timelag==0 && TAC_CNST_switch==0)
           {
             TAC(y-t_TAC-1)=TAC_scalar*sum(TAC_age(y-t_TAC-1));
           }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///// Newton Raphson to fit the TAC //////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  if(SIM_Fmsy_switch==0)
   {
            if(TAC(y-t_TAC-1)==0) //iterations have trouble finding F=0 when target=0; but they work great for >0 values.  This prevents those issues
             {
               Fnew=0;
             }
            if(TAC(y-t_TAC-1)>0)
             {
              Fnew=Fnew_start;
              for(int i=1;i<=NR_iterations;i++)  //newton-raphson iterations
               {
                delt=Fnew*NR_dev;  // NR_dev~0.001
                 for(int s=1;s<=nages;s++)
                  {
                    if(s==1)
                     {
                      fofFvect(s)=wt_at_age(s)*((Fnew*SIM_selectivity(y-t_F,s))/(Fnew*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*(Fnew*SIM_selectivity(y-t_F,s)+M(s))*(1-tspawn)));
                      fprimeFhigh(s)=wt_at_age(s)*(((Fnew+delt)*SIM_selectivity(y-t_F,s))/((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s))*(1-tspawn)));
                      fprimeFlow(s)=wt_at_age(s)*(((Fnew-delt)*SIM_selectivity(y-t_F,s))/((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s))*(1-tspawn)));                                              
                     }
                    if(s>1)
                     {
                      fofFvect(s)=wt_at_age(s)*((Fnew*SIM_selectivity(y-t_F,s))/(Fnew*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*(Fnew*SIM_selectivity(y-t_F,s)+M(s))));
                      fprimeFhigh(s)=wt_at_age(s)*(((Fnew+delt)*SIM_selectivity(y-t_F,s))/((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew+delt)*SIM_selectivity(y-t_F,s)+M(s))));
                      fprimeFlow(s)=wt_at_age(s)*(((Fnew-delt)*SIM_selectivity(y-t_F,s))/((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s)))*SIM_abundance_at_age(y,s)*(1-mfexp(-1*((Fnew-delt)*SIM_selectivity(y-t_F,s)+M(s))));
                     }
                   } 
                  fofF=sum(fofFvect)-TAC(y-t_TAC-1);
                  fprimeF=(sum(fprimeFhigh)-sum(fprimeFlow))/(2.0*delt);
                  Fnew=Fnew-(fofF/fprimeF);
                 if(Fnew<0) Fnew=0.5*(Fnew+(fofF/fprimeF));
                 if(Fnew>max_Fnew) Fnew=max_Fnew;  //At low N, Fnew would sometimes be really high and I'd get an error message.  This prevents those errors.
                }
               } 
           SIM_F_post_TAC(y-t_F_post)=Fnew;
           SIM_F(y-t_F)=SIM_F_post_TAC(y-t_F_post);
     }
  if(SIM_Fmsy_switch==1)
   {
     SIM_F_post_TAC(y-t_F_post)=Fmsy;
     SIM_F(y-t_F)=SIM_F_post_TAC(y-t_F_post);
   }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     for (int a=1;a<=nages;a++)
      {
       if(a==1) //age-0 recruits
       {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn); 
        }
        if (a==2)
        {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }        
        if (a>2 && a<nages)
        {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
        }
        if (a==nages)
        {
          SIM_abundance_spawn(y,a)=SIM_abundance_at_age(y,a)*exp(-(SIM_F_post_TAC(y-t_F_post)*SIM_selectivity(y-t_F,a)+M(a))*tspawn);
         }
       ssb_temp(y,a)=SIM_abundance_spawn(y,a)*wt_mat_mult(a);
        if(TAC_CNST_switch==0)
          {
           if(a==1)
            {
             TAC_age(y-t_TAC,a)=Fmsy*(1-tspawn)*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))*(1-tspawn)))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }
           if(a>1)
            {
             TAC_age(y-t_TAC,a)=Fmsy*SIM_selectivity(y-t_F,a)*wt_at_age(a)*SIM_abundance_at_age(y,a);
             if(TAC_type==1)
              {
                TAC_age(y-t_TAC,a)=wt_at_age(a)*SIM_abundance_at_age(y,a)*(1.-mfexp(-(Fmsy*SIM_selectivity(y-t_F,a)+M(a))))*(Fmsy*SIM_selectivity(y-t_F,a))/(Fmsy*SIM_selectivity(y-t_F,a)+M(a));
              }
            }  
          }     
      }
       SIM_ssb(y)=sum(ssb_temp(y));
       SIM_SPR_yr(y)=SIM_ssb(y)/SIM_SSB_zero;
       
      if(TAC_CNST_switch==1)
      {
        TAC_age(y-t_TAC)=TAC_age(y-t_TAC-1);
      }
        TAC(y-t_TAC)=TAC_scalar*sum(TAC_age(y-t_TAC));
     }
    }

FUNCTION get_biomass
  for (int y=1;y<=nyrs;y++)
    {
      for (int a=1;a<=nages;a++)
      {
        biomass_temp(y,a)=wt_at_age(a)*SIM_abundance_at_age(y,a);
        SIM_biomass(y)=sum(biomass_temp(y));
      }
     }
FUNCTION get_catch
  for (int y=1;y<=nyrs_landings;y++)
    {
      for (int a=1;a<=nages;a++)
      {
       if(a==1)
        {
         SIM_F_at_age(y,a)=SIM_selectivity(y,a)*SIM_F(y);
         //baranov's catch equation
         TRUE_catch_at_age(y,a)=SIM_abundance_at_age(y+t_F,a)*(1.-mfexp(-(SIM_F(y)*SIM_selectivity(y,a)+M(a))*(1-tspawn)))*(SIM_F(y)*SIM_selectivity(y,a))/(SIM_F(y)*SIM_selectivity(y,a)+M(a));
        }
       if(a>1)
        {
         SIM_F_at_age(y,a)=SIM_selectivity(y,a)*SIM_F(y);
         //baranov's catch equation
         TRUE_catch_at_age(y,a)=SIM_abundance_at_age(y+t_F,a)*(1.-mfexp(-(SIM_F(y)*SIM_selectivity(y,a)+M(a))))*(SIM_F(y)*SIM_selectivity(y,a))/(SIM_F(y)*SIM_selectivity(y,a)+M(a));
        }        
      }
        catchsum(y)=sum(TRUE_catch_at_age(y));
      for (int a=1;a<=nages;a++)
      {        
        TRUE_catch_prop(y,a)=TRUE_catch_at_age(y,a)/catchsum(y);
      }
    }
FUNCTION get_rand_CAA
  random_number_generator myrand3(myseed);
  
      for(int y=1;y<=nyrs_landings;y++)
       {
        rand_SIM_catch_prop_temp=0;
        rand_SIM_catch_prop_temp2=0;
        rand_SIM_catch_prop_temp.fill_multinomial(myrand3,value(TRUE_catch_prop(y)));
         for(int n=1;n<=SIM_ncatch;n++) /// look into changing this so can have ncatch change by year (ie different sample sizes for beginning and end of timeseries)
          {
           rand_SIM_catch_prop_temp2(value(rand_SIM_catch_prop_temp(n)))+= 1.0;
          }
        SIM_catch_prop(y)=rand_SIM_catch_prop_temp2;
       }
      for(int y=1;y<=nyrs_landings;y++)
       {
        for(int a=1;a<=nages;a++)
         {
          OBS_catch_prop(y,a)=SIM_catch_prop(y,a)/SIM_ncatch;
         }
       }

FUNCTION get_landings
  for (int y=1;y<=nyrs_landings;y++)
    {
      for (int a=1;a<=nages;a++)
      {
        landings_temp(y,a)=wt_at_age(a)*TRUE_catch_at_age(y,a);
        TRUE_landings(y)=sum(landings_temp(y));
      }
    }

FUNCTION get_rand_landings

  random_number_generator myrand2(myseed+6000);

  for (int y=1;y<=nyrs_landings;y++)
    {
     OBS_landings(y)=TRUE_landings(y)*mfexp(randn(myrand2)*sigma_landings(y)-.5*square(sigma_landings(y)));
    }
    
//FUNCTION get_fish_dep_indices  //need to check if fish dep indices are set up correctly or if they are the same as a fishery independent survey from a simulation standpoint
   //all indices in numbers not weight
   
   //q_HL_E=mfexp(ln_q_HL_E)/(mfexp(ln_q_HL_E)+1);

//      for (int y=1;y<=nyrs_landings;y++)
 //       {
 //        for (int a=1;a<=nages;a++)
//          {
//            TRUE_CPUE_age(y,a)=SIM_selectivity(y,a)*abundance_at_age(y+(t_f-1),a);
//          }
//        TRUE_CPUE(y)=q*sum(TRUE_CPUE_age(y));
 //       }


//FUNCTION get_rand_fish_dep_indices
//  random_number_generator myrand(myseed);

//  for (int y=1;y<=nyrs_landings;y++)
//    {         
//      OBS_CPUE(y)=TRUE_CPUE(y)*mfexp(randn(myrand)*sigma_CPUE(y)-.5*square(sigma_CPUE(y)));
//    }

FUNCTION get_fish_ind_indices  //assume survey occurs Jan 1
  // q_SUMM_E=mfexp(ln_q_SUMM_E)/(mfexp(ln_q_SUMM_E)+1);
 t_survey=nyrs-nyrs_survey;

      for (int y=1;y<=nyrs_survey;y++)
        {
          for (int a=1;a<=nages;a++)
           {
             TRUE_Index_age(y,a)=SIM_selectivity_index(y,a)*SIM_abundance_at_age(y+t_survey,a);
           }
         TRUE_Index(y)=SIM_q_survey*sum(TRUE_Index_age(y));
         Indexsum(y)=sum(TRUE_Index_age(y));
          for (int a=1;a<=nages;a++)
           {
             TRUE_Index_prop(y,a)=TRUE_Index_age(y,a)/Indexsum(y);
           }
         }

FUNCTION get_rand_fish_ind_indices

  random_number_generator myrand4(myseed+8300);
  random_number_generator myrand5(myseed+9400);

  for (int y=1;y<=nyrs_survey;y++)
    {         
      OBS_Index(y)=TRUE_Index(y)*mfexp(randn(myrand4)*sigma_index(y)-.5*square(sigma_index(y)));
    }
    
      for(int y=1;y<=nyrs_survey;y++)
       {
        rand_SIM_index_prop_temp=0;
        rand_SIM_index_prop_temp2=0;
        rand_SIM_index_prop_temp.fill_multinomial(myrand5,value(TRUE_Index_prop(y)));
         for(int n=1;n<=SIM_nindex;n++)
          {
           rand_SIM_index_prop_temp2(value(rand_SIM_index_prop_temp(n))) += 1.0;
          }
        SIM_index_prop(y)=rand_SIM_index_prop_temp2;
       }
      for(int y=1;y<=nyrs_survey;y++)
       {
        for(int a=1;a<=nages;a++)
         {
          OBS_index_prop(y,a)=SIM_index_prop(y,a)/SIM_nindex;
         }
       }

FUNCTION evaluate_the_objective_function
   f=0.0;
REPORT_SECTION
 if(SIM_Fmsy_switch==0)
 {
  report<<"#nages"<<endl;
  report<<nages<<endl;
  report<<"#nyrs"<<endl;
  report<<nyrs<<endl;
  report<<"#nyrs_assessment"<<endl;
  report<<nyrs_landings<<endl;  
  report<<"#nyrs_landings"<<endl;
  report<<nyrs_landings<<endl;  
  report<<"#nyrs_survey"<<endl;
  report<<nyrs_survey<<endl;  
  report<<"#nyrs_TAC"<<endl;
  report<<nyrs_TAC<<endl;  

  report<<"#phase_dummy"<<endl;
  report<<phase_dummy_assessment<<endl;
  report<<"#phase_beta1"<<endl;
  report<<phase_beta1<<endl;
  report<<"#phase_beta2"<<endl;
  report<<phase_beta2<<endl;
  report<<"#phase_beta1_index"<<endl;
  report<<phase_beta1_index<<endl;
  report<<"#phase_beta2_index"<<endl;
  report<<phase_beta2_index<<endl;
  report<<"#phase_ln_sel_age"<<endl;
  report<<phase_ln_sel_age<<endl;
  report<<"#phase_ln_sel_age_index"<<endl;
  report<<phase_ln_sel_age_index<<endl;
  report<<"#phase_F"<<endl;
  report<<phase_F<<endl;
  report<<"#phase_F_devs"<<endl;
  report<<phase_F_devs<<endl;
  report<<"#phase_rec_dev"<<endl;
  report<<phase_rec_dev<<endl;
  report<<"#phase_ln_R_ave"<<endl;
  report<<phase_ln_R_ave<<endl;
  report<<"#phase_h"<<endl;
  report<<phase_h<<endl;
  report<<"#phase_ln_R0_offset"<<endl;
  report<<phase_ln_R0_offset<<endl;
  report<<"#phase_ln_init_abund"<<endl;
  report<<phase_ln_init_abund<<endl;
  report<<"#phase_q"<<endl;
  report<<phase_q<<endl;

  report<<"#beta1_initial"<<endl;
  report<<beta1_initial<<endl;
  report<<"#beta2_initial"<<endl;
  report<<beta2_initial<<endl;
  report<<"#beta1_index_initial"<<endl;
  report<<beta1_index_initial<<endl;
  report<<"#beta2_index_initial"<<endl;
  report<<beta2_index_initial<<endl;
  report<<"#ln_sel_age_initial"<<endl;
  report<<ln_sel_age_initial<<endl;
  report<<"#ln_sel_age_index_initial"<<endl;
  report<<ln_sel_age_index_initial<<endl;
  report<<"#F_initial"<<endl;
  report<<F_initial<<endl;
  report<<"#F_devs_initial"<<endl;
  report<<F_devs_initial<<endl;
  report<<"#rec_dev_initial"<<endl;
  report<<rec_dev_initial<<endl;
  report<<"#ln_R_ave_initial"<<endl;
  report<<ln_R_ave_initial<<endl;
  report<<"#h_initial"<<endl;
  report<<h_initial<<endl;
  report<<"#ln_R0_offset_initial"<<endl;
  report<<ln_R0_offset_initial<<endl;
  report<<"#ln_init_abund_initial"<<endl;
  report<<ln_init_abund_initial<<endl;
  report<<"#q_initial"<<endl;
  report<<q_initial<<endl;

  report<<"#SIM_REC_FXN_type"<<endl;
  report<<SIM_REC_FXN_type<<endl;
  report<<"#ASS_REC_FXN_type"<<endl;
  report<<ASS_REC_FXN_type<<endl;  
  report<<"#est_steep_switch"<<endl;
  report<<est_steep_switch<<endl;
  report<<"#init_abund_switch"<<endl;
  report<<init_abund_switch<<endl;
  report<<"#SIM_survey_sel_switch"<<endl;
  report<<SIM_survey_sel_switch<<endl;
  report<<"#survey_sel_switch"<<endl;
  report<<survey_sel_switch<<endl;
  report<<"#SIM_sel_switch"<<endl;
  report<<SIM_sel_switch<<endl;
  report<<"#sel_switch"<<endl;
  report<<sel_switch<<endl;
  report<<"#TAC_CNST_switch"<<endl;
  report<<TAC_CNST_switch<<endl;
  report<<"#TAC_timelag"<<endl;
  report<<TAC_timelag<<endl;
  report<<"#TAC_type"<<endl;
  report<<TAC_type<<endl;
  report<<"#Fit_OBS"<<endl;
  report<<Fit_OBS<<endl;    
  report<<"#F_devs_pen_switch"<<endl;
  report<<F_devs_pen_switch<<endl;
  report<<"#R_devs_pen_switch"<<endl;
  report<<R_devs_pen_switch<<endl;  
  report<<"#F_devs_pen_mult"<<endl;
  report<<F_devs_pen_mult<<endl;
  report<<"#R_devs_pen_mult"<<endl;
  report<<R_devs_pen_mult<<endl;

  report<<"#weight"<<endl;
  report<<wt_at_age<<endl;
  report<<"#maturity"<<endl;
  report<<maturity<<endl;
  report<<"#sex_scalar"<<endl;
  report<<sex_scalar<<endl;
  report<<"#M"<<endl;
  report<<M<<endl;
  report<<"#tspawn"<<endl;
  report<<tspawn<<endl;

  report<<"#sigma_recruit"<<endl;
  report<<sigma_recruit<<endl;
  report<<"#ASS_sigma_recruit"<<endl;
  report<<ASS_sigma_recruit<<endl;
  report<<"#sigma_F"<<endl;
  report<<sigma_F<<endl;
  report<<"#ncatch"<<endl;
  report<<SIM_ncatch<<endl;
  report<<"#nindex"<<endl;
  report<<SIM_nindex<<endl;
  report<<"#sigma_index"<<endl;
  report<<sigma_index<<endl;
  report<<"#ASS_sigma_index"<<endl;
  report<<ASS_sigma_index<<endl;  
  report<<"#sigma_landings"<<endl;
  report<<sigma_landings<<endl;

  report<<"#wt_index"<<endl;
  report<<wt_index<<endl;
  report<<"#wt_recruit"<<endl;
  report<<wt_recruit<<endl;
  report<<"#wt_index_prop"<<endl;
  report<<wt_index_prop<<endl;
  report<<"#wt_landings"<<endl;
  report<<wt_landings<<endl;
  report<<"#wt_catch_prop"<<endl;
  report<<wt_catch_prop<<endl;
  
 ////////// Simulation Values ///////////////////////////////////////////////////////////
  report<<"#SIM_h"<<endl;
  report<<SIM_h<<endl;
  report<<"#SIM_R_ave"<<endl;
  report<<SIM_R_ave<<endl;
  report<<"#SIM_SSB_zero"<<endl;
  report<<SIM_SSB_zero<<endl;
  report<<"#Fmsy"<<endl;
  report<<Fmsy<<endl;
  report<<"#Fmsy_scalar"<<endl;
  report<<Fmsy_scalar<<endl;
  report<<"#TAC_scalar"<<endl;
  report<<TAC_scalar<<endl;
  
  report<<"#SIM_beta1_index"<<endl;
  report<<SIM_beta1_index<<endl;
  report<<"#SIM_beta2_index"<<endl;
  report<<SIM_beta2_index<<endl;
  report<<"#SIM_beta1"<<endl;
  report<<SIM_beta1<<endl;
  report<<"#SIM_beta2"<<endl;
  report<<SIM_beta2<<endl;
  report<<"#SIM_q_survey"<<endl;
  report<<SIM_q_survey<<endl;

  report<<"#SIM_F"<<endl;
  report<<SIM_F<<endl;
  report<<"#TAC_age"<<endl;
  report<<TAC_age<<endl;
  report<<"#TAC"<<endl;
  report<<TAC<<endl;
  report<<"#SIM_F_at_age"<<endl;
  report<<SIM_F_at_age<<endl;
  report<<"#SIM_selectivity"<<endl;
  report<<SIM_selectivity<<endl;
  report<<"#SIM_selectivity_index"<<endl;
  report<<SIM_selectivity_index<<endl;

  report<<"#SIM_ssb"<<endl;
  report<<SIM_ssb<<endl;
  report<<"#SIM_SPR_yr"<<endl;
  report<<SIM_SPR_yr<<endl;  
  report<<"#SIM_BH_REC"<<endl;
  report<<SIM_BH_REC<<endl;
  report<<"#SIM_recruits"<<endl;
  report<<SIM_recruits<<endl;
  report<<"#SIM_abundance_at_age"<<endl;
  report<<SIM_abundance_at_age<<endl;
  report<<"#SIM_abundance_spawn"<<endl;
  report<<SIM_abundance_spawn<<endl;
  report<<"#SIM_init_abund"<<endl;
  report<<SIM_init_abund<<endl;

  report<<"#SIM_biomass"<<endl;
  report<<SIM_biomass<<endl;
  report<<"#SIM_rec_devs"<<endl;
  report<<SIM_rec_devs<<endl;
  
  report<<"#TRUE_catch_prop"<<endl;
  report<<TRUE_catch_prop<<endl;
  report<<"#OBS_catch_prop"<<endl;
  report<<OBS_catch_prop<<endl;
  
  report<<"#TRUE_landings"<<endl;
  report<<TRUE_landings<<endl;
  report<<"#OBS_landings"<<endl;
  report<<OBS_landings<<endl;
  
  report<<"#TRUE_Index"<<endl;
  report<<TRUE_Index<<endl;
  report<<"#OBS_Index"<<endl;
  report<<OBS_Index<<endl;
  report<<"#TRUE_Index_prop"<<endl;
  report<<TRUE_Index_prop<<endl;
  report<<"#OBS_index_prop"<<endl;
  report<<OBS_index_prop<<endl;

  report<<"#MNconst"<<endl;
  report<<MNconst<<endl;
  report<<"#debug"<<endl;
  report<<debug<<endl;
 }


 if(SIM_Fmsy_switch==1)
 {
  report<<"$nages"<<endl;
  report<<nages<<endl;
  report<<"$nyrs"<<endl;
  report<<nyrs<<endl;
  report<<"$nyrs_assessment"<<endl;
  report<<nyrs_landings<<endl;  
  report<<"$nyrs_landings"<<endl;
  report<<nyrs_landings<<endl;  
  report<<"$nyrs_survey"<<endl;
  report<<nyrs_survey<<endl;  
  report<<"$nyrs_TAC"<<endl;
  report<<nyrs_TAC<<endl;  

  report<<"$phase_dummy"<<endl;
  report<<phase_dummy_assessment<<endl;
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

  report<<"$beta1_initial"<<endl;
  report<<beta1_initial<<endl;
  report<<"$beta2_initial"<<endl;
  report<<beta2_initial<<endl;
  report<<"$beta1_index_initial"<<endl;
  report<<beta1_index_initial<<endl;
  report<<"$beta2_index_initial"<<endl;
  report<<beta2_index_initial<<endl;
  report<<"$ln_sel_age_initial"<<endl;
  report<<ln_sel_age_initial<<endl;
  report<<"$ln_sel_age_index_initial"<<endl;
  report<<ln_sel_age_index_initial<<endl;
  report<<"$F_initial"<<endl;
  report<<F_initial<<endl;
  report<<"$F_devs_initial"<<endl;
  report<<F_devs_initial<<endl;
  report<<"$rec_dev_initial"<<endl;
  report<<rec_dev_initial<<endl;
  report<<"$ln_R_ave_initial"<<endl;
  report<<ln_R_ave_initial<<endl;
  report<<"$h_initial"<<endl;
  report<<h_initial<<endl;
  report<<"$ln_R0_offset_initial"<<endl;
  report<<ln_R0_offset_initial<<endl;
  report<<"$ln_init_abund_initial"<<endl;
  report<<ln_init_abund_initial<<endl;
  report<<"$q_initial"<<endl;
  report<<q_initial<<endl;

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
  report<<"$F_devs_pen_switch"<<endl;
  report<<F_devs_pen_switch<<endl;
  report<<"$R_devs_pen_switch"<<endl;
  report<<R_devs_pen_switch<<endl;  
  report<<"#F_devs_pen_mult"<<endl;
  report<<F_devs_pen_mult<<endl;
  report<<"#R_devs_pen_mult"<<endl;
  report<<R_devs_pen_mult<<endl;

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

  report<<"$sigma_recruit"<<endl;
  report<<sigma_recruit<<endl;
  report<<"$ASS_sigma_recruit"<<endl;
  report<<ASS_sigma_recruit<<endl;
  report<<"$sigma_F"<<endl;
  report<<sigma_F<<endl;
  report<<"$ncatch"<<endl;
  report<<SIM_ncatch<<endl;
  report<<"$nindex"<<endl;
  report<<SIM_nindex<<endl;
  report<<"$sigma_index"<<endl;
  report<<sigma_index<<endl;
  report<<"$ASS_sigma_index"<<endl;
  report<<ASS_sigma_index<<endl;   
  report<<"$sigma_landings"<<endl;
  report<<sigma_landings<<endl;

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
  
 ////////// Simulation Values ///////////////////////////////////////////////////////////
  report<<"$SIM_h"<<endl;
  report<<SIM_h<<endl;
  report<<"$SIM_R_ave"<<endl;
  report<<SIM_R_ave<<endl;
  report<<"$SIM_SSB_zero"<<endl;
  report<<SIM_SSB_zero<<endl;
  report<<"$Fmsy"<<endl;
  report<<Fmsy<<endl;
  report<<"$Fmsy_scalar"<<endl;
  report<<Fmsy_scalar<<endl;
  report<<"$TAC_scalar"<<endl;
  report<<TAC_scalar<<endl;
  
  report<<"$SIM_beta1_index"<<endl;
  report<<SIM_beta1_index<<endl;
  report<<"$SIM_beta2_index"<<endl;
  report<<SIM_beta2_index<<endl;
  report<<"$SIM_beta1"<<endl;
  report<<SIM_beta1<<endl;
  report<<"$SIM_beta2"<<endl;
  report<<SIM_beta2<<endl;
  report<<"$SIM_q_survey"<<endl;
  report<<SIM_q_survey<<endl;

  report<<"$SIM_F"<<endl;
  report<<SIM_F<<endl;
  report<<"$TAC_age"<<endl;
  report<<TAC_age<<endl;
  report<<"$TAC"<<endl;
  report<<TAC<<endl;
  report<<"$SIM_F_at_age"<<endl;
  report<<SIM_F_at_age<<endl;
  report<<"$SIM_selectivity"<<endl;
  report<<SIM_selectivity<<endl;
  report<<"$SIM_selectivity_index"<<endl;
  report<<SIM_selectivity_index<<endl;

  report<<"$SIM_ssb"<<endl;
  report<<SIM_ssb<<endl;
  report<<"$SIM_SPR_yr"<<endl;
  report<<SIM_SPR_yr<<endl;    
  report<<"$SIM_BH_REC"<<endl;
  report<<SIM_BH_REC<<endl;
  report<<"$SIM_recruits"<<endl;
  report<<SIM_recruits<<endl;
  report<<"$SIM_abundance_at_age"<<endl;
  report<<SIM_abundance_at_age<<endl;
  report<<"$SIM_abundance_spawn"<<endl;
  report<<SIM_abundance_spawn<<endl;
  report<<"$SIM_init_abund"<<endl;
  report<<SIM_init_abund<<endl;

  report<<"$SIM_biomass"<<endl;
  report<<SIM_biomass<<endl;
  report<<"$SIM_rec_devs"<<endl;
  report<<SIM_rec_devs<<endl;
  
  report<<"$TRUE_catch_prop"<<endl;
  report<<TRUE_catch_prop<<endl;
  report<<"$OBS_catch_prop"<<endl;
  report<<OBS_catch_prop<<endl;
  
  report<<"$TRUE_landings"<<endl;
  report<<TRUE_landings<<endl;
  report<<"$OBS_landings"<<endl;
  report<<OBS_landings<<endl;
  
  report<<"$TRUE_Index"<<endl;
  report<<TRUE_Index<<endl;
  report<<"$OBS_Index"<<endl;
  report<<OBS_Index<<endl;
  report<<"$TRUE_Index_prop"<<endl;
  report<<TRUE_Index_prop<<endl;
  report<<"$OBS_index_prop"<<endl;
  report<<OBS_index_prop<<endl;

  report<<"$MNconst"<<endl;
  report<<MNconst<<endl;
  report<<"$debug"<<endl;
  report<<debug<<endl;
 }
  
RUNTIME_SECTION
  convergence_criteria .001,.0001, 1.0e-4, 1.0e-7
  maximum_function_evaluations 100000
  


