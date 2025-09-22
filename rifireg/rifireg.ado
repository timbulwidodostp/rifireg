*! rifireg version 0.3  Sept-2020, Gawain Heckley 
*! Verson 0.3 corrects for sign error in SRCI weighting function
*! Verson 0.2 adjusts for small sample correction in states correlate function
program rifireg, eclass sort byable(recall) //sort prop(sw)
   version 9, missing
    syntax [varlist] [if] [in] [aw fw iw] [, ///
      /* Models */            	///
      	AGIni		      	///
	GIni		      	///
	ACIndex 		///
	EINdex 			///
	CINdex 			///
	ARCIndex 		///
	SRCIndex 		///
	WIndex 			///
	/* bivariate rank dependent index options */ ///
	RANKV(varname)		///
	UBound(real 0.0)	///
	LBound(real 0.0)	///
	REtain(string)		/// save the estimated RIF values as a vector
      /* Display Options */   ///
      Level(cilevel)	      ///
      BOotstrap		      ///
      reps(integer 49)        ///
      ]

    if !replay() {
  
    //check to see if we're computing something
      if "`agini'`gini'`acindex'`eindex'`cindex'`arcindex'`srcindex'`windex'" == "0" {
	 di in red "must specify statistic to compute"
	 exit 198
      }
  
	//check to see if ranking variable is specified for bivariate statistics
      if (	"`acindex'" 	== "acindex" | ///
		"`eindex'" 	== "eindex" | ///
		"`cindex'" 	== "cindex" | ///
		"`arcindex'" 	== "arcindex" | ///
		"`srcindex'" 	== "srcindex" | ///
		"`windex'" 	== "windex") & "`rankv'" == "" {
	 di in red "must specify rank variable for bivariate statistic"
	 exit 198
      }
      
      //generate and retain options don't mix for bootstrap 
      if ("`bootstrap'" == "bootstrap") ///
	 & "`norobust'`retain'`generate'" != "" {
	 di in red "can't use norobust, generate or retain options with bootstrap"
	 exit 198
      }

      //get the dependent variable
      tokenize `varlist'
      local y `1'
      //get the rest of the vars
      macro shift
      local rest `*'

      //deal with weights, if set as blank, generate a column of 1's
      if "`weight'" == "" {
	 tempvar eweight
         gen `eweight' = 1.0
         local weight "aweight"
         local exp "=`eweight'"
      }

      //get the weight expression without '=' sign
      local exp_no_eq = regexr("`exp'", "=", "")

      //pick which samples we want to use
      marksample touse
      //do not use any samples that are blank
      markout `touse' `varlist' `exp_no_eq'
      //remove colinear explanatory variables
      _rmcoll `varlist' if `touse'
      local varlist `r(varlist)'

     
//see if they want absolute gini
      if "`agini'" == "agini" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifagini
	 }

	 tempname b v  
	 agini "`regname'" "`y'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap AGINI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"agini `regname' `y' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifagini"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')
	 
	 
	 RIF_out `level' //print
	 exit
      }


//see if they want gini (from rifreg.ado)
      if "`gini'" == "gini" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifgini
	 }

	 tempname b v  
	 gini "`regname'" "`y'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap GINI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"gini `regname' `y' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifgini"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')
	 
	 
	 RIF_out `level' //print
	 exit
      }


//see if they want Absolute concentration index
      if "`acindex'" == "acindex" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifACI
}
	
	 if "`rankv'" != "" { //A ranking variable has been given
	 local rank `rankv'
      }
      else { //Assign ranking variable as dependent variable
	 l
	 local rank `y'
      }

	 tempname b v  
	 aconcentrationindex "`regname'" "`y'" "`rank'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap ACI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"aconcentrationindex `regname' `y' `rank' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifACI"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')


   RIF_out `level'
	exit
}


//see if they want Erreygers index
      if "`eindex'" == "eindex" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifEI
}
	
	 if "`rankv'" != "" { //A ranking variable has been given
	 local rank `rankv'
      }
      else { //Assign ranking variable as dependent variable
	 l
	 local rank `y'
      }
	local ah "`lbound'"
	local bh "`ubound'"
	
	 tempname b v  
	 erreygersindex "`regname'" "`y'" "`ah'" "`bh'" "`rank'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap EI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"erreygersindex `regname' `y' `ah' `bh' `rank' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifEI"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')


   RIF_out `level'
	exit
}

//see if they want attainment relative concentration index
      if "`arcindex'" == "arcindex" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifARCI
	}

	 if "`rankv'" != "" { //A ranking variable has been given
	 local rank `rankv'
      }
      else { //Assign ranking variable as dependent variable
	 l
	 local rank `y'
      }
	local ah "`lbound'"
	local bh "`ubound'"
	
	 tempname b v  
	 arconcentrationindex "`regname'" "`y'" "`ah'" "`bh'" "`rank'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap ARCI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"arconcentrationindex `regname' `y' `ah' `bh' `rank' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifARCI"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')


   RIF_out `level'
	exit
}

//see if they want shortfall relative concentration index
      if "`srcindex'" == "srcindex" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifSRCI
}
	
	 if "`rankv'" != "" { //A ranking variable has been given
	 local rank `rankv'
      }
      else { //Assign ranking variable as dependent variable
	 l
	 local rank `y'
      }
	local ah "`lbound'"
	local bh "`ubound'"
	
	 tempname b v  
	 srconcentrationindex "`regname'" "`y'" "`ah'" "`bh'" "`rank'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap SRCI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"srconcentrationindex `regname' `y' `ah' `bh' `rank' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifSRCI"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')


   RIF_out `level'
	exit
}


//see if they want Wagstaff index
      if "`windex'" == "windex" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifWI
}

	
	 if "`rankv'" != "" { //A ranking variable has been given
	 local rank `rankv'
      }
      else { //Assign ranking variable as dependent variable
	 l
	 local rank `y'
      }
	local ah "`lbound'"
	local bh "`ubound'"
	
	 tempname b v  
	 wagstaffindex "`regname'" "`y'" "`ah'" "`bh'" "`rank'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap WI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"wagstaffindex `regname' `y' `ah' `bh' `rank' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifWI"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')


   RIF_out `level'
	exit
}



//see if they want concentrationindex
      if "`cindex'" == "cindex" {
	 if "`retain'" != "" { //they want to choose a name
	    local regname `retain'
	 }
	 else{
	    local regname rifCI
}
	
	 if "`rankv'" != "" { //A ranking variable has been given
	 local rank `rankv'
      }
      else { //Assign ranking variable as dependent variable
	 l
	 local rank `y'
      }

	 tempname b v  
	 concentrationindex "`regname'" "`y'" "`rank'" "`rest'" "`weight'`exp'" "`exp_no_eq'" "`touse'"
	 mat `b' = e(b)
	 mat `v' = e(V) //may get overwritten if bootstrap
	 if "`retain'" == "" { //they want to choose a name
	    drop `regname'
	 }

	 if "`bootstrap'" == "bootstrap" {
	    doBootstrap "Bootstrap CI" "`reps'" "`v'"  "`regname'" "`touse'" ///
	       `"concentrationindex `regname' `y' `rank' "`rest'" "`weight'`exp'" `exp_no_eq' `touse'"' ""
	    eret loc vcetype "Bootstrap"
	 }
	 eret loc cmd    = "rifCI"
	 eret loc title  = "RIF Regression"
	 //post important matrices, this destroys b, V and touse
	 eret post `b' `v', noclear esample(`touse')
	 
}
}

   RIF_out `level'
end



program RIF_out
  args level
  loc tss  = e(rss)  + e(mss)
  loc df_t = e(df_m) + e(df_r)

  #delimit ;
  di ;
  di in gr "      Source {c |}       SS       df       MS" _c;
  di in gr _col(56) "Number of obs = " in ye %7.0f e(N) ;
  di in gr "{hline 13}{c +}{hline 30}" _c;
  di in gr _col(56) "F(" %3.0f e(df_m) "," %6.0f e(df_r) 
	  ") = " in ye %7.2f e(F) ;
  di in gr "       Model {c |} " in ye %11.0g e(mss) 
	 " " %5.0f e(df_m) " " %11.0g e(mss)/e(df_m) _c;
  di in gr _col(56) "Prob > F      = " 
	  in ye %7.4f Ftail(e(df_m),e(df_r),e(F)) ;
  di in gr "    Residual {c |} " in ye %11.0g e(rss) " " 
   	 %5.0f e(df_r) " " %11.0g e(rss)/e(df_r) _c;
  di in gr _col(56) "R-squared     = " in ye %7.4f e(r2)    ;
  di in gr "{hline 13}{c +}{hline 30}" _c;
  di in gr _col(56) "Adj R-squared = " in ye %7.4f e(r2_a) ;
  di in gr "       Total {c |} " in ye %11.0g `tss' 
	 " " %5.0f `df_t' " " %11.0g `tss'/`df_t' _c;
  di in gr _col(56) "Root MSE      = " in ye %7.5g e(rmse) ;
  #delimit cr
  //display table 
  _coef_table, level(`level')
end


program agini, eclass sort
   args	        ///
      rifagini   ///
      y         ///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   gsort `y'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname cum_wt tot_wt cumwy glp pvar cov_y_r index mean twovm if_index
   qui{
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
 // if-index
		corr `pvar' `y' if `touse' , cov
		sca `cov_y_r'	=r(cov_12) 
		sca `index'	=2*`cov_y_r'
		sca `index'=`index'/`tot_wt'*(`tot_wt'-1)  /* Small sample correction */

		gen `b2x'	=-`y'+`mean' if `touse' 
		gen `c2x'	=2*(`y'*(`pvar')-`glp') if `touse' 
		gen `if_index'	=-2*`index' + `b2x' + `c2x' if `touse' 
// rif-index 
        gen `rifagini' = `index' + `if_index' if `touse'      
   }

   qui reg `rifagini' `rest' [`weightexp'] if `touse'
   eret sca agini = `index'

end


program gini, eclass sort
   args	        ///
      rifgini   ///
      y         ///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   gsort `y'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname cum_wt tot_wt cumwy glp pvar rf gini mean twovm rfvar 
   qui{
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
      integ `glp' `pvar' if `touse', gen(`rfvar')
      qui sum `rfvar' if `touse'
      sca `rf'=r(max)
      gen `a2'=`twovm'*`rf' if `touse'
      sca `gini'=1-`a2'
      gen `b2x'=`a2'*`y'/`mean' if `touse'
      gen `c2x'=`twovm'*(`y'*(`pvar'-1)-`glp') if `touse'
      gen `rifgini' = `b2x'+`c2x'+1 if `touse'
   }

   qui reg `rifgini' `rest' [`weightexp'] if `touse'
   eret sca gini = `gini'

end


program aconcentrationindex, eclass sort
   args	        ///
      rifACI   ///
      y         ///
	rank 	///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   qui sort `rank'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname cum_wt tot_wt cumwy glp pvar cov_y_r index mean twovm if_index
   qui{
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
// if-index
		corr `pvar' `y' if `touse' , cov
		sca `cov_y_r'	=r(cov_12) 
		sca `index'	=2*`cov_y_r'
		sca `index'=`index'/`tot_wt'*(`tot_wt'-1)  /* Small sample correction */
		gen `b2x'	=-`y'+`mean' if `touse' 
		gen `c2x'	=2*(`y'*(`pvar')-`glp') if `touse' 
		gen `if_index'	=-2*`index' + `b2x' + `c2x' if `touse' 
// rif-index 
        gen `rifACI' = `index' + `if_index' if `touse' 
   }

   qui reg `rifACI' `rest' [`weightexp'] if `touse'
   eret sca ACindex = `index'

end

program erreygersindex, eclass sort
   args	        ///
      rifEI   ///
      y         ///
	ah	///
	bh	///
	rank 	///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   sort `rank'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname W_EI cum_wt tot_wt cumwy glp pvar cov_y_r acindex mean twovm if_acindex if_EI rif_acindex EI
   qui{
	di `ah' 
	di `bh'
	sca `W_EI'=4/(`bh'-`ah')
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
// if-acindex
		corr `pvar' `y' if `touse' , cov
		sca `cov_y_r'	=r(cov_12) 
		sca `acindex'	=2*`cov_y_r'
		sca `acindex'=`acindex'/`tot_wt'*(`tot_wt'-1)  /* Small sample correction */

		gen `b2x'	=-`y'+`mean' if `touse' 
		gen `c2x'	=2*(`y'*(`pvar')-`glp') if `touse' 
		gen `if_acindex' =-2*`acindex' + `b2x' + `c2x' if `touse' 
// rif-acindex 
        gen `rif_acindex' = `acindex' + `if_acindex' if `touse'
// if-EI
		sca `EI' 		= `W_EI'*`acindex'
		gen `if_EI' 	= `W_EI'*`if_acindex' if `touse'
// rif-EI
		gen `rifEI' = `W_EI'*`rif_acindex' if `touse' 
   }

   qui reg `rifEI' `rest' [`weightexp'] if `touse'
   eret sca Eindex = `EI'

end

program arconcentrationindex, eclass sort
   args	        ///
      rifARCI   ///
      y         ///
	ah	///
	bh	///
	rank 	///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   sort `rank'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname W_ARCI cum_wt tot_wt cumwy glp pvar cov_y_r acindex mean twovm if_acindex rif_acindex if_W_ARCI if_ARCI rif_ARCI ARCI
   qui{
	di `ah' 
	di `bh'
	sca `W_ARCI'=1/(-`ah')
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
	sca `W_ARCI'=1/(`mean'-`ah')
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
// if-acindex
		corr `pvar' `y' if `touse' , cov
		sca `cov_y_r'	=r(cov_12) 
		sca `acindex'	=2*`cov_y_r'
		sca `acindex'=`acindex'/`tot_wt'*(`tot_wt'-1)  /* Small sample correction */
		gen `b2x'	=-`y'+`mean' if `touse' 
		gen `c2x'	=2*(`y'*(`pvar')-`glp') if `touse' 
		gen `if_acindex' =-2*`acindex' + `b2x' + `c2x' if `touse' 
// rif-acindex 
        gen `rif_acindex' = `acindex' + `if_acindex' if `touse'
// if-ARCI
		sca `ARCI' 	= `W_ARCI'*`acindex'
		gen `if_W_ARCI'	=-(`y'-`mean')/((`mean'-`ah')^2) if `touse'
		gen `if_ARCI'=`if_W_ARCI'*`acindex' + `W_ARCI'*`if_acindex' if `touse'
// rif-ARCI
		gen `rifARCI'=`if_ARCI'+`ARCI' if `touse'
   }

   qui reg `rifARCI' `rest' [`weightexp'] if `touse'
   eret sca ARCindex = `ARCI'

end

program srconcentrationindex, eclass sort
   args	        ///
      rifSRCI   ///
      y         ///
	ah	///
	bh	///
	rank 	///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   sort `rank'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname W_SRCI cum_wt tot_wt cumwy glp pvar cov_y_r acindex mean twovm if_acindex rif_acindex if_W_SRCI if_SRCI rif_SRCI SRCI
   qui{
	di `ah' 
	di `bh'
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
	sca `W_SRCI'=-1/(`bh'-`mean')
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
// if-acindex
		corr `pvar' `y' if `touse' , cov
		sca `cov_y_r'	=r(cov_12) 
		sca `acindex'	=2*`cov_y_r'
		sca `acindex'=`acindex'/`tot_wt'*(`tot_wt'-1)  /* Small sample correction */
		gen `b2x'	=-`y'+`mean' if `touse' 
		gen `c2x'	=2*(`y'*(`pvar')-`glp') if `touse' 
		gen `if_acindex' =-2*`acindex' + `b2x' + `c2x' if `touse' 
// rif-acindex 
        gen `rif_acindex' = `acindex' + `if_acindex' if `touse'
// if-SRCI
		sca `SRCI' 	= `W_SRCI'*`acindex'
		gen `if_W_SRCI'	=(`mean'-`y')/((`bh'-`mean')^2) if `touse'
		gen `if_SRCI'=`if_W_SRCI'*`acindex' + `W_SRCI'*`if_acindex' if `touse'
// rif-SRCI
		gen `rifSRCI'=`if_SRCI'+`SRCI' if `touse'
   }

   qui reg `rifSRCI' `rest' [`weightexp'] if `touse'
   eret sca SRCindex = `SRCI'

end

program concentrationindex, eclass sort
   args	        ///
      rifCI   ///
      y         ///
	rank 	///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   sort `rank'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname cum_wt tot_wt cumwy glp pvar rf index mean twovm rfvar 
   qui{
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
      integ `glp' `pvar' if `touse', gen(`rfvar')
      qui sum `rfvar' if `touse'
      sca `rf'=r(max)
      gen `a2'=`twovm'*`rf' if `touse'
      sca `index'=1-`a2'
      gen `b2x'=`a2'*`y'/`mean' if `touse'
      gen `c2x'=`twovm'*(`y'*(`pvar'-1)-`glp') if `touse'
      gen `rifCI' = `b2x'+`c2x'+1 if `touse'
   }

   qui reg `rifCI' `rest' [`weightexp'] if `touse'
   eret sca index = `index'

end

program wagstaffindex, eclass sort
   args	        ///
      rifWI   ///
      y         ///
	ah	///
	bh	///
	rank 	///
      rest      ///
      weightexp ///
      exp_no_eq ///
      touse     

   qui sort `rank'   /*sort in ascending order */

   tempvar a2 b2x c2x
   tempname W_WI cum_wt tot_wt cumwy glp pvar cov_y_r acindex mean twovm if_acindex rif_acindex z2 z3 if_W_WI if_WI rif_WI WI
   qui{
	di `ah' 
	di `bh'
      sum `y' [`weightexp'] if `touse'
      sca `tot_wt'=r(sum_w)
      sca `mean'= r(mean)
      sca `twovm'=2/`mean'
	sca `W_WI'=(`bh'-`ah')/((`bh'-`mean')*(`mean'-`ah'))
      gen `cum_wt' = sum(`exp_no_eq') if `touse'

      gen double `cumwy' = sum(`y'*`exp_no_eq') if `touse'
      gen double `glp'=`cumwy'/`tot_wt' if `touse'
      gen double `pvar' = `cum_wt'/`tot_wt'  if `touse'
// if-acindex
		corr `pvar' `y' if `touse' , cov
		sca `cov_y_r'	=r(cov_12) 
		sca `acindex'	=2*`cov_y_r'
		sca `acindex'=`acindex'/`tot_wt'*(`tot_wt'-1)  /* Small sample correction */
		gen `b2x'	=-`y'+`mean' if `touse' 
		gen `c2x'	=2*(`y'*(`pvar')-`glp') if `touse' 
		gen `if_acindex' =-2*`acindex' + `b2x' + `c2x' if `touse' 
// rif-acindex 
        gen `rif_acindex' = `acindex' + `if_acindex' if `touse'
// if-WI
		sca `WI' 	= `W_WI'*`acindex'
		gen `z2' 	= ((`bh'+`ah'-2*`mean')*(`y' - `mean')) if `touse'
		gen `z3'	= ((`bh'-`mean')*(`mean'-`ah')) if `touse'
		gen `if_W_WI'	= (-(`bh'-`ah')*`z2')/(`z3'*`z3') if `touse'
		gen `if_WI'	= `if_W_WI'*`acindex' + `W_WI'*`if_acindex' if `touse'
// rif-WI
		gen `rifWI'=`if_WI'+`WI' if `touse'
   }

   qui reg `rifWI' `rest' [`weightexp'] if `touse'
   eret sca Windex = `WI'

end

//small program to perform bootstrap
program doBootstrap, eclass
   args title /// title of the bootstrap, displayed to user
        reps  /// number of repetitions
	v     /// array that will be overwritten w/ variance
	drops /// variables that are dropped *AFTER* each run of cmds
	touse /// variables that are to be used
	cmd1  /// first command that is run, followed by 2nd
	cmd2
   //see if we have enough matrix space to run all the iterations
   if `reps' > `=c(max_matsize)' {
      di in red "matsize is not large enough for the number of repetitions requested"
      di in red "please increase it (ie: set matsize `reps')"
      exit 198
   }
   tempname bootstrap_est sav
   _estimates hold `bootstrap_est'
   _dots 0, title(`title') reps(`reps')
   forval i=1/`reps' {
      preserve
      bsample if `touse'
      `cmd1'
      `cmd2'
      cap drop `drops'
      mat `sav' = nullmat(`sav') \ e(b)
      restore
      _dots `i' 0
   }
   di
   mata : rifmean_variance("`sav'", "`v'") //overwrite the variance matrix
   mat drop `sav'
   _estimates unhold `bootstrap_est'
end

version 9
mata:
mata set matastrict on

void rifmean_variance(string scalar M, string scalar V) 
{
//   st_matrix(V)
   MV = meanvariance(st_matrix(M), 1)
//   means = MV[1,.]
//   means
   MV = diag(diagonal(MV[|2,1\.,.|])')
   //MV
   st_replacematrix(V, MV)
}
end
* $Id: rifreg.ado,v 1.13 2009/04/08 22:34:13 mgevaert Exp mgevaert $
