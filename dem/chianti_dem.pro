;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. See www.chiantidatabase.org
;                   
; Name        : CHIANTI_DEM
;     		          
; Purpose     : Calculates the Differential Emission Measure DEM(T) using 
;		the CHIANTI database, from a given set of observed lines.
;		Constant pressure or density can be used.
;
; Category    : diagnostic analysis
;               
; Method :      This routine has several options, all in the form of keywords.
; 		First, the input file with the observed fluxes is read.
;		The first time you use this routine
;		you'll have to do the calculation of the contribution
;		functions G(T), so the routine GET_CONTRIBUTIONS will come 
;		into play. You'll have to specify the value of the 
;		pressure or density, and  you'll be asked to select an  
;		ionization equilibrium file.
;		GET_CONTRIBUTIONS searches the CHIANTI database 
;		for all the lines corresponding to the observed 
;		lines, i.e. that lie in a OBS_WVL(i) +/- DELTA_LAMBDA_OBS(i) 
;		interval centered on the observed wavelength OBS_WVL(i).
;		The routine calculates the C(T) values (G(T)=Ab(element)*C(T))
;		for the temperature interval and steps listed in the
;		ionization equilibrium file you have chosen.
;
;		You can either select a constant pressure OR a constant 
;		density for all the lines; if you select a constant pressure,
;               for each ion the contribution function is calculated at an 
;               electron density N_e equal to the ratio of the pressure 
;               and the temperature of maximum ionization fraction:  
;               C=C( T, N_e= P/T_ion )  
;               The C(T) values are stored by GET_CONTRIBUTIONS in the output 
;		file OUTPUT.CONTRIBUTIONS that can be used later to calculate
;		the DEM, changing various parameters,
;		without having to start again and read the CHIANTI database,
;		which can take some time. 
;
;               Additionally, an IDL structure with the G(T) of the
;               lines is also saved in a genx file.
;
;		In the case no CHIANTI lines corresponding to an observed
;		line are found, the routine writes the wavelength of the line
;		(to be excluded from the fit) in the array
;		EXCLU_OBS_WVL_NO_TEO. The lines with no theoretical 
;		counterparts are then automatically  excluded from the fit by 
;	 	CHIANTI_DEM. You might consider the possibility to start again
;		incrementing the DELTA_LAMBDA_OBS, to see if there are 
;		CHIANTI lines in the vicinity.
;
;		Note: if you want to exclude some of the observed lines from 
;		the fit, you just have to use the keyword EXCLUDE_OBS_WVL, 
;		BUT  GET_CONTRIBUTIONS will store anyway the results (if any)
;		in the C(T) file.
;
;		After having excluded the lines in EXCLUDE_OBS_WVL, 
;		any *.abund file present in the CHIANTI database or in 
;		the working directory can be selected, and eventually edited,
;		if you like to change some abundances. 
;		Then the $G(T)$ are calculated, multiplying each theoretical
;		line by the abundance factor.
;
;               Next, the G(T) of the lines are interpolated over a
;               grid of temperatures given as input with the keywords
;
;               MIN_logT:   minimum temperature (log T) for the DEM
;                            calculation. If not defined, the minimum
;                            value is taken from the ionisation
;                            equilibrium table. 
;               MAX_logT:   maximum temperature (log T) for the DEM
;                            calculation.  If not defined, the minimum
;                            value is taken from the ionisation
;                            equilibrium table.
;               DT_logt:      temperature step (log T) for the DEM
;                            calculation. If not defined, a
;                            default=0.05 is assumed.
;
;               Make sure that the requeste minimum and maximum are
;               within the temperature array in the ioneq file. 
;
;               Then the
;		theoretical lines contributing to each blend are sorted by
;		intensity and then their G(T) can be plotted if the keyword
;		PLOT_GT is activated. It is recommended to do this the first 
;		time, to check if there are some observed lines  
;		blended with lines of other elements, in which case it is
;		better to exclude them with a second run (if you are not 
;		sure about the abundances).
; 
;               Then the total G(T) are interpolated over the
;               requested temperature grid.
;
;		Then the G(T) for each blend are summed and optionally
;		plotted. The EM Loci curves are also plotted
;		optionally. 
;
;		Then  the fit starts.
; 
;               NEW OPTIONS (as of 2017):
; 
;               1) DO_XRT_DEM:  keyword to ask the routine to run this
;                            inversion (default)

;               By default, the routine runs the routine MPFIT_DEM,
;               which is a modification of the routine
;               XRT_DEM_ITER_NOWIDGET, written by M. Weber, and which
;               is part of the XRT_DEM  package, and available within
;               the SSW XRT path. 
;
;               The G(T) are resampled using the IDL routine INTERPOL
;               onto the chosen temperature grid.
;
;
;               As an input, one can pass a CHIANTI DEM file name with
;               the  spline nodes and values, OR provide the spline
;               nodes and values:
;
;               SPL_LOGT:    logt T values for the spline nodes. Note
;               that the number should not be larger than the number of
;               input lines -1.
;
;               SPL_LOGDEM: log DEM input values. If not defined, a
;               value of 20. will be used.
;
;               Note: if the SPL_LOGT values are not given as input,
;               the program will select an equally-spaced number of
;               spline nodes within the given minimum and maximum
;               temperatures. 
;
;               Some checks have been included, to avoid problems.
;
;               **** IMPORTANT NOTE **** 
;
;               As the line intensities are calculated with a sum,
;               their values depend somewhat on the temperature
;               step. A value of 0.05 in log T (the default) or less
;               is recommended.  
;
;
;               MIN_LIMITS,  MAX_LIMITS: minimum and maximum limits for
;               the log DEM values at the spline temperatures.
;               Note: if no limits are given, by default the routine
;               applies a minimum value for the log DEM=17.
;
;
;
;               OUTPUTS:
;
;                1) output+'_mpfit_dem.dem': the DEM as a CHIANTI ascii 
;                file.
;                2) output+'_mpfit_dem.save': an IDL save file with:
;                input: the input structure for XRT_DEM_ITER_NOWIDGET
;                logT_out: the log T
;                log_dem_out: the log DEM 
;                log_dem_mciter: (optional) the results of the Monte
;                Carlo runs
;                Users can restore the save files and re-run the
;                XRT_DEM_ITER_NOWIDGET later on independently of
;                CHIANTI_DEM. 
;
;
;               2) As an option, the routine can run DATA2DEM_REG, a routine
;               that recovers the DEM using a GSVD approach, detailed
;               in Hannah & Kontar A&A 539, A146 2012. 
;               Users must download the IDL suite of routines, found
;               in http://www.astro.gla.ac.uk/~iain/demreg/
;               and add them to the IDL path.
;               A series of keywords are passed to the code, see the
;               header of  DATA2DEM_REG. 
; 
;		A series of input parameters can change the 
;		result (DEM), especially the number of temperatures
;		and the temperature range.
;
;               INPUT:
;
;               DO_DEMREG: keyword to ask the routine to run
;                          DATA2DEM_REG. 
;
;               NT_DEMREG: Number of temperatures  for the DEM
;                          calculation (default=20). Note that this
;                          number must be larger than the number of
;                          input lines+1.
;
;               OUTPUT:
;                1) output+'_demreg.dem': the DEM as a CHIANTI ascii 
;                file.
;
;                2) output+'_demreg.save': an IDL save file with 
;                REG: a structure containing all the input and output
;                results of the routine.
;                Users can restore the save file and re-run 
;                DATA2DEM_REG later on independently of CHIANTI_DEM.  
;
; 
;               3) As an option, the routine can run the the PINTofALE
;               command-line function MCMC_DEM(), which runs a
;               Markov-Chain Monte-Carlo algorithm on a set of line
;               fluxes and returns an estimate of the DEM that
;               generates the observed fluxes.  
;               Users should download the package and add the routines
;               to the IDL path. For more information, see
;               http://hea-www.harvard.edu/PINTofALE/
;               Note: the MCMC_DEM() has many keywords. Only some are
;               passed to this routine.
;
;               INPUT:
;
;               DO_MCMC: keyword to ask the routine to run MCMC
;
;               Note: all the keywords accepted by MCMC_DEM() are
;               accepted. they are passed to the routine via the
;               _extra keyword. 
;
;               OUTPUT:
;               1) output+'_mcmc.dem': the DEM as a CHIANTI ascii 
;                file.
;               2) output+'mcmc.save' an IDL save file with all the
;               input and output results, keyword parameters, etc.
;               Users can restore the save file and re-run 
;               MCMC_DEM()  later on independently of CHIANTI_DEM.  
;
;
; Use         : IDL> chianti_dem,output='test_obs',file_input='test_obs',$
;				pressure=3.e15
;
;
; Examples    : 
;		Assume you have a file input 'test_obs' like this:
;
; 171.114    4811.0    1443.0 0.25   Fe IX
; 174.604    4005.0    1202.0 0.25    Fe X
; 180.448    3877.0    1163.0 0.25 Fe XI bl Fe X
; 195.149    3443.0    1033.0 0.25  Fe XII
; 201.246    1091.0     327.2 0.25 Fe XIII
; 211.345    2100.0     630.1 0.25  Fe XIV
; 284.153    4221.0    1266.0 0.25 Fe XV bl
; 315.173     416.4     124.9 0.25 Mg VIII
; 319.906     419.1     125.7 0.25 Si VIII
; 403.370      99.4      29.8 0.25 Mg VI bl Ne VI 
; 434.894     127.0      38.1 0.25 Mg VII bl ? 
; 445.742      41.5      12.4 0.25  S  XIV
; 465.259     349.5     104.9 0.25  Ne VII
; 661.834 3.932e+01 1.179e+01 0.25   unid.  
; 685.790      73.2      22.0 0.25  N  III
; 703.840     176.5      53.0 0.25  O  III
; 760.300      60.1      18.0 0.25    O  V
; 765.110     143.4      43.0 0.25   N  IV
;		
;		IDL> chianti_dem,output='test_obs',file_input='test_obs',$
;		   pressure=1.e15,cut_gt=1e-30,/plot_gt
;
;		After having selected the  ionization file,
;		the C(T) (with MAX(C(T)) gt 1e-30)  are stored in the file
;		'test_obs.contributions'. Then select one of the abundance 
;		files. 
;		Have a look at the plots of the  G(T), and annotate
;		if there is a line you want to exclude, let's say the second.
;		Have a look at the DEM obtained ('test.dem') and at 
;		the details contained in the file 'test.general'. 
;		Maybe there is another line you want to exclude, let's say 
;		the last one. To constrain e.g. the temperature bounds run:
;
;		So run
;		IDL> chianti_dem,output='test_2',file_input='test_obs',$
;		file_gt='test_obs.contributions', xrt_min_t=5.5,xrt_max_t=6.6
;
;               2) To run the DATA2DEM_REG, do e.g. 
;               IDL> chianti_dem,output='test',file_input='test_obs',$
;               file_gt='test.contributions',/do_demreg,demreg_logt_min=5.5,$ 
;               demreg_logt_max=6.6, nt_demreg=20 
;
;               3) to run MCMC_DEM, do e.g.:
;               IDL> chianti_dem,output='test',file_input='test_obs',$
;               file_gt='test.contributions', /do_mcmc, mcmc_logt_step=0.1,$
;               mcmc_logt_max=6.6, mcmc_logt_min=5.5
;
;    
; Inputs      : many all in form of keywords, se above. The required ones are 
;		OUTPUT and FILE_GT (or  PRESSURE/DENSITY)
;               
;               
; Opt. Inputs : various... see above
;               
; Outputs     : OUTPUT.CONTRIBUTIONS  
;
;		Created only if the keyword FILE_GT is NOT set. 
;		Is the file where all the contribution  functions G(T) are 
;		stored. In the first two lines  the ionization equilibrium 
;		file name, and the constant value of pressure or density 
;		adopted are reported. Then for each line you have reported  
;               the observed wavelength, the theoretical one, the element and
;		ionization stage, then the C(T) values. At the end the 
;		specification for each transition.
;
;		OUTPUT.DEM
;		Is the file where the log T and log DEM values are 
;		written, with a format suitable 
;		as input for the CH_SS procedure,that calculates the 
;		synthetic spectrum. At the end some info on how it was 
;		calculated are printed.
;
;		OUTPUT.GENERAL
;		Is the file where general information is stored.
;		The abundance file, the ionization equilibrium file and the
;		constant value of pressure or density  used are reported. 
;		Then there is one line for each
;		observed line, with the provisional identification, the 
;		observed wavelength, the observed flux, the theoretical one
;		(corresponding to the DEM), the error on the flux,
;		the square of the difference between the theoretical and the 
;		observed fluxes divided by the error (this number should be 
;		close to zero if the line is well reproduced), and finally
;		the ratio of the theoretical flux versus the observed one 
;		(which should be close to 1).
;		After this line, there is one line per each theoretical line
;		contributing to the blend, with the identification, the 
;		theoretical wavelength, the configuration and terms, and the
;		contribution to the total theoretical flux (in percentage) 
;		of each line in the blend.
;
;	
; Opt. Outputs:
;		An abundance file with the modifications inserted.
;	
;		Postscript files of the G(T).
;	
;		A postscript file with the DEM (OUTPUT.DEM.PS)
;		
;
; Other Keywords    : 
;
;
;	ARCSEC: 
;		optional. If set, it means that the intensities in the input
;		file are per arcsec-2 .
;		These intensities are then  converted to 
;		 sterad-1 .
;
;	CUT_GT:	
;		optional. If set, only the those theoretical lines that
;		have a MAX(C(T))  greater than the value set, are kept; 
;		it is useful to set this value in order to reduce the number 
;		of lines in the file where the C(T) are stored.;  
;		if not set, a default value of 1e-30 is adopted.
;
;	DEM_FILE:
;		optional.If set (,/DEM_FILE) you have to choose a DEM file to 
;		be used as a start, instead of the default constant value of 
;		10.^22.
;		You can either choose one of the files in the CHIANTI database
;		or any you have in the working directory. 
;		The values in the file are marked as crosses, the mesh points
;		are marked with triangles.
;
;	DENSITY : 
;		the value of the density (Ne). Required if you do NOT have
;		already the contribution  functions G(T). 
;
;	EXCLUDE_OBS_WVL:
;		optional.
;		If set, you can  exclude some of the observed lines from 
;               the fit. Note that even if you set this keyword and run 
;		GET_CONTRIBUTIONS all the theoretical lines found corresponding
;		to all the lines in the input file are written in the C(T) 
;		file. It is only in the fit that the lines are excluded.
;
;	FILE_GT:
;		optional.
;		If NOT set, the routine GET_CONTRIBUTIONS is called.
;               If set, it has to specify the name of the file created by 
;		GET_CONTRIBUTIONS, where all the contribution  functions G(T) 
;               are stored. In the first two lines the ionization equilibrium 
;		file name, and the value of the pressure or density 
;               adopted is reported. Then for each line you have reported  
;               the observed wavelength, the theoretical one, the element and
;               ionization stage, then the C(T) values. At the end the 
;               specification for each transition.
;
;	FILE_INPUT:
;		optional.
; 		if set, you are not requested to select the observation file
;		using a widget-type search.
;		The input file  must contain 5 columns, unformatted:
;		1)the observed wavelength (A)
;		2)the observed flux in erg cm-2 s-1 st-1
;		3)the corresponding error on the flux in erg cm-2 s-1 st-1
;		4)half the width (A) of the range (centered on the observed 
;		  wavelength) where you want to look for the corresponding 
;		  theoretical lines. A value of HWHM or more would do.
;		5)The identification, written as string (max 20 characters)
;
;
;	N_MATCHES:   
;		optional.          
;		In the unlikely event that more than 1000 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;
;	OUTPUT  :
;		required.
;	  	It is the output name. Suffixes will be added when creating 
;		the various outputs.
;
;	PHOT:
;		optional.
;		If set, it means that in the input file the intensities
;		are in photons instead of ergs. 
;
;	PLOT_GT:
;		optional.
;		If set (,/PLOT_GT),  plots of the  G(T) for each 
;		observed line not excluded are created.
;
;	PRESSURE:     		
;		the value of the pressure (Ne T). Required if you do NOT have
;		already the contribution  functions G(T).
;
;	VERBOSE:
;		optional. Set to get various messages and the details of the 
;		result.
;
;       ABUND_NAME:
;		optional. The name of the CHIANTI abundance file.
;
;       RADTEMP   The blackbody radiation field temperature (default 6000 K).
;
;       RPHOT    Distance from the centre of the star in stellar radius units.
;                I.e., RPHOT=1 corresponds to the star's surface. (Default is
;                infinity, i.e., no photoexcitation.)
;
;       RADFUNC         The name of a user-defined function that will generate
;                       a radiation spectrum as a function of temperature. 
;                       This radiation field will replace the black-body that
;                       is assumed when using the RADTEMP keyword in the call
;                       to pop_solver.
;
; Calls       : GET_CONTRIBUTIONS
;		ZION2SPECTROSCOPIC
;		PRINT2D_PLOT
;               MPFIT_DEM
;               Other IDL routines external to the CHIANTI software.
;		
;
; Restrictions: 
;		In the unlikely event that more than 1000 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;               
; Side effects: None known yet.
;               
; Category    : spectrum
;               
; Prev. Hist. :
;       Written by Ken Dere (NRL) as part of the CHIANTI package 
;       in collaboration with Brunella Monsignori Fossi, Enrico Landi
;       (Arcetri Observatory, Florence), Helen Mason and Peter Young
;       (DAMTP, Cambridge Univ.). Incorporated into the CDS software.  
;
; Written     : 
;       V. 1.0  5 November  1997 Giulio Del Zanna (GDZ), 
;	UCLAN  (University of Central Lancashire, UK)
;
;
; Modified    : Removed the print2d_plot subroutine. Increased the default value
;               of N_MATCHES from 20 to 50.  Changed way to deal with xuvtop.
;               GDZ, 31-Oct-2000
;
; Version     : 2.0 GDZ, DAMTP, University of Cambridge, 31-Oct-2000
;
;              V.3, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;              v.4, 11 Mar 2014, Giulio Del Zanna, major re-write,
;              including three different inversion methods.
;              Test version.
;
;              v.5, 9 Apr, GDZ,
;              fixed the definition of the G(T)=AxC(T), line 1006 
;
;              v.6, 25 Jun 2014, GDZ, 
;              added two wrapper routines to call two of the inversion
;              methods, so the routine compiles even if the programs
;              are not available.
;
;              v.7, 14 Jul 2014, GDZ
;              changed output of DEMREG so the positive solutions are
;              printed in the dem file.
;
;             v8, 29 Apr 2015, GDZ
;             added EM loci plot and various outputs.
;
;             Version 9, 3 June 2017, GDZ. Replaced the SSW XRT_DEM routine
;             with a simplified version, which additionally can take
;             as input limits to the DEM values, the temperatures and
;             input values of the spline nodes.
;
;             Version 10, 6-Jul-2018, GDZ.
;             Added some error checking and changed names of
;             some keywords; added a plot of the DEM and ask if you want to
;             redo the fit.
;
;             v.11, 10 Oct 2018, GDZ.
;             Changed the default  DT_logt, now 0.05
; 
;             v.12, 10-Nov-2020, GDZ
;             Significant modifications. Added checks on
;             inputs. Removed the COMMON blocks. added rphot, radtemp,
;             RADFUNC and some output. Changed default n_matches and
;             fixed a few minor bugs.
;
;
; VERSION     :  12, 10-Nov-2020
;-

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;                              MAIN PROCEDURE:
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;--------------

pro chianti_dem,output=output,file_input=file_input,pressure=pressure,$
                density=density,cut_gt=cut_gt,plot_gt=plot_gt,$
                n_matches=n_matches,file_gt=file_gt,n_iter=n_iter,dchisq_m=dchisq_m,$
                exclude_obs_wvl=exclude_obs_wvl,dem_file=dem_file,arcsec=arcsec,phot=phot,$
                do_xrt_dem=do_xrt_dem, dt_logt=dt_logt, min_logt=min_logt,max_logt=max_logt,$ 
                spl_logt=spl_logt, spl_logdem=spl_logdem, $
                min_limits=min_limits, max_limits=max_limits,$
                do_mcmc=do_mcmc,  do_demreg=do_demreg,$
                abund_name=abund_name, _extra=_extra,verbose=verbose,$
                rphot=rphot, radtemp=radtemp, RADFUNC=RADFUNC

;-------------

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;	ON_ERROR, 2

; by default do the XRT_DEM inversion unless the other options are set

if not keyword_set(do_mcmc) and not keyword_set(do_demreg) then begin 
   if not keyword_set(do_xrt_dem) then do_xrt_dem=1
endif 

!p.background=255
!p.color=0
!p.charsize=1
!p.multi=0


;some keyword checking......
;-----------------------------------

if n_elements(output) eq 0  then begin
   print,' '
   print,'You have to define an output core name (e.g. output="test") !!'
   print,'Returning.. '
   return
endif

if  keyword_set (dem_file) and n_elements(spl_logt) gt 0 then begin
   print,' '
   print,'You have to define either the spline values or ask for a DEM file which has the spline nodes.'
   print,'Returning.. '
   return
endif

if n_elements(spl_logt) gt 0 then begin 

if n_elements(spl_logdem) gt 0 then begin
 if  n_elements(spl_logdem) ne n_elements(spl_logt) then $
    message,'Error, No. of spline DEM values not equal than the No. of Temperatures'

endif 
endif 


if n_elements(file_gt) eq 0  then begin

   if n_elements(pressure) eq 0 and n_elements(density) eq 0 then begin
      print,' You have  to define either a  constant pressure'+$
            '(e.g. pressure=1.e16)'
      print,' or a constant electron density (e.g. density=1.e-9) !!!'
      print,' '
      return
   endif

   if keyword_set(pressure)  and keyword_set(density)  then begin
      print,' You  have  to decide if you want the G(T) calculated '+$
            'for constant pressure or not ! '
      print,'(e.g.  pressure=1.e16) or '+$
            'for constant electron density (e.g.  density=1.e-9) !!!'
      print,' '
      return
   endif

   if n_elements(cut_gt) eq 0 then cut_contrib=1.e-30 else $
      cut_contrib=cut_gt

   if  n_elements(exclude_obs_wvl) ne 0 then $
      print,'If you DO NOT want to calculate the G(T) functions for some '+ $
            'of the observed lines, you have to cancel them from the '+$
            'input file and start again !'
   print,' '

endif else begin
   
   if keyword_set(pressure) then $
      print,'You defined the pressure,but this value will not be '+$
            'considered.The value stored in the save file is the correct one'
   print,' '
   if keyword_set(density)  then $
      print,'You defined the density, but this value will not be '+$
            'considered.The value stored in the save file is the correct one'
   print,' '

   if keyword_set(cut_gt) then begin
      print,'You cannot select a min G(T) at this point !'
      print,'your selection is ignored ....'
      print,' '
   endif

endelse

; expected maximum number of wavelength matches:
;-----------------------------------------------
if  keyword_set(n_matches) then n_ch=n_matches else n_ch=10000

IF n_elements(rphot) EQ 0 THEN BEGIN 
   photoexcitation =0

   IF n_elements(RADFUNC) eq 1 THEN BEGIN
      err_msg ='% Error ! if you define RADFUNC you need to define RPHOT   - EXIT'
   print, err_msg
   return
   END
  
ENDIF ELSE BEGIN
   
   photoexcitation = 1 
IF n_elements(radtemp) EQ 0 and n_elements(RADFUNC) eq 0 THEN BEGIN 
   radtemp=6d3
   ;this is used in the common block:
radt=radtemp
   IF  keyword_set(verbose) THEN print,'setting a default radiation temperature of 6000 K'
ENDIF ELSE if n_elements(radtemp) EQ 1 and n_elements(RADFUNC) eq 0 THEN begin
   radtemp=float(radtemp)
   ;this is used in the common block:
radt=radtemp
endif else if n_elements(radtemp) EQ 1 and n_elements(RADFUNC) eq 1 THEN begin
  err_msg ='% Error ! you cannot define both RADTEMP and RADFUNC - EXIT'
   print, err_msg
   return
ENDIF 
ENDELSE  



;------------------------------------------------------------
;define the CHIANTI top directory:
;--------------------------------
defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
   message, 'system variable !xuvtop must be set '
xuvtop =!xuvtop


;pick up the file with the input intensities:
;--------------------------------------------
if  keyword_set(file_input) then file_input_obs=file_input else begin
   file_input_obs=pickfile(title='Select the  data file with the fluxes :')
endelse

;read the input file
;-------------------

openr,lur,file_input_obs,/get_lun
;
input=fltarr(4) & str=' '

readf,lur,input,str
obs_wvl=input(0)
obs_int=input(1)
obs_sig=input(2)
obs_delta_lambda=input(3)
obs_id=strtrim(str,2)

while not eof(lur) do begin

   readf,lur,input,str
   obs_wvl=[obs_wvl,input(0)]
   obs_int=[obs_int,input(1)]
   obs_sig=[obs_sig,input(2)]
   obs_delta_lambda=[obs_delta_lambda,input(3)]
   obs_id=[obs_id,strtrim(str,2)]
;
endwhile
;
free_lun,lur
;
n_obs=n_elements(obs_wvl)

obs_wvl=1000.*obs_wvl & obs_wvl=round(obs_wvl)
obs_wvl=obs_wvl/1000.

if  keyword_set(phot) then begin 

   obs_int=obs_int *1.9866e-8/  obs_wvl
   obs_sig=obs_sig *1.9866e-8/  obs_wvl

   print,' '
   print,' The intensities have been converted from photons to ergs '
   print,' '

endif


if  keyword_set(arcsec) then begin 
;***********************************

   obs_int=obs_int /(!PI/(180.*60.^2))^2
   obs_sig=obs_sig /(!PI/(180.*60.^2))^2

   print,' '
   print,' The intensities have been converted from  arcsec-2 to  sterad-1 '
   print,' '
endif

IF  keyword_set(verbose) THEN BEGIN 
   print,' '
   print,'This are the values in the input file : '
   print,' '

;
   for i=0,n_obs-1 do begin
      print,i,obs_wvl(i),obs_id(i),obs_int(i),obs_sig(i),obs_delta_lambda(i),$
            format='(i5,f10.3,a20,2f10.1,1x,f4.2)'
   endfor

   print,' '
ENDIF


;In case you want to exclude a line from the fit and you typed a 
;---------------------------------------------------------------
; wrong wavelength:
;-----------------

if  n_elements(exclude_obs_wvl) ne 0 then begin
   for i=0,n_elements(exclude_obs_wvl)-1 do begin
      b=where (exclude_obs_wvl(i) eq obs_wvl)
      if b(0) eq -1 then begin
         print,'You have typed a wrong wavelength : ',$
               exclude_obs_wvl(i),' START AGAIN'
         return
      endif
   endfor
endif  


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(file_gt) then begin

;start reading the contribution file
;-----------------------------------
   ioneq_name=' '
   const_net='' 
   lambda_obs=fltarr(100000)

   openr,lucon,file_gt,/get_lun
   readf,lucon,ioneq_name
   readf,lucon,const_net

   a=0.& i=0L 
   while(not eof(lucon)) do begin
      readf,lucon,a
      lambda_obs[i]=a
      i=i+1
   endwhile
   free_lun,lucon
   
   
; read the input grid of temperatures:
   read_ioneq,ioneq_name,ioneq_logt,ioneq,ioneq_ref
   
   n_dem_temp=n_elements(ioneq_logt)
   
; inizialise the arrays:                
;----------------------
   
   ch_wvl=fltarr(n_ch,n_obs)
   ch_term=strarr(n_ch,n_obs)
   
   ch_z=intarr(n_ch,n_obs)
   ch_ion=intarr(n_ch,n_obs)
   ch_id=strarr(n_ch,n_obs)

;; ;this array will have the C(T) 
;; ;-----------------------------
   ch_contr_wa=fltarr(n_dem_temp,n_ch,n_obs)

;; ;this array will have the G(T)=C(T)*abundances !!!
;; ;------------------------------------------------- 
   ch_contr=fltarr(n_dem_temp,n_ch,n_obs)

   ch_n_contr=intarr(n_obs)
   
   nlines_teo=i

   lambda_obs= lambda_obs(0:nlines_teo-1)

;check if there are observed lines NOT present in the C(T) file
;--------------------------------------------------------------
   
   exclu_obs_wvl_no_teo=fltarr(n_obs)
   i_ex=0
   for iobs=0,(n_obs-1) do begin
      index=where (lambda_obs eq obs_wvl(iobs) )
      if  index(0) ne -1 then begin
         ch_n_contr(iobs)=n_elements(index)
      endif else begin
         exclu_obs_wvl_no_teo(i_ex)=obs_wvl(iobs)
         i_ex=i_ex+1
      endelse
   endfor

;clean up the zeros:
   n=where(exclu_obs_wvl_no_teo ne 0.)
   if n(0) eq -1 then exclu_obs_wvl_no_teo=0. else $
      exclu_obs_wvl_no_teo=exclu_obs_wvl_no_teo(n)
   m=where(ch_n_contr ne 0)

;this array is essential for reading the file of the C(T)
;--------------------------------------------------------
   ch_n_contr=ch_n_contr(m)


;EXCLUDE THE LINES WITH NO CHIANTI COUNTERPARTS FROM THE FIT:
;----------------------------------------------------------------

   if exclu_obs_wvl_no_teo(0) ne 0. then begin

      print,'The following lines do not have CHIANTI lines '+$
            'in the C(T) file and have been excluded : '
      print,exclu_obs_wvl_no_teo
      print,' '

;create the list of indexes corresponding to the excluded lines:
      out=intarr(n_elements(exclu_obs_wvl_no_teo))
      i=0
      for j=0,n_elements(exclu_obs_wvl_no_teo)-1 do begin
         out(i)=where( (exclu_obs_wvl_no_teo(j) ne obs_wvl) eq 0)
         i=i+1
      endfor

      in=indgen(n_obs)
;create the list of indexes corresponding to the lines we keep for the fit:
      remove,out,in

;rearrange the arrays:
      n_obs=n_obs-n_elements(exclu_obs_wvl_no_teo)
      obs_wvl=obs_wvl(in)
      obs_id=obs_id(in)
      obs_int=obs_int(in)
      obs_sig=obs_sig(in)
      obs_delta_lambda=obs_delta_lambda(in)

      ch_wvl=fltarr(n_ch,n_obs)
      ch_id=strarr(n_ch,n_obs)
      ch_z=intarr(n_ch,n_obs)
      ch_ion=intarr(n_ch,n_obs)
      ch_contr_wa=fltarr(n_dem_temp,n_ch,n_obs)
      ch_term=strarr(n_ch,n_obs)


   endif
;---------------------------------------------------------------------------
;read the contribution file
;-----------------------------------

   openr,lucon,file_gt,/get_lun

;readf,lucon,abund_name
   readf,lucon,ioneq_name                
   readf,lucon,const_net

   b0=0.&b1=0.&b2=' '& b3=1 &b4=1 &b5=fltarr(n_dem_temp) &b6=' '
   
   format='(f10.3,1x,f10.3,a15,2i5,1x,'+trim(n_dem_temp)+'e11.3,a40)'

   for iobs=0,n_obs-1 do begin
      n_list=ch_n_contr(iobs)
      if n_list ne 0 then begin 
         for ilist=0,n_list-1 do begin
            readf,lucon,b0,b1,b2,b3,b4,b5,b6,$
                  format=format

            obs_wvl(iobs)=b0 &ch_wvl(ilist,iobs)=b1 
            ch_id(ilist,iobs)=b2
            ch_z(ilist,iobs)=b3 &ch_ion(ilist,iobs)=b4 
            ch_contr_wa(*,ilist,iobs)=b5
            ch_term(ilist,iobs)=b6

         endfor
      endif
   endfor
   free_lun,lucon

;;end of the if you  have the contribution functions already 
endif    else begin 

;----------------------------------------------------------------------
;--- NOW WE START WITH GET_CONTRIBUTIONS ------------------------------
;----------------------------------------------------------------------

input={obs_int:obs_int,obs_sig:obs_sig,obs_wvl:obs_wvl,$
     obs_id: obs_id,obs_delta_lambda:obs_delta_lambda}


   get_contributions,input, out,output_name=output,  density=density, pressure=pressure,$
                     cut_contrib=cut_contrib,n_ch=n_ch, rphot=rphot, radtemp=radtemp,RADFUNC=RADFUNC

ioneq_name = out.ioneq_name 
ioneq_logt=out.ioneq_logt
const_net= out.const_net
exclu_obs_wvl_no_teo=out.exclu_obs_wvl_no_teo
ch_wvl= out.ch_wvl
ch_id= out.ch_id
ch_z= out.ch_z
ch_ion= out.ch_ion
ch_contr_wa= out.ch_contr_wa
ch_term= out.ch_term
ch_n_contr= out.ch_n_contr

; read the input grid of temperatures:
;   read_ioneq,ioneq_name,ioneq_logt,ioneq,ioneq_ref
   
   n_dem_temp=n_elements(ioneq_logt)
   
   
;EXCLUDE THE LINES WITH NO THEORETICAL COUNTERPARTS FROM THE FIT:
;----------------------------------------------------------------------

   if exclu_obs_wvl_no_teo(0) ne 0. then begin

;create the list of indexes corresponding to the excluded lines:
;--------------------------------------------------------------
      out=intarr(n_elements(exclu_obs_wvl_no_teo))
      i=0
      for j=0,n_elements(exclu_obs_wvl_no_teo)-1 do begin
         out(i)=where( (exclu_obs_wvl_no_teo(j) ne obs_wvl) eq 0)
         i=i+1
      endfor

      in=indgen(n_obs)
;create the list of indexes corresponding to the lines we keep for the fit:
;-------------------------------------------------------------------------
      remove,out,in

      n_obs=n_obs-n_elements(exclu_obs_wvl_no_teo)

      obs_wvl=obs_wvl(in)
      obs_id=obs_id(in)
      obs_int=obs_int(in)
      obs_sig=obs_sig(in)
      obs_delta_lambda=obs_delta_lambda(in)

      ch_wvl=ch_wvl(*,in)
      ch_id=ch_id(*,in)
      ch_z=ch_z(*,in)
      ch_ion=ch_ion(*,in)
      ch_contr_wa=ch_contr_wa(*,*,in)
      ch_term=ch_term(*,in)
      ch_n_contr=ch_n_contr(in)

   endif
endelse

;----------------------------------------------------------------------
; NOW WE HAVE THE G(T), PROVIDED BY READING THE FILE OR BY  GET_CONTRIBUTIONS 
;----------------------------------------------------------------------



;NOW WE EXCLUDE THE LINES NOT WANTED FROM THE FOLLOWING:
;(There might be the case that you run without a contribution file, then 
;with a contribution file, without excluding lines)
;---------------------------------------------------------------------

if n_elements (exclude_obs_wvl) gt 0 then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;check the case in which you want to exclude a line that has been already
;automatically excluded because it doesn't have theoretical counterparts:
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if  exclu_obs_wvl_no_teo(0) ne 0. then begin 

      wvl_out=fltarr(n_elements(exclude_obs_wvl))
      coinc=0
      for j=0,n_elements(exclu_obs_wvl_no_teo)-1 do begin
         l=where( (exclude_obs_wvl eq exclu_obs_wvl_no_teo(j)))
         if l(0) ne -1 then begin 
            wvl_out(coinc)=exclude_obs_wvl(l)
            coinc=coinc+1
         endif
      endfor


      if coinc gt 0 then begin


         n=where(wvl_out ne 0.)
         wvl_out=wvl_out(n)

         print,'the line(s) ',wvl_out,' was (were) excluded already !'

;if there are lines that you want to exclude and were not excluded:
;------------------------------------------------------------------

         if n_elements(exclude_obs_wvl)-coinc gt 0 then begin

            for j=0,n_elements(wvl_out)-1 do begin
               pp=where (wvl_out(j) ne  exclude_obs_wvl)
               if pp(0) ne -1 then begin
                  exclude_obs_wvl=exclude_obs_wvl(pp)
               endif
            endfor


;in case you wanted to exclude only lines already excluded, skip.....
;-----------------------------------------------------------------------
         endif else goto,endexclusion

      endif

   endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;end of the case in which you want to exclude a line that has been already
;automatically excluded because it doesn't have theoretical counterparts
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;create the list of indexes corresponding to the excluded lines:
   out=intarr(n_elements(exclude_obs_wvl))
   i=0
   for j=0,n_elements(exclude_obs_wvl)-1 do begin
      out(i)=where( (exclude_obs_wvl(j) ne obs_wvl) eq 0)
      i=i+1
   endfor

   in=indgen(n_obs)
;create the list of indexes corresponding to the lines we keep for the fit:
   remove,out,in

   n_obs=n_obs-n_elements(exclude_obs_wvl)

   obs_wvl=obs_wvl(in)
   obs_id=obs_id(in)
   obs_int=obs_int(in)
   obs_sig=obs_sig(in)
   obs_delta_lambda=obs_delta_lambda(in)

   ch_wvl=ch_wvl(*,in)
   ch_id=ch_id(*,in)
   ch_z=ch_z(*,in)
   ch_ion=ch_ion(*,in)
   ch_contr_wa=ch_contr_wa(*,*,in)
   ch_term=ch_term(*,in)
   ch_n_contr=ch_n_contr(in)

endif
endexclusion:

if n_elements(exclude_obs_wvl) ge 1 then $
   print,'the lines ',exclude_obs_wvl,' will be excluded from the fit !'



;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;                  INSERT THE ABUNDANCE FACTOR  
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_elements(abund_name) gt 0 then begin 

if not file_exist(abund_name) then message,'Error, abundance file does not exist ! '

endif else get_abund_data, abund_name


;read the abundance file 
;------------------------------------
; the abundances are converted from  
; logaritmic values....abund(g)=10.^(abund(g)-12.)
;----------------------------------------------------------
read_abund,abund_name,abund,abund_ref


;;MULTIPLY FOR THE ABUNDANCE FACTOR THE G(T)
;*********************************************

; GDZ
ch_contr=ch_contr_wa            ;fltarr(n_dem_temp,n_ch,n_obs)

for iobs=0,n_obs-1 do begin
   n_list=ch_n_contr(iobs)
   for ilist=0,n_list-1 do begin
      ch_contr(*,ilist,iobs)=abund(ch_z(ilist,iobs)-1)*ch_contr_wa(*,ilist,iobs)
   endfor
endfor

;SORT FOR EACH OBSERVED LINE THE CORRESPONDING THEORETICAL LINES ,BY THEIR G(T)
;-----------------------------------------------------------------------------

;first find out where the maximum of the G(T) is for each line:
;-------------------------------------------------------------
ch_contr_max=fltarr(n_ch,n_obs)
for iobs=0,n_obs-1 do begin
   n_list=ch_n_contr(iobs)
   for ilist=0,n_list-1 do begin
      ch_contr_max(ilist,iobs)=max(ch_contr(*,ilist,iobs))
   endfor
endfor

for iobs=0,n_obs-1 do begin
   ilist_sort=sort(ch_contr_max(0:ch_n_contr(iobs)-1,iobs))

   ch_wvl(0:ch_n_contr(iobs)-1,iobs)=ch_wvl(ilist_sort,iobs)
   ch_id(0:ch_n_contr(iobs)-1,iobs)=ch_id(ilist_sort,iobs)
   ch_z (0:ch_n_contr(iobs)-1,iobs)   =ch_z(ilist_sort,iobs)
   ch_ion (0:ch_n_contr(iobs)-1,iobs) =ch_ion(ilist_sort,iobs)
   ch_contr(*,0:ch_n_contr(iobs)-1,iobs) = ch_contr(*,ilist_sort,iobs)
   ch_term (0:ch_n_contr(iobs)-1,iobs) =ch_term(ilist_sort,iobs)
endfor
;----------------------------------------------------------------------



;----------------------------------------------------------------------------
; Define the main temperature grid (temperatures are logT),
; derived from the grid in the ionization eq file by default
;----------------------------------------------------------------------------

d_dem_temp= IONEQ_LOGT[1]-IONEQ_LOGT[0] ;  
log_dem_temp=IONEQ_LOGT
n_dem_temp=n_elements(log_dem_temp)
dlnt=alog(10.^d_dem_temp)
dem_temp=10.^log_dem_temp

; OPTIONAL INPUTS:
;    max_t  -  maximum temperature (log T) for DEM calculation
;    min_t  -  minimum temperature (log T) for DEM calculation
;    dt     -  temperature step (log T) for DEM calculation
      
   if n_elements(min_logT) eq 0 then min_logT=min(IONEQ_LOGT)
   if n_elements(max_logT) eq 0 then max_logT=max(IONEQ_LOGT)

   
   if min_logT lt min(log_dem_temp)  then $
      message,'requested  MIN_T  outside the G(T) temperature range, Aborting.'
   
   if max_logT gt max(log_dem_temp) then  $
      message,'requested  MAX_T  outside the G(T) temperature range, Aborting.'
    
  if n_elements(dt_logT) ne 0 then begin 

 if dt_logT lt (IONEQ_LOGT[1]-IONEQ_LOGT[0]) then begin 
    print,'Warning: requested bin in log T is smaller than the bin in the IONEQ file:'+$
    trim((IONEQ_LOGT[1]-IONEQ_LOGT[0]))
   end 

   endif else begin 
 print, 'Assuming by default a bin in log T [K] =0.05'
wait,1
dt_logT =0.05
end  



;  plot out all contribution functions
;-------------------------------------

if keyword_set(plot_gt) then begin 


;cmax=max(ch_contr)
;plot_io,fltarr(2),xr=[min_logT,max_logT],yr=[cmax/1.e+3,cmax]

   for iobs=0,n_obs-1 do begin

      n_list=ch_n_contr(iobs)

      print,' '
      print,'Observed line (wavelength , identification) : '
      print,' '
      print,	trim (obs_wvl(iobs))+'  '+STRING(197b)+'  '+obs_id(iobs)
      print,' '
      print,'CHIANTI  lines (wavelength (A), ion,  G(T) max  ) :' 
      print,' '

      y_min=1.e-30
      y_max=1.e-22
      x_min=min_logT
      x_max=max_logT

      device,window_state=ws
      if(ws(11) ne 1) then  window,11,ysize=425,xsize=525  else wset,11

      pr=''
begin_plot_gt:

      plot_io,fltarr(2),xr=[x_min,x_max],xstyle=1,ystyle=1,$
              yr=[y_min,y_max],xtitle = 'log T [ !eo!nK ]',$
              ytitle = ' G(T) [ erg cm!e3!n s!e-1!n sr!e-1!n ]',$
              title='CHIANTI SPECTRAL CODE ',chars=1.2

      for ilist=0,n_list-1 do begin

         if total(ch_contr(*,ilist,iobs)) eq 0. then begin 
            print,'NO  LINES !!!'
         endif else begin

            tmax=where(ch_contr(*,ilist,iobs) eq max(ch_contr(*,ilist,iobs)))
            print,ch_wvl(ilist,iobs),' ',ch_id(ilist,iobs),' ',$
                  max(ch_contr(*,ilist,iobs)),' '

            oplot,log_dem_temp, ch_contr(*,ilist,iobs) 

            comment=strtrim(string(ch_id(ilist,iobs)),2)+$
                    string(ch_wvl(ilist,iobs),'(f8.3)')+' '+string("305B) ; "

            if  max(ch_contr(*,ilist,iobs)) gt y_min then $
               xyouts,log_dem_temp(tmax), max(ch_contr(*,ilist,iobs)),$
                      comment,chars=1.2 
         endelse
      endfor


      break_file,abund_name, disk,dir,f,ext
      xyouts,0.500,0.90,/normal, 'Abund.  :  '+f,chars=1.2 ,alignment=0.5
      break_file,ioneq_name, disk,dir,f,ext
      xyouts,0.500,0.87,/normal, 'Ion Eq.  : '+f,chars=1.2 ,alignment=0.5
      xyouts,0.500,0.84,/normal, const_net,      chars=1.2   ,alignment=0.5

      print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                    go_to_line=go_to_line,out_name=out_name, /ask_name  

      if go_to_line eq 'y' then goto,begin_plot_gt



   endfor                       ;iobs

endif                           ;,/plot_gt

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; Sum up the G(T) within each blend, and find where is the peak:
;---------------------------------------------------------------

ch_tot_contr=fltarr(n_dem_temp,n_obs)

for iobs=0,n_obs-1 do begin
   n_list=ch_n_contr(iobs)
   for ilist=0,n_list-1 do begin
      ch_tot_contr(0,iobs)=ch_tot_contr(*,iobs)+ch_contr(*,ilist,iobs)
   endfor
;  print,ch_tot_contr(*,iobs)
endfor




;=== Prepare input Temperatures
   
   ntemp = long(((max_logT-min_logT)/dt_logT)+1)
   logT_interpolated = findgen(ntemp) * dt_logT + min_logT
   
   ch_tot_contr_interpolated=fltarr(ntemp,n_obs)

print, 'Interpolating the G(T) ... '
; interpolate the emissivities onto the log T grid:
   for iobs=0,n_obs-1 do $
      ch_tot_contr_interpolated[0:ntemp-1,iobs]=interpol(ch_tot_contr[*,iobs],log_dem_temp,logT_interpolated)

; avoid zeros:
ch_tot_contr_interpolated=ch_tot_contr_interpolated>0.


;get the maximum and the temperature of the maximum of the summed G(T)
;----------------------------------------------------------------------
ch_tot_contr_max=fltarr(n_obs)
temp_max_tot_contr=fltarr(n_obs)

for iobs=0,n_obs-1 do begin
   ch_tot_contr_max(iobs)=max(ch_tot_contr(*,iobs))
   temp_max_tot_contr(iobs)=min_logT+d_dem_temp*$
                            ( where (ch_tot_contr_max(iobs) eq (ch_tot_contr(*,iobs)) ))
endfor



;  plot out all summed contribution functions
;-------------------------------------
device,window_state=ws
if(ws(13) ne 1) then  window,13,ysize=425,xsize=525  else wset,13

y_min=1.e-27
y_max=1.e-22
x_min=min_logT
x_max=max_logT

pr=''
begin_plot_summed_gt:

plot_io,fltarr(2),xr=[x_min,x_max],xstyle=1,ystyle=1,$
        yr=[y_min,y_max],xtitle = 'log T [ !eo!nK ]',$
        ytitle = 'summed G(T) [ erg cm!e3!n s!e-1!n sr!e-1!n ]',$
        title='CHIANTI SPECTRAL CODE ',chars=1.2

for iobs=0,n_obs-1 do begin
   
   oplot,logT_interpolated, ch_tot_contr_interpolated(*,iobs)

   comment=strtrim(string(obs_id(iobs)),2)
   
   xyouts,temp_max_tot_contr(iobs),$
          ch_tot_contr_max(iobs),$
          comment,chars=1.2  
endfor

break_file,abund_name, disk,dir,f,ext
xyouts,0.600,0.90,/normal, 'Abund.  :  '+f,chars=1.2 ,alignment=0.5
break_file,ioneq_name, disk,dir,f,ext
xyouts,0.600,0.87,/normal, 'Ion Eq.  : '+f,chars=1.2 ,alignment=0.5
xyouts,0.600,0.84,/normal, const_net,      chars=1.2 ,alignment=0.5


print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
              go_to_line=go_to_line,out_name=output+'_gt.ps'   

if go_to_line eq 'y' then goto,begin_plot_summed_gt


;---------------------------------------------------------------
; plot the EM loci curves 

yes_no,'Do you want to plot the EM loci curves ?',yesno

if yesno then begin 
   
   order_cut=4
   print, 'using an  order='+trim(order_cut)+' to cut out in the plot small G(T) values'

   em_lines_gt = fltarr(n_elements(logT_interpolated), n_obs)
   
   
   FOR  iobs=0,n_obs-1 do begin
      
      em_lines_GT(*, iobs) = obs_int(iobs)/ch_tot_contr_interpolated(*,iobs)
      
      domain=where(ch_tot_contr_interpolated(*,iobs) LE  max(ch_tot_contr_interpolated(*,iobs))*10.^(-order_cut) )

      em_lines_GT(domain, iobs) = 0.
      
      
   endfor 

   window, 0

   x_min= 4.                     
   x_max=6.5                     
   y_min = 23.
   y_max = 30.
   title = 'EM loci '
   
   pr = ''
begin_plot_igt:

   plot ,logT_interpolated , fltarr(n_elements(logT_interpolated))  ,PSYM=6, symsize=1.2,$
         xr=[x_min,x_max],yr=[y_min,y_max], /nodata, $
         xstyle=1,ystyle=1, xtitle = 'log T [K]',$
         ytitle ='log EM [cm!E-5!N]', $
         title=title,chars=1.5

   for iobs=0,n_obs-1 do begin

      good = where(em_lines_GT(*, iobs) NE 0, nn)
      
      IF nn GT 0 THEN oplot, logT_interpolated[good], alog10(em_lines_GT[good, iobs])
      
      inn = where(alog10(em_lines_GT[*, iobs]) EQ max(alog10(em_lines_GT[*, iobs])))
      inn = inn(0)

      xyouts, logT_interpolated(inn), alog10(em_lines_GT(inn, iobs))+0.0 ,$ 
              ' '+strtrim(obs_id(iobs),2), $ ;+' '+strtrim(obs_wvl,2) ,$
              chars=1.2, align=0.5           ;     Orientat
      
   endfor 

   print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                 go_to_line=go_to_line,out_name=output+'_i_gt.ps' ;out_name, /ask
   if go_to_line eq 'y' then goto,begin_plot_igt

   !p.thick=1
   
endif                           ; EM loci


;---------------------------------------------------------------

;---------------------------------------------------------------
if keyword_set(do_xrt_dem) then begin 
;---------------------------------------------------------------

   

if  keyword_set (dem_file)  then begin

   dem_name = ch_get_file( path=!xuvtop+'/dem', filter='*.dem',  tit=' Select a DEM file with spline values') 
   
   read_dem,dem_name,log_t_f,log_dem_f,tref

;THIS is required: if you don't have DEM values at the extremes of the 
;range, the spline produces values that then create problems...

; check the T ranges in the DEM file:

if min_logT  lt min(log_t_f) then $
  message,'ERROR ! - minimum temperature is lower than minimum in  DEM file.. EXIT.' 

if max_logT gt max(log_t_f) then $
  message,'ERROR ! - maximum temperature is higher than maximum in  DEM file.. EXIT.' 

; OK, we have an overlapping interval.
;define the mesh points of the spline:
;-------------------------------------

ind=where(log_t_f ge  min_logT and log_t_f le max_logT)
spl_logt=log_t_f[ind]

; add extreme points:
if min_logT lt min(spl_logt) then spl_logt=[min_logT, spl_logt]
if max_logT gt max(spl_logt) then spl_logt=[spl_logt, max_logT]

   spl_logdem=spline(log_t_f,log_dem_f, spl_logT)

   device,window_state=ws
   if(ws(10) ne 1) then  window,10,ysize=425,xsize=525  else wset,10

   x_min=min_logT
   x_max=max_logT
   y_min=float(min(spl_logdem))
   y_max=float(max(spl_logdem))


   break_file,dem_name, disk,dir,f,ext

   plot,   log_t_f,log_dem_f ,xr=[x_min,x_max],yr=[y_min,y_max],$
           psym=1,charsize=0.9,$ ;   xstyle=1,ystyle=1,$
           title='file: '+f,xtit='log T', ytit='log DEM(T)'

   oplot,  spl_logt, spl_logdem


endif  else if  n_elements(spl_logt) gt 0 then begin 

; check the input values:

if min_logT ne  min(spl_logt) then message,'Error ! min(spl_logt) must be = min_logT'
if max_logT ne  max(spl_logt) then message,'Error ! max(spl_logt) must be = max_logT'

; define the spline nodes log DEM values: 
if n_elements(spl_logdem) eq 0 then begin 

   spl_logdem=22.+fltarr(n_elements(spl_logt))

end 
end 
   
   
; define main arrays:
   
   logT_out=logT_interpolated
   log_dem_out= spline(spl_logt , spl_logdem, logT_out)
   dem_out= 10.^log_dem_out
   
   OUT_LOGT=logT_interpolated
   OUT_LOGDEM=log_dem_out

;--------------------------------   
   
; predicted intensity
   exp_int = fltarr(n_obs) 
; effective temperature:
   t_eff= fltarr(n_obs)
   
   for iobs=0,n_obs-1 do BEGIN
; approximate values:
      exp_int[iobs]=total(ch_tot_contr_interpolated[0:ntemp-1,iobs] *(10.^logT_out)*dem_out*alog(10.^dt_logT))
      
      t_eff[iobs]=total(ch_tot_contr_interpolated[0:ntemp-1,iobs]*$
                        10.^(2*logT_out) *dem_out*alog(10.^dt_logT))/ exp_int[iobs]

   endfor 
   
   print, '    wvl     Iobs   log Teff Icalc/Iobs   ID  '

   sort_t=sort(t_eff) 
   
   for iobs=0,n_obs-1 do $
      print, obs_wvl(sort_t[iobs]), obs_int(sort_t[iobs]),  $
             alog10(t_eff[sort_t[iobs]]), $
             exp_int(sort_t[iobs])/obs_int(sort_t[iobs]) ,obs_id(sort_t[iobs]),$
             format='(f10.3,2x,e8.2,2x,f4.2,2x,f5.2,2x,a)'
   window,1
   
   x_min=min(logT_out)
   x_max=max(logT_out)
   y_min=min(alog10(dem_out))
   y_max=max(alog10(dem_out))

   pr = ''
   begin_plot_xrt_dem_i:

   plot,  logT_out, alog10(dem_out),chars=1.4, $
          xr=[x_min,x_max],yr=[y_min,y_max],$
          xstyle=1,xtitle = ' log Teff [ !eo!nK ]',$
          ytitle ='log DEM [ cm!S!E-5 !NK!S!E-1!N ] ',ystyle=1,$
          title='CHIANTI DEM INVERSION '

; oplot the  observed/expected ratio * DEM 
;----------------------------------------
   
   point=fltarr(n_obs)     

   for iobs=0,n_obs-1 do begin         
      point[iobs]=spline(logT_out, alog10(dem_out), alog10(t_eff[iobs]))         
      xyouts, alog10(t_eff[iobs]), alog10(obs_int[iobs]/exp_int[iobs]*$
                                          10.^point[iobs]), $
              ' '+strtrim(obs_id[iobs],2),$
                                ;+' '+$
                                ;strtrim(obs_wvl[iobs],2)+' '+STRING(197b) ,$
              charsize=0.8, Orientation=90
   endfor     

   oplot, alog10(t_eff), alog10(obs_int/exp_int* 10.^point), psym=6

   print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                 go_to_line=go_to_line,out_name=output+'_mpfit_dem_teff.ps'   

   if go_to_line eq 'y' then goto, begin_plot_xrt_dem_i
    
   
;----------------------------------------------------   
   
   yes_no,'Run the MPFIT ? ', yesno_run, 'Y'
   
   if yesno_run then begin 
   
;       SOLV_FACTOR - [Optional] (float scalar)
;                     The least-squares solver is not completely 
;                     insensitive to the order of magnitude of the numbers
;                     it is manipulating. SOLV_FACTOR is used to 
;                     normalize the inputs to move the solver into a 
;                     "sweet spot". The default choice is arbitrary but 
;                     seems to work well (default = 1e21). Of course,
;                     the outputs are un-normalized at the end.
   
   default, solv_factor, 21
   
;       MAXITER     - [Optional] (float scalar)
;                     This program works by iterating a least-squares
;                     search. This keyword may be used to specify the 
;                     maximum number of iterations for each DEM solution.
;
   default, maxiter,     4000


   delvarx, logT_out, dem_out, out_logt, out_logdem, log_dem_out
   
   
   mpfit_dem, obs_int, ch_tot_contr_interpolated, logt_interpolated,$
              obs_err=obs_sig,$              
              dt=dt_logT,$       ;
              min_logt=min_logT, max_logt=max_logT,$
              spl_logt=spl_logt, spl_logdem=spl_logdem, $
              min_limits=min_limits, max_limits=max_limits,$
              out_logt=out_logt, out_logdem=out_logdem, $ ; output                            
              error=error, solv_factor=solv_factor,  maxiter=maxiter ,verbose=verbose
   
   if error eq 1 then begin 
      print,'ERROR in  MPFIT_DEM ! - rerun with different input parameters. '
   goto, this_end_mpfit
   endif 
   
endif 
        
   logT_out = out_logt
   dem_out = 10.^out_logdem
   log_dem_out=out_logdem
   
                                ;Create the output DEM file
;---------------------------

   dem_name=output+'_mpfit.dem'
   openw,lun_dem,dem_name,/get_lun

   for i=0, n_elements(logT_out)-1 do $
      printf,lun_dem,logT_out[i], alog10(dem_out[i])

   printf,lun_dem,'-1'
   printf,lun_dem,'%file:  ',output+'.dem'
   printf,lun_dem,'% DEM:  Produced by CHIANTI_DEM '
   
   printf,lun_dem,  '% DEM obtained with the MPFIT  program, in the log T ='+$
          trim(min_logT)+'-'+trim(max_logT)+' range and step (log T)='+trim(dt_logT)
   
   printf,lun_dem,'% With the ionization equilibrium file ',  ioneq_name
   
   printf,lun_dem,'% With the abundance file ',abund_name
   printf,lun_dem,'% the observation file "',file_input_obs
   printf,lun_dem,'% calculated at '+const_net
   printf,lun_dem,'%  '
   printf,lun_dem,'-1'

   free_lun,lun_dem


   
; predicted intensity
   exp_int = fltarr(n_obs) 
; effective temperature:
   t_eff= fltarr(n_obs)
   
   for iobs=0,n_obs-1 do BEGIN
; approximate values:
      exp_int[iobs]=total(ch_tot_contr_interpolated[0:ntemp-1,iobs] *(10.^logT_out)*dem_out*alog(10.^dt_logT))
      
      t_eff[iobs]=total(ch_tot_contr_interpolated[0:ntemp-1,iobs]*$
                        10.^(2*logT_out) *dem_out*alog(10.^dt_logT))/ exp_int[iobs]

   endfor 
   
   
; save the results:
   save, file=output+'_mpfit_dem.save',/ver,ch_tot_contr_interpolated, logT_interpolated, logT_out,log_dem_out,$
         obs_int,obs_id, obs_wvl, exp_int,t_eff,temp_max_tot_contr, /compress
   
   print, '    wvl     Iobs   log Teff Icalc/Iobs   ID  '

   sort_t=sort(t_eff) 
   
   for iobs=0,n_obs-1 do $
      print, obs_wvl(sort_t[iobs]), obs_int(sort_t[iobs]),  $
             alog10(t_eff[sort_t[iobs]]), $
             exp_int(sort_t[iobs])/obs_int(sort_t[iobs]) ,obs_id(sort_t[iobs]),$
             format='(f10.3,2x,e8.2,2x,f4.2,2x,f5.2,2x,a)'
   
;--------------------------------------------------------------------------------

openw,2, output+'_info.tex'

 for iobs=0,n_obs-1 do begin 

ratio=exp_int(sort_t[iobs])/obs_int(sort_t[iobs])
this_int=obs_int(sort_t[iobs])

   if keyword_set(phot) then this_int=this_int /1.9866e-8*obs_wvl[sort_t[iobs]]
     if keyword_set(arcsec) then  this_int=this_int* (!PI/(180.*60.^2))^2
 
    if this_int ge 1. then this_int_s= trim(string(this_int , format='(f12.1)')) else $
       if this_int lt 1.  and this_int gt 0.05 then this_int_s= trim(string(this_int , format='(f12.2)')) else $
      this_int_s= trim(string(this_int , format='(f13.3)'))
     

  
pp= str_sep(obs_id[sort_t[iobs]], ' ',/trim) 
good=where(trim(pp) ne '')
pp=pp[good]

if n_elements(pp) eq 2 then $
   ion_string=' \ion{'+trim(pp[0])+'}{'+strlowcase(trim(pp[1]))+'}' else $
      ion_string=pp[0]

   text2= ion_string+' & '+string(obs_wvl(sort_t[iobs]), format='(f7.2)')+' & '+$        
           this_int_s+' & '+$
           string(temp_max_tot_contr[sort_t[iobs]] , format='(f5.2)')+' & '+$
           string(alog10(t_eff(sort_t[iobs])), format='(f5.2)')+' & '+$           
           string(ratio, format='(f5.2)')+' &  &  & \\'

printf,2, text2

 endfor 

close,2


;--------------------------------------------------------------------------------


   
   window,1
   
   x_min=min(logT_out)
   x_max=max(logT_out)
   y_min=min(alog10(dem_out))
   y_max=max(alog10(dem_out))

   pr = ''
   begin_plot_xrt_dem:

   plot,  logT_out, alog10(dem_out),chars=1.4, $
          xr=[x_min,x_max],yr=[y_min,y_max],$
          xstyle=1,xtitle = ' log Teff [ !eo!nK ]',$
          ytitle ='log DEM [ cm!S!E-5 !NK!S!E-1!N ] ',ystyle=1,$
          title='CHIANTI DEM INVERSION '

; oplot the  observed/expected ratio * DEM 
;----------------------------------------

   point=fltarr(n_obs)     

   for iobs=0,n_obs-1 do begin         
      point[iobs]=spline(logT_out, alog10(dem_out), alog10(t_eff[iobs]))         
      xyouts, alog10(t_eff[iobs]), alog10(obs_int[iobs]/exp_int[iobs]*$
                                          10.^point[iobs]), $
              ' '+strtrim(obs_id[iobs],2),$
                                ;+' '+$
                                ;strtrim(obs_wvl[iobs],2)+' '+STRING(197b) ,$
              charsize=0.8, Orientation=90
   endfor     

   oplot, alog10(t_eff), alog10(obs_int/exp_int* 10.^point), psym=6

   print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                 go_to_line=go_to_line,out_name=output+'_mpfit_dem_teff.ps'   

   if go_to_line eq 'y' then goto, begin_plot_xrt_dem
   

; over-plot the  observed/expected ratio * DEM at the temperature of the maximum of the G(T):
   window,2
   pr=''      
   begin_plot_xrt_dem2:

   plot,  logT_out, log_dem_out, chars=1.4, $
          xr=[x_min,x_max],yr=[y_min,y_max],$
          xstyle=1,xtitle = ' log Tmax [ !eo!nK ]',$
          ytitle ='log DEM [ cm!S!E-5 !NK!S!E-1!N ] ',ystyle=1,$
          title='CHIANTI DEM MPFIT'
   

   for iobs=0,n_obs-1 do begin        
      point=spline(logT_out, log_dem_out, temp_max_tot_contr[iobs])
      oplot, [temp_max_tot_contr[iobs]], [alog10(obs_int[iobs]/exp_int[iobs]* 10.^point)], psym=6
      xyouts, temp_max_tot_contr[iobs], alog10(obs_int[iobs]/exp_int[iobs]*$
                                               10.^point), $
              ' '+strtrim(obs_id[iobs],2), charsize=0.8, Orientation=90
   endfor


   print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                 go_to_line=go_to_line,out_name=output+'_mpfit_dem_tmax.ps'   

   if go_to_line eq 'y' then goto,begin_plot_xrt_dem2

;Open the output file where all the info will be stored:
;-------------------------------------------------------
   output_general=output+'.general'
   openw,luo,output_general,/get_lun
;
   printf,luo,' abundance file = '+abund_name
   printf,luo,' ionization equilibrium file = '+ioneq_name
   printf,luo,const_net

   IF  keyword_set(verbose) THEN BEGIN
      print,' '
      print,'------------------------------------------------ '
      print,'      not  sorted '
      print,'------------------------------------------------ '
      print,' '
      print,' Row     ident. obs.lambda  obs.flux  calc.flux   I_calc/I_obs'
      print,'       Ion  wavelength  terms  fractional contribution to  I_calc '

   ENDIF
;
;
   printf,luo,' '
   printf,luo,'------------------------------------------------ '
   printf,luo,'     not sorted  '
   printf,luo,'------------------------------------------------ '
   printf,luo,' '

   printf,luo,' Row     ident. obs.lambda  obs.flux  calc.flux   I_calc/I_obs'
   printf,luo,'  Ion  wavelength  terms  fractional contribution to  I_calc '

;
   for iobs=0,n_obs-1 do begin

      IF  keyword_set(verbose) THEN BEGIN
         print, ' '
         print, ' '
         print,iobs,obs_id[iobs],obs_wvl[iobs],obs_int[iobs],exp_int[iobs],$
;               obs_sig[iobs],(exp_int[iobs]-obs_int[iobs])^2/obs_sig[iobs]^2,$
               exp_int[iobs]/obs_int[iobs], $
               format='(i5,1x,a20,1f10.3,1x,2e10.3,f5.2)'
      ENDIF

      printf,luo, ' '
      printf,luo, ' '
      printf,luo,iobs,obs_id[iobs],obs_wvl[iobs],obs_int[iobs],exp_int[iobs],$
;             obs_sig[iobs],(exp_int[iobs]-obs_int[iobs])^2/obs_sig[iobs]^2, $
             exp_int[iobs]/obs_int[iobs], $
             format='(i5,1x,a20,1f10.3,1x,2e10.3,f5.2)'
;
      n_list=ch_n_contr[iobs]

      for ilist=0,n_list-1 do begin
         
         this_contr=interpol(ch_contr[*,ilist,iobs],IONEQ_LOGT,logT_interpolated)
         this_contr = total(this_contr*(10.^logT_out)*dem_out*alog(10.^dt_logT))
         this_contr=this_contr/exp_int[iobs]

         if this_contr gt 0.01 then begin 
            IF  keyword_set(verbose) THEN BEGIN

               if ilist eq 0 then $
                  print, '------------------------------------------------------------------------------------------------'
               
               print,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
                     ch_term(ilist,iobs),this_contr,$
                     format='(10x,a15,f10.3,a40,f10.2)'
;
            ENDIF  

            if ilist eq 0 then $
               printf,luo, '------------------------------------------------------------------------------------------------'

            printf,luo,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
                   ch_term(ilist,iobs),this_contr,$
                   format='(10x,a15,f10.3,a40,f10.2)'
         endif 

      endfor 

      printf,luo,' '

;
   endfor                       ; iobs  



   IF  keyword_set(verbose) THEN BEGIN

      print,'-------------------------------------- '
      print,' sorted by observed wavelength'
      print,'-------------------------------------- '
      print,' '
      print,' Row     ident. obs.lambda  obs.flux  calc.flux     I_calc/I_obs'
      print,' '

   ENDIF

;
   printf,luo,' '
   printf,luo,'-------------------------------------- '
   printf,luo,' sorted by observed wavelength'
   printf,luo,'-------------------------------------- '
   printf,luo,' '
   printf,luo,' Row     ident. obs.lambda  obs.flux  calc.flux   I_calc/I_obs'
   printf,luo,' '
;
   measure=obs_wvl
;
;
   isort=sort(measure)
;
   for j=0,n_obs-1 do begin

      iobs=isort(j)

      IF  keyword_set(verbose) THEN BEGIN

         print, ' '
         print, ' '
         print,iobs,obs_id[iobs],obs_wvl[iobs],obs_int[iobs],exp_int[iobs],$
               exp_int[iobs]/obs_int[iobs], $
               format='(i5,1x,a20,1f10.3,1x,2e10.3,1f10.2)'

      ENDIF
;
      printf,luo, ' '
      printf,luo, ' '
      printf,luo,iobs,obs_id[iobs],obs_wvl[iobs],obs_int[iobs],exp_int[iobs],$
             exp_int[iobs]/obs_int[iobs], $
             format='(i5,1x,a20,1f10.3,1x,2e10.3,1f10.2)'
;
      n_list=ch_n_contr[iobs]

      for ilist=0,n_list-1 do begin
         
                  this_contr=interpol(ch_contr[*,ilist,iobs],IONEQ_LOGT,logT_interpolated)
         this_contr = total(this_contr*(10.^logT_out)*dem_out*alog(10.^dt_logT))
         this_contr=this_contr/exp_int[iobs]

         if this_contr gt 0.01 then begin 
            IF  keyword_set(verbose) THEN BEGIN

               if ilist eq 0 then $
                  print, '------------------------------------------------------------------------------------------------'

               print,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
                     ch_term(ilist,iobs),this_contr,$
                     format='(10x,a15,f10.3,a40,f10.2)'
;
            ENDIF

            if ilist eq 0 then $
               printf,luo, '------------------------------------------------------------------------------------------------'

            printf,luo,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
                   ch_term(ilist,iobs),this_contr,$
                   format='(10x,a15,f10.3,a40,f10.2)'
;
         endif 

      endfor
;
   endfor
;

   printf,luo,' '
   printf,luo,' '

   if exclu_obs_wvl_no_teo(0) ne 0. then begin
      print,'The following lines have been excluded  from the fit'+ $
            ' because they DO NOT have theoretical counterparts :'
      printf,luo,'The following lines have been excluded from the fit'+$
             ' because they DO NOT have theoretical counterparts :'
      print,exclu_obs_wvl_no_teo
      printf,luo,exclu_obs_wvl_no_teo
   endif

   printf,luo,' '
   printf,luo,' '

   if n_elements(exclude_obs_wvl) gt 0 then begin	
      print,'The following lines have been excluded  from the fit:'
      printf,luo,'The following lines have been excluded  from the fit:'
      print,exclude_obs_wvl
      printf,luo,exclude_obs_wvl
   endif


   free_lun,luo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;write the output with the useful information for user-written plotting 
;routines.

   openw,lu,output+'.out',/get_lun

   for iobs=0,n_obs-1 do begin

      printf,lu, obs_id[iobs], obs_wvl[iobs], obs_int[iobs], exp_int[iobs],$
             obs_sig[iobs], $
             format='(a20,1x, 1f10.3,1x, 3e10.3)'

   endfor

   free_lun,lu

 
   
   !p.thick=1

;=== Run the Monte Carlo loops ===========================

   yes_no, 'Run the Monte Carlo simulations ?', yesno 

   if yesno then begin 
      
      MC_iter=400
      read,'type the number of iterations (e.g. 400):',MC_iter
      
      log_dem_mciter=fltarr(n_elements(logT_out), MC_iter)
      log_dem_mciter[*,0]=log_dem_out

      seed = systime(1)

      for ii = 1,MC_iter-1 do begin
         print,' M.C. loop ' +              $
               strcompress(ii,/rem) + ' of ' + strcompress(MC_iter,/rem) + '.'
         
         mpfit_dem, (obs_int+ randomn(seed, n_obs)*obs_sig) > 0.0, $
                    ch_tot_contr_interpolated, logt_interpolated,$
                    obs_err=obs_sig,$              
                    dt=dt_logT,$ ;
                    min_logt=min_logT, max_logt=max_logT,$
                    spl_logt=spl_logt, spl_logdem=spl_logdem, $
                    min_limits=min_limits, max_limits=max_limits,$
                    out_logt=out_logt, out_logdem=out_logdem, $ ; output
                    error=error, solv_factor=solv_factor,  maxiter=maxiter ,verbose=0
         
         if not error then $
         log_dem_mciter[*,ii] = out_logdem
         
      endfor
      
      save, file=output+'_xrt_dem.save',/ver, logT_out, log_dem_mciter,$
            obs_int, obs_sig, obs_id, obs_wvl, exp_int,t_eff,temp_max_tot_contr,/compress
      
      window,3
      pr=''      
      begin_plot_xrt_dem3:

      
      plot, logT_out,log_dem_out, psym=10,th=th, col=0, chars=1.4, $
            xr=[x_min,x_max],yr=[y_min,y_max],$
            xstyle=1,xtitle = ' log T [ !eo!nK ]',$
            ytitle ='log DEM [ cm!S!E-5 !NK!S!E-1!N ] ',ystyle=1,$
            title='XRT DEM INVERSION TECHNIQUE'

; over plot the other solutions:
      for ii=1, MC_iter-1 do oplot, logT_out, log_dem_mciter[*,ii], th=th, col=100, psym=10
      oplot, logT_out,log_dem_out, psym=10,th=th, col=0
      
      
      print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                    go_to_line=go_to_line,out_name=output+'_xrt_dem_mc.ps'   

      if go_to_line eq 'y' then goto,begin_plot_xrt_dem3
      !p.thick=1       
      
   endif                        ; run the  Monte Carlo loops 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


this_end_mpfit:   
   
endif                           ; xrt_dem

;--------------------------------------------------------------------------------

if keyword_set(do_demreg) then begin 
   
   
   
; ;order of regularization, default is 0th
   if n_elements( order_demreg) eq 0 then  order_demreg=0
; ;control the regularization parameter/chisq of result in DEM space: reg_tweak=1
   if n_elements( reg_tweak_demreg) eq 0 then reg_tweak_demreg=1
; ;Use guess solution in final regularization? default is no, guess=0.
   if n_elements( guess_demreg) eq 0 then guess_demreg=0.
;; Use the min of the EM loci curves as the initial guess solution
;; used to weight/create the constraint matrix and possibly in the regularization itself (if guess=1)
   if n_elements( gloci_demreg) eq 0 then gloci_demreg=0


   if n_elements(demreg_logt_min) eq 0 then demreg_logt_min=min(logT_interpolated) 
   if n_elements(demreg_logt_max) eq 0 then demreg_logt_max=max(logT_interpolated) 

   if n_elements(nt_demreg) eq 0 then  nt_demreg=20 
; nt must be bigger than nf.
   nt_demreg=nt_demreg > n_elements(obs_int)+1

   reg=run_data2dem_reg(logT_interpolated, ch_tot_contr_interpolated, obs_int, obs_sig,$
                        mint=demreg_logt_min,maxt=demreg_logt_max,nt=nt_demreg,$
                        order=orderi_demreg , guess=guess_demreg ,reg_tweak=reg_tweak_demreg,$
                        channels=channels,debug=debug,gloci=gloci_demreg, /pos, error=error)
   
   if error ne 1 then begin 
                                ; now to plot the results....

; save the output !

      save, file=output+'_demreg.save',reg
      print, 'IDL save file '+output+'_demreg.save  saved !'


      openw,lun_dem,output+'_demreg.dem',/get_lun

      for ii=0,  n_elements(reg.logt )-1 do $
         printf, lun_dem, reg.logt[ii], alog10(reg.dem_pos[ii]) 


      printf,lun_dem,'-1'
      printf,lun_dem,'% DEM:   Produced by CHIANTI_DEM '  
      printf,lun_dem,  '% DEM obtained with the DATA2DEM_REG  program, in the log T ='+$
             trim(demreg_logt_min)+'-'+trim(demreg_logt_max)+' range and number of temperatures='+trim(nt_demreg)
      printf,lun_dem,'% With the ionization equilibrium file ',  ioneq_name  
      printf,lun_dem,'% With the abundance file ',abund_name
      printf,lun_dem,'% the observation file "',file_input_obs
      printf,lun_dem,'% calculated at '+const_net
      printf,lun_dem,'%  '
      printf,lun_dem,'-1'

      free_lun,lun_dem
      
      x_min=min(reg.logt)
      x_max=max(reg.logt)
      y_min=2d19
      y_max=2d23
      
      pr=''      
      begin_plot_demreg:
      
      linecolors
      !p.charsize=1.5
; plot the regularized DEM and both vertical and horizontal errors
      window,1,xsize=650,ysize=500,title='Regularized DEM'
      !p.multi=0
      ploterr,reg.logt,reg.dem_pos,reg.elogt_pos,reg.edem_pos,$
              /nohat,errcolor=9,  xr=[x_min,x_max],yr=[y_min,y_max],$
              xstyle=17,ystyle=17,/ylog,title='Regularized DEM', $
              xtitle='log!D10!N T',ytitle='DEM(T) [cm!U-5!N K!U-1!N]'
      
      
      print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                    go_to_line=go_to_line,out_name=output+'_demreg.ps'   

      if go_to_line eq 'y' then goto,begin_plot_demreg
      !p.thick=1

      window,2,xsize=500,ysize=500,title='Line Intensities and Residuals'     
      !p.multi=[0,1,2]
      loadct,0,/silent
      linecolors
      nf=n_elements(reg.channels)
      plot,indgen(nf),reg.data,/ylog,psym=6,$
           xrange=[-1,nf],xtickf='(a1)',xticks=nf+1,ystyle=16,thick=3,$
           ytit='Line Intensities',xtit=' ',/nodata,$
           yrange=[0.9*min(reg.data),1.1*max(reg.data)]
      oplot,indgen(nf),reg.data_reg,psym=6,color=5,thick=1
      oplot,indgen(nf),reg.data,psym=7,color=2,thick=2
      for i=0, nf-1 do oplot, [i,i],reg.data[i]+[-reg.edata[i],reg.edata[i]],thick=5,color=2


      maxr=1.1*max(abs(reg.residuals))
      plot,indgen(nf),reg.residuals,xrange=[-1,nf],xtickn=[' ',reg.channels,' '],$
           xticks=nf+1,ystyle=17,thick=1,yrange=maxr*[-1,1],psym=6,$
           ytit='Residuals',xtit='Line'
      oplot,[-2,nf],[0,0],lines=1
      xyouts,-0.5,.75*maxr,'chisq='+string(reg.chisq,format='(f4.1)'),/data
      
      
      !p.charsize=1.
      !p.multi=0.

   endif                        ; we have the DATA2DEM_REG  routine.
endif                           ; DEMREG


;-----------------------------------------------------------                          

if keyword_set(do_mcmc) then begin 
   
      
      ;;   convert from  ergs to photons '
      obs_int_phot=obs_int/1.9866e-8*obs_wvl
      obs_sig_phot=obs_sig/1.9866e-8*obs_wvl
      
      demrng=dblarr(nt,2) & demrng[*,0]=1e21 & demrng[*,1]=1e28
            
      dem_out=run_mcmc_dem(obs_wvl,obs_int_phot,1d23*ch_tot_contr_interpolated,$
                           Z=1.+intarr(n_elements(obs_int)),$
                           logt=logt_interpolated, demrng=demrng,$
                           fsigma=obs_sig_phot,ulim=ulim,$
                           diffem=diffem, nsim=nsim,nbatch=nbatch,nburn=nburn,smoot=smoot,$
                           storpar=storpar,storidx=storidx,$
                           simprb=simprb,simdem=simdem,demerr=demerr,simflx=simflx,$
                           simprd=simprd,nosrch=nosrch,softlim=softlim,$
                           sampenv=sampenv,smooscl=smooscl, _extra=_extra, $      
                           savfil=output+'_mcmc.save')
      
      
      if dem_out[0] eq -1 then return 
      
      
; the output of mcmc_dem is a DEM in dln T
      
      dem_out=dem_out/10.d^logt_interpolated
      
;  expected intensity
      exp_int = fltarr(n_obs)
      
; effective temperature:
      t_eff= fltarr(n_obs)
      
      for iobs=0,n_obs-1 do BEGIN
         
; This is  approximate but is what lineflx.pro, called by mcmc_dem,
; does to obtain the predicted intensities
         exp_int[iobs]=total(ch_tot_contr_interpolated[*,iobs]*(10.^logt_interpolated)*dem_out*alog(10.^step))
         
         
         t_eff[iobs]=total(ch_tot_contr_interpolated[*,iobs]*$
                           10.^(2*logt_interpolated) *dem_out*alog(10.^step))/exp_int[iobs]
         
      endfor 
      
; save the above quantities for later on:
      
      save, file=output+'_mcmc2.save', $
            logt_interpolated, dem_out, obs_int,obs_id, obs_wvl, exp_int,t_eff,temp_max_tot_contr
      
      
      print, '    wvl    Iobs   log Teff Iobs/Icalc   ID  '

      sort_t=sort(t_eff) 
      
      for iobs=0,n_obs-1 do $
         print, obs_wvl(sort_t[iobs]), obs_int(sort_t[iobs]),  $
                alog10(t_eff[sort_t[iobs]]), $
                obs_int(sort_t[iobs])/ exp_int(sort_t[iobs]) ,obs_id(sort_t[iobs]),$
                format='(f10.3,2x,e8.2,2x,f4.2,2x,f5.2,2x,a)'

;---------------------------------------------------------------
      
      
      openw,lun_dem,output+'_mcmc.dem',/get_lun

      for ii=0,  n_elements(logT_interpolated)-1 do $
         printf, lun_dem, logt_interpolated[ii],alog10(dem_out[ii])

      printf,lun_dem,'-1'
      printf,lun_dem,'% DEM:   Produced by CHIANTI_DEM '  
      printf,lun_dem,  '% DEM obtained with the MCMC_DEM program, part of PINTofALE, in the log T ='+$
             trim(min_logT)+'-'+trim(max_logT)+' range and step (log T)='+trim(dt_logT)  
      printf,lun_dem,'% With the ionization equilibrium file ',  ioneq_name  
      printf,lun_dem,'% With the abundance file ',abund_name
      printf,lun_dem,'% the observation file "',file_input_obs
      printf,lun_dem,'% calculated at '+const_net
      printf,lun_dem,'%  '
      printf,lun_dem,'-1'

      free_lun,lun_dem 
      
;--------------------------------------------------------------------
      
      window,4
      
      x_min=min(logt_interpolated)
      x_max=max(logt_interpolated)
      y_min=min(alog10(dem_out))
      y_max=max(alog10(dem_out))

      pr = ''
      begin_plot_dem_mcmc:

      plot,  logt_interpolated, alog10(dem_out),chars=1.4, $
             xr=[x_min,x_max],yr=[y_min,y_max],$
             xstyle=1,xtitle = ' log Teff [ !eo!nK ]',$
             ytitle ='log DEM [ cm!S!E-5 !NK!S!E-1!N ] ',ystyle=1, $
             title='MCMC DEM INVERSION TECHNIQUE', psym=10


;oplot the  observed/expected ratio * DEM 
;----------------------------------------

      point=fltarr(n_obs)         

      for iobs=0,n_obs-1 do begin            
         point[iobs]=spline(logt_interpolated, alog10(dem_out), alog10(t_eff[iobs]))            
         xyouts, alog10(t_eff[iobs]), alog10(obs_int[iobs]/exp_int[iobs]*$
                                             10.^point[iobs]), $
                 ' '+strtrim(obs_id[iobs],2),$
                                ;+' '+$
                                ;strtrim(obs_wvl[iobs],2)+' '+STRING(197b) ,$
                 charsize=0.8, Orientation=90

      endfor     

      oplot, alog10(t_eff), alog10(obs_int/exp_int* 10.^point), psym=6
      
      
      print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
                    go_to_line=go_to_line,out_name=output+'_mcmc_teff.ps'   
      
      if go_to_line eq 'y' then goto,begin_plot_dem_mcmc
      !p.thick=1     
      
;      endelse                    ; we have the routine
      
;  mcmc_plot,logt,simdem,demerr,simprb,'PROB',storidx,sampct=nsim/25,$
;        slect=1,col_tabl=1,subtitle='(!4v!X!u2!n='+strtrim(2*min(simprb),2)+'/'+strtrim(n_elements(flx)-total(ulim),2)+')',$
;        ps_fil=output+'_mcmc.ps'
      
      
      
   endif                        ; mcmc


return
end

