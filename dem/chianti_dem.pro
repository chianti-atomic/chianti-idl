;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;                   
; Name        : CHIANTI_DEM
;     		          
; Purpose     : Calculates the Differential Emission Measure DEM(T) using 
;		the CHIANTI database, from a given set of observed lines.
;		Constant pressure or density can be used.
;
; Category    : diagnostic analysis
;               
; Explanation : This routine has several options, all in the form of keywords.
; 		First, the input file with the observed fluxes is read.
;		THE FIRST TIME YOU USE THIS ROUTINE
;		you'll have to do the calculation of the contribution
;		functions G(T), so the routine GET_CONTRIBUTIONS will come 
;		into play. You'll have to specify the value of the 
;		pressure or density, and  you'll be asked to select an  
;		ionization equilibrium file and an abundance file.
;		GET_CONTRIBUTIONS searches the CHIANTI database (ion per ion) 
;		for all the theoretical lines corresponding to the observed 
;		lines, i.e. that lie in a OBS_WVL(i) +/- DELTA_LAMBDA_OBS(i) 
;		interval centered on the observed wavelength OBS_WVL(i).
;		The routine calculates the C(T) values (G(T)=Ab(element)*C(T))
;		for the temperature interval log(T)= 4.0 - 8.0  
;		with steps of log(T) = 0.1 .
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
;		which can take long time.
;
;		In the case no theoretical lines corresponding to an observed
;		line are found, the routine writes the wavelength of the line
;		(to be excluded from the fit) in the array
;		EXCLU_OBS_WVL_NO_TEO. The lines with no theoretical 
;		counterparts are then automatically  excluded from the fit by 
;	 	CHIANTI_DEM. You might consider the possibility to start again
;		incrementing the DELTA_LAMBDA_OBS, to see if there are 
;		theoretical lines in the vicinity.
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
;		line by the abundance factor. Then the
;		theoretical lines contributing to each blend are sorted by
;		intensity and then their G(T) can be plotted if the keyword
;		PLOT_GT is activated. It is recommended to do this the first 
;		time, to check if there are some observed lines terribly 
;		blended with lines of other elements, in which case it is
;		better to exclude them with a second run (if you are not 
;		sure about the abundances).
;		Then the G(T) for each blend are summed and plotted.
;		Then  the fit starts calling DEM_FIT. 
;		A series of parameters can change the 
;		result (DEM), especially the number and position of the mesh
;		points of the spline that represents the DEM. The keyword
;		MESH_POINTS serves for this purpose. 
;		The other keywords that control the fit are N_ITER, DCHISQ_M.
;   		At the end of the fit, the files OUTPUT.DEM and OUTPUT.GENERAL
;		are created.
;
;
;
; Use         : IDL>chianti_dem,output='serts89',file_input='serts89.obs',$
;				pressure=3.e15
;
;
; Examples    : 
;		Assume you have a file input 'serts89.obs' like this:
;		
;		243.031   491.    97.    0.1  He II
;		256.323   1580.   186.   0.1  He II b
;		315.024   253.    31.    0.1  Mg VIII
;		335.401   10400.  1650.  0.1  Fe XVI
;		319.839   113.    14.    0.1  Si VIII
;		356.027   218.    25.    0.1  Si X 
;		
;		IDL>chianti_dem,output='serts89',file_input='serts89.obs',$
;		   pressure=3.e15,cut_gt=1e-30,/plot_gt
;
;		After having selected the  ionization file,
;		the C(T) (with MAX(C(T)) gt 1e-30)  are stored in the file
;		'serts89.contributions'. Then select one of the abundance 
;		files. 
;		Have a look at the plots of the  G(T), and annotate
;		if there is a line you want to exclude, let's say the second.
;		Have a look at the DEM obtained ('serts89.dem') and at 
;		the details contained in the file 'serts89.general'. 
;		Maybe there is another line you want to exclude, let's say 
;		the last one. Maybe you want to change the mesh points, too.
;		So run
;		IDL>chianti_dem,output='serts89_2',file_input='serts89.obs',$
;		file_gt='serts89.contributions',$
;		exclude_obs_wvl=[243.031,356.027 ],$
;		mesh_points= [4.5,5.,5.5,6.2,7.5]
;		
;		The files 'serts89_2.dem' and 'serts89_2.general' will be
;		created. They have the essential information about what you 
;		did.
;
;    
; Inputs      : various, all in form of keywords. The required ones are 
;		OUTPUT and FILE_GT (or  PRESSURE/DENSITY)
;               
;               
; Opt. Inputs : various... see the software note.
;               
; Outputs     : OUTPUT.CONTRIBUTIONS  
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
;		as input for the DMM_SS procedure,that calculates the 
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
;		OUTPUT.OUT
;		This file , toghether with OUTPUT.DEM , 
;		can be used to reproduce the results  using 
;		user-written software. See the software notes.
;		The ouput has this format: 
;		format='(a20,1x, 1f10.3,1x, 3e10.3, 1x,  f4.2,1x,f6.3)'	
;	
; Opt. Outputs:
;		An abundance file with the modifications inserted.
;	
;		Postscript files of the G(T).
;	
;		A postscript file with the DEM (OUTPUT.DEM.PS)
;		
;               A postscript file with other plots too (OUTPUT_4PLOTS.PS)
;
; Keywords    : 
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
;	DCHISQ_M:
;		optional. If not set, a default value of DCHISQ_MIN=1.e-5 
;		is assumed. For each iteration, the CHISQ and it's variation 
;		are calculated. As long as the iteration achieves an
;		improvement in CHISQ greater than  DCHISQ_MIN , another 
;		iteration will be performed.
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
;	MESH_POINTS:    	
;		optional. It is a vector that specifies the mesh points for the
;		spline that represent the fitted DEM, in log(T).
;		If not set, the default values 
;		[4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.0] are assumed.
;
;	N_ITER:
;		optional.It is the number of iterations of the fitting routine.
;		If not set, a default value of 20 is assumed. 
;		Changing this value alone might not affect the fit, since 
;		also the value of DCHISQ_MIN is checked during the fit.
;
;	N_MATCHES:   
;		optional.          
;		In the unlikely event that more than 50 (default value for 
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
;	QUIET:
;		optional. Set to avoid various messages and the details of the 
;		result.
;
; Calls       : GET_CONTRIBUTIONS
;		DEM_FIT
;		ZION2SPECTROSCOPIC
;		print2d_plot
;		
; Common      : obs, 	obs_int,obs_sig,n_obs
;		obs_o,	obs_wvl,obs_id,obs_delta_lambda
; 		dem, 	d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
;		contr,	ch_tot_contr
;		ab,	abund_name,abund_info,xuvtop,ioneq_name
;
;		these are the commons with GET_CONTRIBUTIONS.PRO:
;
;		various,	exclu_obs_wvl_no_teo,const_net,$
;		 dem_temp_min,dem_temp_max,n_dem_temp,$
;		 ch_wvl,ch_l1,ch_l2,ch_id,ch_z,ch_ion,ch_contr_wa,$
;		 ch_pop,ch_contr_list,	 ch_term,ch_n_contr
;
; Restrictions: 
;		In the unlikely event that more than 50 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;
;		Also, if the starting DEM values are not proper, or you 
;		don't have enough constraints at lower and higher temperatures,
;		you might get "strange" results, and should consider using 
;		different starting values.
;
;		Of course you need to have the enviroment variable CDS_SS_DERE
;		pointing to the CHIANTI database top directory.
;
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
; Version     : 2.0 GDZ, DAMTP,  31-Oct-2000
;
;              V.3, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
; VERSION     :  3, 21-May-2002, GDZ 
;-

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;The following is the event processor for the widget that allows you to select
;-----------------------------------------------------------------------------
;any *.abund file present in the CHIANTI database AND in the working directory
;-----------------------------------------------------------------------------

pro get_ab_event,ev

common ab,	abund_name,abund_info,xuvtop,ioneq_name
;
widget_control,ev.id,get_uvalue=uvalue
name = strmid(tag_names(ev,/str),7,1000)
if n_elements(uvalue) eq 0 then uvalue = ''
if n_elements(value)  eq 0 then  value = ''
if uvalue eq 'D' then widget_control ,/destroy,ev.top 

if uvalue eq 'ABUND' then begin

          if tag_exist(ev,'VALUE') then begin
             if ev.value eq '-----' then return
             abund_name = $
 concat_dir(concat_dir(xuvtop, 'abundance'), ev.value+'.abund')
             if not file_exist(abund_name) then abund_name = ev.value+'.abund'
             print,'Using abundance file: ',abund_name
             widget_control, abund_info, set_val=ev.value
          endif else begin
             print,'No abundance file name'
          endelse
endif
end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;                              MAIN PROCEDURE:
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;--------------

pro chianti_dem,output=output,file_input=file_input,pressure=pressure,$
density=density,cut_gt=cut_gt,plot_gt=plot_gt,mesh_points=mesh_points,$
n_matches=n_matches,file_gt=file_gt,n_iter=n_iter,dchisq_m=dchisq_m,$
exclude_obs_wvl=exclude_obs_wvl,dem_file=dem_file,arcsec=arcsec,phot=phot,$
quiet=quiet

;-------------

common obs, 	obs_int,obs_sig,n_obs
common dem, 	d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
common contr,	ch_tot_contr

;these are the commons with GET_CONTRIBUTIONS.PRO:

common obs_o,	obs_wvl,obs_id,obs_delta_lambda
common ab,	abund_name,abund_info,xuvtop,ioneq_name
common various,	exclu_obs_wvl_no_teo,const_net,$
		dem_temp_min,dem_temp_max,n_dem_temp,$
		ch_wvl,ch_l1,ch_l2,ch_id,ch_z,ch_ion,ch_contr_wa,$
		ch_pop,ch_contr_list,	ch_term,ch_n_contr

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	ON_ERROR, 2


!p.background=255
!p.color=0
!p.charsize=1
!p.multi=0


;some boring keyword checking......
;-----------------------------------

if n_elements(output) eq 0  then begin
print,' '
print,'You have to define an output core name (e.g. output="test") !!'
print,' '
return
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

;------------------------------------------------------------
;define the CHIANTI top directory:
;--------------------------------
defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
  message, 'system variable !xuvtop must be set '
xuvtop =!xuvtop


;pick up the file with the fluxes:
;---------------------------------
if  keyword_set(file_input) then file_input_obs=file_input else begin
file_input_obs=pickfile(title='Select the  data file with the fluxes :')
endelse


; expected maximum number of wavelength matches:
;-----------------------------------------------
if  keyword_set(n_matches) then n_ch=n_matches else n_ch=50
;if n_elements(n_matches) eq 0 then n_ch=20  else n_ch=n_matches 


if n_elements(n_iter) eq 0 then niter=20 else niter=n_iter
if  keyword_set(dchisq_m)  then dchisq_min=dchisq_m  else dchisq_min=1.e-5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;read the input file with the observed intensities:
;--------------------------------------------------

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

IF NOT keyword_set(quiet) THEN BEGIN 
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

;Various definitions (temperatures are logT) in common with GET_CONTRIBUTIONS
;----------------------------------------------------------------------------
!p.multi=0
dem_temp_min=4.0
dem_temp_max=8.0
d_dem_temp=0.1
;***************

dlnt=alog(10.^d_dem_temp)
;*************************

n_dem_temp=fix((dem_temp_max-dem_temp_min)/d_dem_temp)+1
;************************************

log_dem_temp=dem_temp_min+findgen(n_dem_temp)*d_dem_temp
dem_temp=10.^log_dem_temp



step=0.01
;*******

n_log_t_dem_step=fix( (dem_temp_max-dem_temp_min)/step) +1    ;  (400/1)

log_t_dem_step=dem_temp_min+findgen(n_log_t_dem_step)* step    ;log_dem_temp con passo step
t_dem_step=10.^log_t_dem_step


;define the mesh points:
;-----------------------
if keyword_set(mesh_points) then begin
 log_t_mesh=mesh_points 
endif else log_t_mesh=[4.0,4.5,5.,5.5,6.,6.5,7.,7.5,8.0]

n_dem_mesh=n_elements(log_t_mesh) 


if  keyword_set (dem_file)  then begin


	ww=''
	print, ' Do you want a CHIANTI DEM (*.dem) standard file [c] '
	read,' or one (*.dem) in your working directory [w] ? ',ww

	if ww ne 'w' and ww ne 'c' then return

	if ww eq 'c' then  dem_name =$
                  pickfile(path=concat_dir(!xuvtop, 'dem'),filter='*.dem',$
        	title='Select DEM File')

	if ww eq 'w' then begin
	cd, current=dir
 	dem_name =pickfile(path=dir,filter='*.dem',$
        	title='Select DEM File')
	endif


        read_dem,dem_name,log_t_f,log_dem_f,tref


;THIS is required: if you don't have DEM values at the extremes of the 
;range, the spline produces values that then create problems...

        log_dem=spline(log_t_f,log_dem_f,log_dem_temp)

	missing_high=where(log_dem_temp gt max(log_t_f), nc )

;replace with the last value....
;--------------------------------
	if nc gt 0 then log_dem(missing_high)= log_dem_f(n_elements(log_t_f)-1)

	missing_low=where( log_dem_temp lt min(log_t_f),nc)

;replace with the first value....
;--------------------------------
	if nc gt 0 then log_dem(missing_low)= log_dem_f(0)

        dem=10.^log_dem

	device,window_state=ws
	if(ws(10) ne 1) then  window,10,ysize=425,xsize=525  else wset,10

        x_min=dem_temp_min
        x_max=dem_temp_max
        y_min=float(min(log_dem))
        y_max=float(max(log_dem))


	break_file,dem_name, disk,dir,f,ext

        plot,   log_t_f,log_dem_f ,xr=[x_min,x_max],yr=[y_min,y_max],$
                psym=1,charsize=0.9,$   ;   xstyle=1,ystyle=1,$
                title='file: '+f,xtit='log T', ytit='log DEM(T)'
        oplot,  log_dem_temp,log_dem


;log_t_mesh have been defined already......

        n_dem_mesh=n_elements(log_t_mesh)

;THIS is required: if you don't have DEM values at the extremes of the 
;range, the spline produces values that then create problems...

       log_dem_mesh=spline(log_t_f,log_dem_f ,log_t_mesh)

	missing_high=where(log_t_mesh gt max(log_t_f), nc )

;replace with the last value....
;--------------------------------
	if nc gt 0 then log_dem_mesh(missing_high)=$
		 log_dem_f(n_elements(log_t_f)-1)

	missing_low=where( log_t_mesh lt min(log_t_f),nc)

;replace with the first value....
;--------------------------------
	if nc gt 0 then log_dem_mesh(missing_low)= log_dem_f(0)

        dem_mesh=10.^log_dem_mesh

	oplot,log_t_mesh,log_dem_mesh,psym=5

endif else begin 

	log_dem_mesh=fltarr(n_elements(log_t_mesh))
	log_dem_mesh(*)=22.               ;default starting value for the DEM
	dem_mesh=10.^log_dem_mesh

	log_dem=spline(log_t_mesh,log_dem_mesh,log_dem_temp)
        dem=10.^log_dem

endelse


ch_wvl=fltarr(n_ch,n_obs)
ch_l1=intarr(n_ch,n_obs)
ch_l2=intarr(n_ch,n_obs)
ch_term=strarr(n_ch,n_obs)
ch_z=intarr(n_ch,n_obs)
ch_ion=intarr(n_ch,n_obs)
ch_id=strarr(n_ch,n_obs)

;this array will have the C(T) 
;-----------------------------
ch_contr_wa=fltarr(n_dem_temp,n_ch,n_obs)

;this array will have the G(T)=C(T)*abundances !!!
;------------------------------------------------- 
ch_contr=fltarr(n_dem_temp,n_ch,n_obs)

ch_pop=fltarr(n_dem_temp,n_ch,n_obs)
ch_n_contr=intarr(n_obs)
ch_contr_list=intarr(n_ch,n_obs)
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if keyword_set(file_gt) then begin

;start reading the contribution file
;-----------------------------------
	ioneq_name=' '
	const_net='' 
	lambda_obs=fltarr(1000)

	openr,lucon,file_gt,/get_lun
	readf,lucon,ioneq_name
	readf,lucon,const_net

	a=0.& i=0  
		while(not eof(lucon)) do begin
		readf,lucon,a
		lambda_obs(i)=a
		i=i+1
		endwhile
		free_lun,lucon

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


;EXCLUDE THE LINES WITH NO THEORETICAL COUNTERPARTS FROM THE FIT:
;----------------------------------------------------------------

	if exclu_obs_wvl_no_teo(0) ne 0. then begin

	print,'The following lines do not have theoretical counterparts '+$
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

;print,obs_wvl(iobs),ch_wvl(ilist,iobs),ch_id(ilist,iobs),$
;ch_z(ilist,iobs),ch_ion(ilist,iobs),ch_contr_wa(*,ilist,iobs),$
;ch_term(ilist,iobs), format=format 
    			endfor
  		endif
	endfor
	free_lun,lucon

;;end of the if you  have the contribution functions already 
endif    else begin 

;----------------------------------------------------------------------
;--- NOW WE START WITH GET_CONTRIBUTIONS ------------------------------
;----------------------------------------------------------------------

get_contributions,output=output,density=density, pressure=pressure,$
		cut_contrib=cut_contrib,n_ch=n_ch

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
;----------------------------------------------------------------------

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


;get the abundance file
;----------------------

get_ab_base=widget_base(title='Select the abundance file ',xsize=800,ysize=100)
 ;, /column,/frame,x_scroll=sz(0),y_scroll=sz(1))

;  the TOP widget is where all the user/numerical/text input takes place
;
top = widget_base(get_ab_base,/row,space=50)

; the top middle widget contains...
;
top_middle= widget_base(top,/column,space=50)

;  information/selection of abundance data
;
cs1 = widget_base(top_middle,/row)
junk = {CW_PDMENU_S,FLAGS:0,NAME:''}
af = FINDFILE(concat_dir(concat_dir(!xuvtop, 'abundance'), '*.abund'))

break_file,af,disk,dir,f1,ext
af = findfile('*.abund')
break_file,af,disk,dir,f2,ext

if f2(0) ne '' then f = [f1,'-----',f2] else f = f1

desc = [{CW_PDMENU_S,FLAGS:1,NAME:'Select Abundance file'}]
for i=0,n_elements(f)-1 do begin
  desc = [desc,{CW_PDMENU_S,FLAGS:0,NAME:f(i)}]
endfor
desc(n_elements(desc)-1).flags=2
menu = cw_pdmenu(cs1,desc,/return_name,uvalue='ABUND',font=font)
cs2 = widget_base(cs1,/column)
abund_info = cw_field(cs2,title='Current abundance file:',value='    ',$
            /row,xsize=20)

done=widget_button(cs2,UVALUE='D',VALUE='OK')

;  make the whole thing happen
;
widget_control,get_ab_base,/realize

xmanager,'get_ab',get_ab_base  , modal=modal, group_leader=group_leader

;--------END OF THE WIDGET  STUFF-----------------------------------------

;read the abundance file and print it
;------------------------------------
; the abundances are converted from  
; logaritmic values....abund(g)=10.^(abund(g)-12.)
;----------------------------------------------------------
read_abund,abund_name,abund,abund_ref
;
gz=where(abund gt 0.)
nz=n_elements(gz)

IF NOT keyword_set(quiet) THEN BEGIN 

 print,' '
 print,' abundances'
 print,' '

 for kz=0,nz-1 do begin
   jz=gz(kz)
   z2element,jz+1,ele
   print,jz+1,strpad(ele,12,/after),abund(jz),format='(i5,2x,a12,e10.2)'
 endfor
 print,' '

ENDIF

;

;this is necessary to write down the name of the ab file used for the dem:
;--------------------------------------------------------------------------
yesno=' '
read,'Do you want to change the abundances? [y/N] ',yesno
if yesno eq 'y' then begin
	new_abund_name=''
	read,'Enter the new core file name (a suffix ".abund" will be added) : ',$
		new_abund_name
	new_abund_name=new_abund_name+'.abund'
	spawn,'cp '+abund_name+' '+new_abund_name
	print,' '
	print,'The file " ',new_abund_name,'"  has been created in the working directory !'
	print,' '
	spawn,'chmod u+w '+new_abund_name
	editor=''
	read,'Which editor do you want to use ? (e.g. emacs) ',editor
	print,' '
	spawn,editor+' '+new_abund_name

;overwrite the string definition of the abundance file:
;-----------------------------------------------------
	abund_name=new_abund_name


; the abundances are converted from
; logaritmic values....abund(g)=10.^(abund(g)-12.)
;----------------------------------------------------------
	read_abund,abund_name,abund,abund_ref
;
	gz=where(abund gt 0.)
	nz=n_elements(gz)

IF NOT keyword_set(quiet) THEN BEGIN 

	print,' '
	print,' The new abundances are: ' 
	print,' '

	for kz=0,nz-1 do begin
   	jz=gz(kz)
   	z2element,jz+1,ele
   	print,jz+1,strpad(ele,12,/after),abund(jz),format='(i5,2x,a12,e10.2)'
	endfor

	print,' '
ENDIF

endif

;;MULTIPLY FOR THE ABUNDANCE FACTOR THE G(T)
;*********************************************

        ch_contr=fltarr(n_dem_temp,n_ch,n_obs)

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


;  plot out all contribution functions
;-------------------------------------

if keyword_set(plot_gt) then begin 


;cmax=max(ch_contr)
;plot_io,fltarr(2),xr=[dem_temp_min,dem_temp_max],yr=[cmax/1.e+3,cmax]

  for iobs=0,n_obs-1 do begin

     n_list=ch_n_contr(iobs)

	print,' '
	print,'Observed line (wavelength , identification) : '
	print,' '
	print,	trim (obs_wvl(iobs))+'  '+STRING(197b)+'  '+obs_id(iobs)
	print,' '
	print,'Theoretical lines (wavelength (A), ion,  G(T) max  ) :' 
	print,' '

	y_min=1.e-30
	y_max=1.e-22
	x_min=dem_temp_min
	x_max=dem_temp_max

	device,window_state=ws
	if(ws(11) ne 1) then  window,11,ysize=425,xsize=525  else wset,11

pr=''
begin_plot_gt:

	plot_io,fltarr(2),xr=[x_min,x_max],xstyle=1,ystyle=1,$
	yr=[y_min,y_max],xtitle = 'log T [ !eo!nK ]',$
	ytitle = ' G(T) [ erg cm!e3!n s!e-1!n st!e-1!n ]',$
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
		string(ch_wvl(ilist,iobs),'(f8.3)')+' '+string("305B)

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



  endfor   ;iobs

endif      ;,/plot_gt

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

;get the maximum and the temperature of the maximum of the summed G(T)
;----------------------------------------------------------------------
ch_tot_contr_max=fltarr(n_obs)
temp_max_tot_contr=fltarr(n_obs)

for iobs=0,n_obs-1 do begin
	ch_tot_contr_max(iobs)=max(ch_tot_contr(*,iobs))
	temp_max_tot_contr(iobs)=dem_temp_min+d_dem_temp*$
	( where (ch_tot_contr_max(iobs) eq (ch_tot_contr(*,iobs)) ))
endfor

;  plot out all summed contribution functions
;-------------------------------------
device,window_state=ws
if(ws(13) ne 1) then  window,13,ysize=425,xsize=525  else wset,13

	y_min=1.e-27
	y_max=1.e-22
	x_min=dem_temp_min
	x_max=dem_temp_max

pr=''
begin_plot_summed_gt:

	plot_io,fltarr(2),xr=[x_min,x_max],xstyle=1,ystyle=1,$
	yr=[y_min,y_max],xtitle = 'log T [ !eo!nK ]',$
	ytitle = 'summed G(T) [ erg cm!e3!n s!e-1!n st!e-1!n ]',$
	title='CHIANTI SPECTRAL CODE ',chars=1.2

for iobs=0,n_obs-1 do begin
 
 	        oplot,log_dem_temp, ch_tot_contr(*,iobs)

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;ch_int=fltarr(n_obs)
;for iobs=0,n_obs-1 do ch_int(iobs)=total(ch_tot_contr(*,iobs)*dem*dem_temp)*dlnt
;
;  summarize contributions to each line
;----------------------------------------
;for iobs=0,n_obs-1 do begin
;     print, obs_id(iobs),obs_wvl(iobs),obs_int(iobs),obs_sig(iobs),$
;	   ch_int(iobs), obs_int(iobs)/ch_int(iobs)
;     n_list=ch_n_contr(iobs)
;     for ilist=0,n_list-1 do begin
;        	print,ch_wvl(ilist,iobs),ch_id(ilist,iobs),$
;			ch_z(ilist,iobs),ch_ion(ilist,iobs), $
;               		format='(10x,f10.3,a15,2i5)'
;
;     endfor
;endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
flambda=10.
scale=1.
failed=0

   dem_fit,y,chisqr,flambda=flambda,scale=scale,$
	 niter=niter,dchisq_min=dchisq_min,$
	failed=failed ,quiet=quiet

if failed eq 1 then begin 
message,' Sorry, the fit failed, you might have to start again with a new selection of DEM/points ! ',/continue

yesno=''
read,' Do you want to see the results anyway ? [Y/n] ',yesno
if strlowcase(yesno) eq 'n' then return

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; overplot the final DEM:
;------------------------


;treat negative values:
;***********************

log_dem_mesh=log_dem_mesh >0

  log_dem=spline(log_t_mesh,log_dem_mesh,log_dem_temp)
   dem=10.^log_dem 

  bad= where( log_dem lt 0.,nbad)
  if nbad gt 0 then begin 

  log_dem(bad)=0
  dem(bad)=0
  endif


;oplot,dem_temp,dem,thick=3



;Create the output DEM file to be used as input by CHIANTI_SS.PRO (DEM(T)):
;---------------------------------------------------------***************

dem_name=output+'.dem'
openw,lun_dem,dem_name,/get_lun

	for i=0,n_dem_temp-1 do $
  	  printf,lun_dem,log_dem_temp(i),log_dem(i)    ; + log_dem_temp(i)    

  printf,lun_dem,'-1'
  printf,lun_dem,'%file:  ',output+'.dem'
  printf,lun_dem,'%dem:   Produced by CHIANTI_DEM.PRO '
  printf,lun_dem,'%With the abundance file ',abund_name
  printf,lun_dem,'%the observation file "',file_input_obs
  printf,lun_dem,'% calculated at '+const_net
  printf,lun_dem,'% ADD YOUR NOTES HERE'
  printf,lun_dem,'-1'

free_lun,lun_dem

IF NOT keyword_set(quiet) THEN BEGIN
;
print, ' '
print,' Row   ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
print, ' '

for i=0,n_obs-1 do begin
   print,i,obs_id(i),obs_wvl(i),obs_int(i),y(i),obs_sig(i),$
	(y(i)-obs_int(i))^2/obs_sig(i)^2, obs_int(i)/y(i),$
                format='(i5,1x,a20,1f10.3, 3e10.3, 2f10.3 )'
endfor

print, ' '

ENDIF


;
;  -----------------------------------------------------
;ch_z(ch_n_contr(iobs),iobs)=iz
;            ch_ion(ch_n_contr(iobs),iobs)=ion

;Open the output file where all the info will be stored:
;-------------------------------------------------------

output_general=output+'.general'
openw,luo,output_general,/get_lun
;
printf,luo,' abundance file = '+abund_name
printf,luo,' ionization equilibrium file = '+ioneq_name
printf,luo,const_net

IF NOT keyword_set(quiet) THEN BEGIN
print,' '
print,'------------------------------------------------ '
print,'      not  sorted '
print,'------------------------------------------------ '
print,' '
print,'                       I_obs      I_theo    sigma_obs'
print,' Row     ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
print,' '

ENDIF
;
;
printf,luo,' '
printf,luo,'------------------------------------------------ '
printf,luo,'     not sorted  '
printf,luo,'------------------------------------------------ '
printf,luo,' '

printf,luo,'                          I_obs      I_theo    sigma_obs'
printf,luo,' Row     ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
printf,luo,' '

;
for iobs=0,n_obs-1 do begin

IF NOT keyword_set(quiet) THEN BEGIN
   print, ' '
   print, ' '
   print,iobs,obs_id(iobs),obs_wvl(iobs),obs_int(iobs),y(iobs),$
         obs_sig(iobs),(y(iobs)-obs_int(iobs))^2/obs_sig(iobs)^2,$
         y(iobs)/obs_int(iobs), $
        format='(i5,1x,a20,1f10.3,1x,3e10.3,3f10.3)'
;         format='(i5,2x,a20,5f10.3,20x,f10.3)'
ENDIF

   printf,luo, ' '
   printf,luo, ' '
   printf,luo,iobs,obs_id(iobs),obs_wvl(iobs),obs_int(iobs),y(iobs),$
          obs_sig(iobs),(y(iobs)-obs_int(iobs))^2/obs_sig(iobs)^2, $
          y(iobs)/obs_int(iobs), $
        format='(i5,1x,a20,1f10.3,1x,3e10.3,3f10.3)'
;          format='(i5,2x,a20,5f10.3,20x,f10.3)'
;
   n_list=ch_n_contr(iobs)

    for ilist=0,n_list-1 do begin

        this_contr=total(ch_contr(*,ilist,iobs)*dem*dem_temp)*dlnt
        this_contr=this_contr/y(iobs)


	IF NOT keyword_set(quiet) THEN BEGIN

   	if ilist eq 0 then $
	print, '------------------------------------------------------------------------------------------------'
	  
          print,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
          ch_term(ilist,iobs),this_contr,$
	  format='(10x,a15,f10.3,a40,f10.3)'
;
	ENDIF

	if ilist eq 0 then $
	printf,luo, '------------------------------------------------------------------------------------------------'

        printf,luo,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
        ch_term(ilist,iobs),this_contr,$
	format='(10x,a15,f10.3,a40,f10.3)'

    endfor

printf,luo,' '

;
endfor


log_dem_step=spline(log_t_mesh,log_dem_mesh ,log_t_dem_step)
dem_step  =10.^log_dem_step


ch_dom_z=intarr(n_obs)
ch_dom_ion=intarr(n_obs)

temp_gt_dem=fltarr(n_obs)
dem_temp_gt_dem=fltarr(n_obs)

dom_ions_spec=strarr(n_obs)
ch_dom_wave=fltarr(n_obs)

IF NOT keyword_set(quiet) THEN BEGIN

;-------------------------------------------------------------
print,' '
print,'------------------------------------------------ '
print,' sorted by dominant contributing element and ion'
print,'------------------------------------------------ '
print,' '
print,'			      I_obs	 I_theo	   sigma_obs'
print,' Row     ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
print,' '

ENDIF

;
printf,luo,' '
printf,luo,'------------------------------------------------ '
printf,luo,' sorted by dominant contributing element and ion'
printf,luo,'------------------------------------------------ '
printf,luo,' '

printf,luo,'			      I_obs	 I_theo	   sigma_obs'
printf,luo,' Row     ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
printf,luo,' '
;

;
for iobs=0,n_obs-1 do begin
   n_list=ch_n_contr(iobs)
;
   if n_list gt 1 then begin
      dom=fltarr(n_list)
      for ilist=0,n_list-1 do begin
         dom(ilist)=total(ch_contr(*,ilist,iobs)*dem*dem_temp)*dlnt
      endfor

      idom=where(max(dom) eq dom)

      ch_dom_z(iobs)=ch_z(idom,iobs)
      ch_dom_ion(iobs)=ch_ion(idom,iobs)

	ch_dom_wave(iobs)=ch_wvl(idom,iobs)

	prod_gt_dem=  spline ( dem_temp, ch_contr(*,idom,iobs)*dem ,t_dem_step)

	temp_gt_dem (iobs) = t_dem_step(where (prod_gt_dem eq max (prod_gt_dem)))
	dem_temp_gt_dem(iobs) = dem_step(where (prod_gt_dem eq max (prod_gt_dem)))


;ch_dom_temp(iobs)=where(max(ch_contr(*,idom,iobs)*dem) eq ch_contr(*,idom,iobs)*dem)

;
   endif else begin

      ch_dom_z(iobs)=ch_z(0,iobs)
      ch_dom_ion(iobs)=ch_ion(0,iobs)
	ch_dom_wave(iobs)=ch_wvl(0,iobs)

      if max(ch_contr(*,0,iobs)) gt 0. then begin

	prod_gt_dem=  spline ( dem_temp, ch_contr(*,0,iobs)*dem ,t_dem_step)

	temp_gt_dem (iobs) = t_dem_step(where (prod_gt_dem eq max (prod_gt_dem)))
	dem_temp_gt_dem(iobs) = dem_step(where (prod_gt_dem eq max (prod_gt_dem)))

; ch_dom_temp(iobs)=where(max(ch_contr(*,0,iobs)*dem) eq ch_contr(*,0,iobs)*dem)


      endif
   endelse
;

endfor   ; iobs


;
measure=intarr(n_obs)
;
;
for iobs=0,n_obs-1 do begin
    measure(iobs)=ch_dom_z(iobs)*100+ch_dom_ion(iobs)
endfor
isort=sort(measure)
;

for j=0,n_obs-1 do begin

   iobs=isort(j)

IF NOT keyword_set(quiet) THEN BEGIN

   print, ' '
   print, ' '
   print,iobs,obs_id(iobs),obs_wvl(iobs),obs_int(iobs),y(iobs),$
	 obs_sig(iobs),(y(iobs)-obs_int(iobs))^2/obs_sig(iobs)^2,$
	 y(iobs)/obs_int(iobs), $
        format='(i5,1x,a20,1f10.3,1x,3e10.3,3f10.3)'
;         format='(i5,2x,a20,5f10.3,20x,f10.3)'
;
ENDIF


   printf,luo, ' '
   printf,luo, ' '
   printf,luo,iobs,obs_id(iobs),obs_wvl(iobs),obs_int(iobs),y(iobs),$
	 obs_sig(iobs),(y(iobs)-obs_int(iobs))^2/obs_sig(iobs)^2, $
	 y(iobs)/obs_int(iobs), $
        format='(i5,1x,a20,1f10.3,1x,3e10.3,3f10.3)'
;         format='(i5,2x,a20,5f10.3,20x,f10.3)'

   n_list=ch_n_contr(iobs)

   for ilist=0,n_list-1 do begin

      	this_contr=total(ch_contr(*,ilist,iobs)*dem*dem_temp)*dlnt
      	this_contr=this_contr/y(iobs)



	IF NOT keyword_set(quiet) THEN BEGIN

   	if ilist eq 0 then $
	  print, '------------------------------------------------------------------------------------------------'

      	  print,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
	      ch_term(ilist,iobs),this_contr,$
	      format='(10x,a15,f10.3,a40,f10.3)'
;
	ENDIF

   	if ilist eq 0 then $
	printf,luo, '------------------------------------------------------------------------------------------------'

	printf,luo,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
	      ch_term(ilist,iobs),this_contr,$
	format='(10x,a15,f10.3,a40,f10.3)'


   endfor   ; ilist
;

zion2spectroscopic,ch_dom_z(iobs),ch_dom_ion(iobs),dom_ion

dom_ions_spec(iobs)=dom_ion


;
endfor ; j:n_obs
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;observed/theoric ratio :
 ot_ratio = obs_int /y

; this is    alog10 (DEMeff * ( (Iobs + err_obs ) / I_theo)) 
erru=  alog10 ( dem_temp_gt_dem *  (obs_int + obs_sig )/ y )

errl= alog10 ( dem_temp_gt_dem *  (obs_int - obs_sig )/ y )


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


device,window_state=ws
if(ws(2) ne 1) then  window,2,ysize=500,xsize=600 else wset,2
;
dem=10.^spline(log_t_mesh,log_dem_mesh,log_dem_temp)

x_min=(dem_temp_min)
x_max=(dem_temp_max)
y_min=19.
y_max=24.

print, ' ' 
print,'identification from the input file ---- and dominant element :'
;---------------------------------------------------------------
print, ' ' 

for iobs=0,n_obs-1 do $

print, string ( trim(obs_id(iobs))+' '+strtrim(obs_wvl(iobs),2)+' '+STRING(197b), format='(a40)')+$
	'   ----   '+$
    trim (dom_ions_spec(iobs))+' '+ strtrim(ch_dom_wave(iobs),2)+' '+STRING(197b)

print, ' ' 

pr=''
begin_plot_dem:


	plot, log_t_dem_step, log_dem_step ,/ynozero,xr=[x_min,x_max],yr=[y_min,y_max],$
	xstyle=1,ystyle=1,$
	xtitle='log T!Dmax!n ( G(T)*DEM(T) ) [ !eo!nK ]',$
	ytitle='log DEM [ cm!S!E-5 !NK!S!E-1!N ]',$
	title='CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',chars=1.2

   oplot,alog10( temp_gt_dem ),alog10( dem_temp_gt_dem * ot_ratio),psym=1

   errplot,alog10( temp_gt_dem ), erru, errl



if pr ne 'y' then begin 
ll='' 
read, 'Do you want to overplot the labels from the input file [i] or the dominant ion [l] ? ' ,ll
if ll ne 'i' and ll ne 'l' then goto,begin_plot_dem
endif


if ll eq 'i' then  begin

  for iobs=0,n_obs-1 do begin
	if obs_id(iobs) ne '' then $
	xyouts,alog10( temp_gt_dem (iobs)  ),$
		alog10( dem_temp_gt_dem (iobs) * ot_ratio(iobs) ),$
	'  '+trim(obs_id(iobs)),$  
; +' '+strtrim(obs_wvl(iobs),2)+' '+STRING(197b),$
		charsize=1.2  ;,orientation=90
  endfor
endif

if ll eq 'l' then begin

  for iobs=0,n_obs-1 do begin
	if obs_id(iobs) ne '' then $
	xyouts,alog10( temp_gt_dem (iobs)  ),$
		alog10( dem_temp_gt_dem (iobs) * ot_ratio(iobs) ),$
		'  '+trim (dom_ions_spec(iobs)) ,$
;  +' '+ strtrim(ch_dom_wave(iobs),2)+' '+STRING(197b),$
		charsize=1  ;,orientation=90
  endfor
endif


       break_file,dem_name, disk,dir,f,ext
       xyouts,0.500,0.90,/normal,'DEM file  : '+f,chars=1.2
       break_file,abund_name, disk,dir,f,ext
       xyouts,0.500,0.87,/normal, 'Abund.  : '+f,chars=1.2 
       break_file,ioneq_name, disk,dir,f,ext
       xyouts,0.500,0.84,/normal, 'Ion Eq. : '+f,chars=1.2 
       xyouts,0.500,0.81,/normal, const_net ,   chars=1.2
       xyouts,0.500,0.78,/normal,              chars=1.2,$
		' !4v!3!e2!n  = '+string(chisqr) ;format='(f4.1)')

print2d_plot, pr=pr, x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
	  go_to_line=go_to_line,out_name=output+'_dem.ps'   

if go_to_line eq 'y' then goto,begin_plot_dem


IF NOT keyword_set(quiet) THEN BEGIN

print,'-------------------------------------- '
print,' sorted by observed wavelength'
print,'-------------------------------------- '
print,' '
print,'			      I_obs	 I_theo	   sigma_obs'
print,' Row     ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
print,' '

ENDIF

;
printf,luo,' '
printf,luo,'-------------------------------------- '
printf,luo,' sorted by observed wavelength'
printf,luo,'-------------------------------------- '
printf,luo,' '
printf,luo,'			      I_obs	 I_theo	   sigma_obs'
printf,luo,' Row     ident. obs.lambda  obs.flux  calc.flux obs.err     (I_theo-I_obs )/sigma_obs I_theo/I_obs'
printf,luo,' '
;
measure=obs_wvl
;
;
isort=sort(measure)
;
for j=0,n_obs-1 do begin

   iobs=isort(j)

IF NOT keyword_set(quiet) THEN BEGIN

   print, ' '
   print, ' '
   print,iobs,obs_id(iobs),obs_wvl(iobs),obs_int(iobs),y(iobs),$
	obs_sig(iobs),(y(iobs)-obs_int(iobs))^2/obs_sig(iobs)^2, $
	y(iobs)/obs_int(iobs), $
        format='(i5,1x,a20,1f10.3,1x,3e10.3,3f10.3)'

ENDIF

;
   printf,luo, ' '
   printf,luo, ' '
   printf,luo,iobs,obs_id(iobs),obs_wvl(iobs),obs_int(iobs),y(iobs),$
	obs_sig(iobs),(y(iobs)-obs_int(iobs))^2/obs_sig(iobs)^2, $
	y(iobs)/obs_int(iobs), $
        format='(i5,1x,a20,1f10.3,1x,3e10.3,3f10.3)'
;
   n_list=ch_n_contr(iobs)

   for ilist=0,n_list-1 do begin

      	this_contr=total(ch_contr(*,ilist,iobs)*dem*dem_temp)*dlnt
      	this_contr=this_contr/y(iobs)



	IF NOT keyword_set(quiet) THEN BEGIN

   	if ilist eq 0 then $
	  print, '------------------------------------------------------------------------------------------------'

      	  print,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
	     ch_term(ilist,iobs),this_contr,$
		format='(10x,a15,f10.3,a40,f10.3)'
;
	ENDIF

   	if ilist eq 0 then $
	printf,luo, '------------------------------------------------------------------------------------------------'

     	printf,luo,strpad(ch_id(ilist,iobs),12,/after),ch_wvl(ilist,iobs),$
	     ch_term(ilist,iobs),this_contr,$
	format='(10x,a15,f10.3,a40,f10.3)'
;
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

printf,luo,' '
printf,luo,' chisqr = ',chisqr
print,' dem chisqr = ',chisqr

free_lun,luo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

device,window_state=ws
if(ws(3) ne 1) then  window,3,ysize=650,xsize=900 else wset,3

pr=''
begin_post4:

!p.multi=[0,2,2,0,0]
!p.title=' '
;!p.thick=0.6
;!p.charsize=0.8

;x_min=(dem_temp_min)
;x_max=(dem_temp_max)
;y_min=19. & y_max=24.


plot, log_t_dem_step , log_dem_step , /ynozero,xrange=[x_min,x_max],$
	xtitle='log T!Dmax!n ( G(T)*DEM(T) ) [ !eo!nK ]',$
	ytitle='log DEM  [cm!S!E-5 !NK!S!E-1!N]',xstyle=1,ystyle=1,$
	yrange=[y_min,y_max],thick=1,$
;        title='CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',$
	xmargin=[10,4],ymargin=[3,3],chars=1.2


   	oplot,alog10( temp_gt_dem ) ,$
		alog10(dem_temp_gt_dem * ot_ratio),psym=1

	errplot,alog10(temp_gt_dem ), erru, errl


for iobs=0,n_obs-1 do begin

	xyouts,alog10( temp_gt_dem (iobs)),$
		alog10(dem_temp_gt_dem(iobs)  * ot_ratio(iobs)),$
	'  '+trim(obs_id(iobs)), charsize=1  ;,orientation=90

endfor

;;;;;;;;;;;;;;;;;;;;
;II plot
;--------
;x_min=(dem_temp_min)   &x_max=(dem_temp_max)
x1=min(temp_max_tot_contr) & x2=max(temp_max_tot_contr)
y1=min(alog10( y/obs_int )) &y2=max(alog10( y/obs_int ))

plot,[0], [0],symsize=.8,thick=1,/nodata,$
	xrange=[x1,x2],yr=[y1,y2],$
	xtitle="log T!Dmax!N ( G(T) ) [ !eo!nK ]",$
	ytitle=" log (I!Dtheo!N / I!Dobs!N)",psym=1,$
;	title='CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',$
	xmargin=[10,4],ymargin=[3,3],chars=1.2


for iobs=0,n_obs-1 do begin
    oplot,[temp_max_tot_contr(iobs)], $
	[alog10( y(iobs))- alog10(obs_int(iobs))],symsize=0.8, thick=1,psym=1
endfor

;;;;;;;;;;;;;;;;;;;;
;III plot
;--------
x1=min(obs_wvl) &x2=max(obs_wvl)
y1=min(alog10( y/obs_int )) &y2=max(alog10( y/obs_int ))

plot,[obs_wvl(0)],[0.],/nodata,/ynozero,$
	symsize=.3,thick=1,psym=1,$
	xrange=[x1,x2],yr=[y1,y2],$
	xtitle=" wavelength [ "+string("305B)+" ]",$
	ytitle=" log (I!Dtheo!N / I!Dobs!N)",$
;	title='CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',$
	xmargin=[10,4],ymargin=[6,2],chars=1.2


for iobs=0,n_obs-1 do begin
      oplot, [obs_wvl(iobs)],[alog10( y(iobs))- alog10(obs_int(iobs))],$
	     symsize=0.8, thick=1,psym=1
endfor


;;;;;;;;;;;;;;;;;;;;
;IV plot
;-------
x1=min(alog10( y )) &x2=max(alog10( y ))
y1=min(alog10(obs_int-obs_sig )) &y2=max(alog10( obs_int+obs_sig ))

plot,[alog10( y(0)) ],[alog10(obs_int(0))],/nodata,$
	xrange=[x1,x2],yr=[y1,y2],$ 
	xtitle='log I!Dtheo !N [ergs cm!S!E-2 !Ns!E-1 !Nst!E-1 !N]',$ 
	ytitle=" log I!Dobs !N [ergs cm!S!E-2 !Ns!E-1 !Nst!E-1 !N]",$
;	title='CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',$
	xmargin=[10,4],ymargin=[6,2],chars=1.2

for iobs=0,n_obs-1 do begin
	oplot, [alog10( y(iobs))], [alog10(obs_int(iobs))],$
		symsize=0.8, thick=1,psym=1
	errplot, [alog10( y(iobs))], [alog10(obs_int(iobs)-obs_sig(iobs))],$
		 [alog10(obs_int(iobs)+obs_sig(iobs))]
endfor

xyouts,0.5,0.97,/norm,'CHIANTI SPECTRAL CODE AND INVERSION TECHNIQUE',$
	charsize=1.4 , align=0.5,orientation=0

       break_file,dem_name, disk,dir,f,ext
       xyouts,0.2,0.91,/norm,'DEM file  : '+f,chars=1
       break_file,abund_name, disk,dir,f,ext
       xyouts,0.2,0.89,/norm, 'Abund.  : '+f,chars=1
       break_file,ioneq_name, disk,dir,f,ext
       xyouts,0.2,0.87,/norm, 'Ion Eq. : '+f,chars=1
       xyouts,0.2,0.85,/norm,chars=1, const_net
       xyouts,0.2,0.83,/norm,chars=1,$
                ' !4v!3!e2!n  = '+string(chisqr)  ;,format='(f4.1)')



print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
	  go_to_line=go_to_line, out_name=output+'_4plots.ps'   

if go_to_line eq 'y' then goto,begin_post4


;write the output with the useful information for user-written plotting 
;routines.

openw,lu,output+'.out',/get_lun

for iobs=0,n_obs-1 do begin

   printf,lu, obs_id(iobs), obs_wvl(iobs), obs_int(iobs), y(iobs),$
          obs_sig(iobs), alog10(temp_gt_dem(iobs) ), $
	alog10( dem_temp_gt_dem (iobs)),$
          format='(a20,1x, 1f10.3,1x, 3e10.3, 1x,  f4.2,1x,f6.3)'

endfor

free_lun,lu



return
end


