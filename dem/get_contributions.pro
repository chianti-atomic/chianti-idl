;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas.  See www.chiantidatabase.org 
;
;     		          
; Purpose     : calculates the contribution functions G(T) at constant 
;		pressure or density  of the lines present in the CHIANTI
;		database, corresponding to  a given set of observed lines.
;
; Category    :
;               
; Explanation : This routine is called by CHIANTI_DEM. It cannot be used as
;		a stand-alone routine.
;		The observation file is read by CHIANTI_DEM.
;		GET_CONTRIBUTIONS starts reading 
;		the ionization equilibrium file and the masterlist of the
;		ions present in the CHIANTI database.
;		GET_CONTRIBUTIONS  then searches the 
;		CHIANTI database (ion per ion) for all the 
;		theoretical lines corresponding to the observed	lines, i.e. 
;		that lie in a OBS_WVL(i) +/- DELTA_LAMBDA_OBS(i) interval
;		centered on the observed wavelength OBS_WVL(i).
;
;		A constant pressure OR a constant density for all the lines
; 		is used. If you select a constant pressure,
;               for each ion the contribution function is calculated at an 
;               electron density N_e equal to the ratio of the pressure 
;               and the temperature of maximum ionization fraction:  
;               C=C( T, N_e = P/T_ion)
;               The C(T) values are stored by GET_CONTRIBUTIONS in the output 
;		file OUTPUT.CONTRIBUTIONS. The temperature grid
;		follows the ionization equilibrium file.
;
;		In the case no theoretical lines corresponding to an observed
;		line are found, the routine writes the wavelength of the line
;		to be excluded from the fit in the array EXCLU_OBS_WVL_NO_TEO;
;		these lines are then excluded from the fitting  by 
;	 	CHIANTI_DEM. You might consider the possibility to start again
;		incrementing the DELTA_LAMBDA_OBS, to see if there are 
;		theoretical lines in the vicinity.
;
;
; Use         : called by CHIANTI_DEM to calculate the contribution functions
;
; Examples    : 
;
;    
; Inputs      :	Various, in form of keywords.
;               
;               
; Opt. Inputs : none
;               
; Outputs     : OUTPUT.CONTRIBUTIONS  
;		Is the file where all the contribution  functions G(T) are 
;		stored. In the first two lines  the ionization equilibrium 
;		file name, and the constant value of pressure or density 
;		adopted are reported. Then for each line you have reported  
;               the observed wavelength, the theoretical one, the element and
;		ionization stage, then the C(T) values. At the end the 
;		specification for each transition.
;
;
; Opt. Outputs: None
;               
; Keywords    : (all passed by CHIANTI_DEM)
;
;
;	CUT_GT:			if set, only those
;				theoretical lines that have a MAX(G(T)) greater
;      			        than the value set, are kept; it is useful to 
;				set this value in order to reduce the number 
;				of lines in the file where the G(T) are stored.
;
;	DENSITY : 		the value of the density (Ne).
;
;	FILE_INPUT:		if set, you are not 
;				requested to select the observation file.
;
;  
;	N_MATCHES:              
;		In the unlikely event that more than 100 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;
;	PRESSURE:     	the value of the pressure (Ne T)
;
;	OUTPUT  :  	the core name for the output
;
;
; Calls       : 
;		ch_synthetic
;
; Common      : 	      
;		obs, 	obs_int,obs_sig,n_obs
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
;		 ch_pop,ch_contr_list, ch_term,ch_n_contr
;
;
; Restrictions: ;
;		THIS IS NOT A STAND-ALONE PROCEDURE. 
;		It is called by CHIANTI_DEM,
;		and has a lot of common blocks with other procedures.
;
;		In the unlikely event that more than 100 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;
;               
; Side effects: None known.
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
;       Version 1, Giulio Del Zanna (GDZ)  5 November  1997 
;	UCLAN (University of Central Lancashire, UK) 
;
; Modified    :  
;       Version 2, 31-Oct-2000, GDZ, DAMTP.  Rewritten completely the routine,
;       to make it compatible with CHIANTI v.3. Based the core calculations on
;       new implementations due to Peter Young, CfA. 
;
;       Version 3, 5-Dec-2000, GDZ, DAMTP. Fixed a bug when checking the 
;       values in the .splups files.
;
;
;       Ver. 4,  25-Apr-02, GDZ  
;              Revised to account for v.4 variations. By default the proton
;              rates are included in the calculation of the level population.
;
;       V.5, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V. 6,  10-July-2002  GDZ
;                  Corrected a bug. It now properly includes by default the
;                  proton rates in the population solver. 
;
;       V. 7,  4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure. 
;
;       V. 8, 3-Nov-03  GDZ
;                  Modified format e8.2 to e9.2 for Windows compatibility.
;
;       V. 9, 4-May-05 Enrico Landi (EL)
;                  Modified in order to include ionization and recombination
;                  data in the input to POP_SOLVER
;
;       V.10, 12-Jun-2009, Enrico Landi
;                  Changed the definition of the temperature array for ion fractions
;                  in the IONREC variable, now taken directly from the output of
;                  READ_IONEQ.PRO
;
;       V.11, 10-Mar-2014 Giulio Del Zanna (GDZ)
;                  Modified substantially by calling ch_synthetic.
;
;       V.12,  5-Dec-2018, GDZ
;                  Added another output, a savegen file with a
;                  structure having all the G(T).
;
; VERSION     :    V.12
;
;
;-        
;---------------------------

pro get_contributions,output=output,density=density, pressure=pressure,$
                      cut_contrib=cut_contrib,n_ch=n_ch


;these are the commons with CHIANTI_DEM:

common obs, 	obs_int,obs_sig,n_obs
common obs_o,	obs_wvl,obs_id,obs_delta_lambda
common dem, 	d_dem_temp,dem_temp,log_dem_temp,log_t_mesh,log_dem_mesh
common contr,	ch_tot_contr
common ab,	abund_name,abund_info,xuvtop,ioneq_name


common various,	exclu_obs_wvl_no_teo,const_net,$
   dem_temp_min,dem_temp_max,n_dem_temp,$
   ch_wvl,ch_l1,ch_l2,ch_id,ch_z,ch_ion,ch_contr_wa,$
   ch_pop,ch_contr_list, ch_term,ch_n_contr


;ON_ERROR, 2

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
   message, 'system variable !xuvtop must be set  '
xuvtop = !xuvtop


;create an array  where eventual observed wavelength that DO NOT have 
;------------------------------------------------------------
;corresponding theoretical lines (at maximum can have an n_obs length)
;------------------------------------------------------------
;will be stored.
;---------------

exclu_obs_wvl_no_teo=fltarr(n_obs)


;pick up the Ionization Equilibrium File:
;----------------------------------------


ioneq_name=pickfile(path= concat_dir(!xuvtop, '/ioneq'),filter='*.ioneq',title='Select Ionization Equilibrium File')

ff = findfile(ioneq_name)
IF  ff(0)  NE ''  THEN $
   read_ioneq,ioneq_name,ioneq_t,ioneq,ioneq_ref ELSE BEGIN 
   message, 'Error,  no ioneq file found !'
   err_msg = 'Error,  no ioneq file found !'
   return
END 


n_dem_temp=n_elements(ioneq_t)

; inizialise the arrays:                
;----------------------

ch_wvl=fltarr(n_ch,n_obs)
; ch_l1=intarr(n_ch,n_obs)
; ch_l2=intarr(n_ch,n_obs)
ch_term=strarr(n_ch,n_obs)

ch_z=intarr(n_ch,n_obs)
ch_ion=intarr(n_ch,n_obs)
ch_id=strarr(n_ch,n_obs)

;; ;this array will have the C(T) 
;; ;-----------------------------
ch_contr_wa=fltarr(n_dem_temp,n_ch,n_obs)

;; ;this array will have the G(T)=C(T)*abundances !!!
;; ;------------------------------------------------- 
;  ch_contr=fltarr(n_dem_temp,n_ch,n_obs)

ch_n_contr=intarr(n_obs)


; get the minimum and maximum wavelengths

iobs=0
w1=obs_wvl[iobs]-obs_delta_lambda[iobs]
w2=obs_wvl[iobs]+obs_delta_lambda[iobs]

FOR iobs=1, n_obs-1 DO BEGIN
   
   w1 = w1 < (obs_wvl[iobs]-obs_delta_lambda[iobs])
   w2 = w2 > (obs_wvl[iobs]+obs_delta_lambda[iobs])

ENDFOR

print,' getting CHIANTI contribution functions G(T)  for ALL the  observed lines'

ch_synthetic, w1, w2, output=output_ch, err_msg=err_msg, msg=msg, $
              pressure=pressure, density=density, $
              all=all, ioneq_name=ioneq_name, $
              noprot=noprot, rphot=rphot, radtemp=radtemp, /goft



;      check to see if there any wavelength matches -  loop over the obs lines.
;---------------------------------------------------------------
anylines = -1


FOR iobs=0, n_obs-1 do BEGIN
   
   w1=obs_wvl[iobs]-obs_delta_lambda[iobs]
   w2=obs_wvl[iobs]+obs_delta_lambda[iobs]

; select the lines 
   
   anylines=where(output_ch.lines.wvl ge w1 and output_ch.lines.wvl le w2, nn)
   if nn eq 0 then begin 
      print, 'no CHIANTI lines for the observed '+string(obs_wvl[iobs])+' - remove the line! '
      stop
      return
   endif else if nn  GE  n_ch then begin
      message,'Error,  increase n_matches to at least '+string(nn)
   ENDIF else begin 
      
; GDZ Addition 
        IF iobs EQ 0 THEN lines = {obs0:output_ch.lines[anylines]} else $
               lines = create_struct(lines, 'obs'+trim(iobs), output_ch.lines[anylines])

      n_lines_this_iobs = 0

      FOR i=0,nn-1 DO BEGIN     ; for each contributing line within iobs       
         

         n_lines_this_iobs =n_lines_this_iobs +1

         ch_wvl(ch_n_contr(iobs),iobs)=output_ch.lines[anylines[i]].wvl
         
         ch_id(ch_n_contr(iobs),iobs)= output_ch.lines[anylines[i]].snote
         ch_z(ch_n_contr(iobs),iobs)=output_ch.lines[anylines[i]].iz
         ch_ion(ch_n_contr(iobs),iobs)=output_ch.lines[anylines[i]].ion


         l1 = output_ch.lines[anylines[i]].lvl1
         l2 = output_ch.lines[anylines[i]].lvl2

         ch_term(ch_n_contr(iobs),iobs)=output_ch.lines[anylines[i]].ident 


         ch_contr_wa(*,ch_n_contr(iobs),iobs)=output_ch.lines[anylines[i]].goft

; ch_n_contr is an  intarr(n_obs) defined within CHIANTI_DEM:

         ch_n_contr(iobs)=ch_n_contr(iobs)+1

      ENDFOR                    ; for each contributing line within iobs

; ch_contr_list is an   intarr(n_ch,n_obs) defined within CHIANTI_DEM, but not needed.

   ENDELSE                      ;nn gt 0

ENDFOR                          ;  iobs


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; GDZ Addition 

out = rem_tag(output_ch, 'LINES')

out =  add_tag(out,  obs_wvl, 'all_obs_wvl')
out =  add_tag(out, obs_delta_lambda, 'all_obs_delta_lambda')

out =  add_tag(out, lines, 'LINES')

savegen, file=output+'_gt.genx', st=out
delvarx,out
print, 'savefile '+output+'_gt.genx  saved ! '


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;

output_contributions=output+'.contributions'
openw,lucontr,output_contributions,/get_lun

printf,lucontr,ioneq_name
if keyword_set(pressure) then  begin
   const_net= 'constant pressure: '+$
              string(format='(e9.2)',pressure) +' [ cm!e-3!n !eo!nK ]' 
   printf,lucontr,const_net
endif else begin
   if keyword_set(density)  then  begin
      const_net=' constant density: '+$
                string(format='(e9.2)',density)+' [ cm!e-3!n ]'
      printf,lucontr,const_net
   endif
endelse
;

ch_contr_max=fltarr(n_ch,n_obs)
for iobs=0,n_obs-1 do begin
   n_list=ch_n_contr(iobs)
   if n_list ge 1 then begin
      for ilist=0,n_list-1 do begin
         ch_contr_max(ilist,iobs)=max(ch_contr_wa(*,ilist,iobs))
      endfor
   endif
   
endfor

;Cut out those lines that have max(G(T)) less than cut_contrib
;------------------------------------------------------------

if cut_contrib ne 0 then begin
   for iobs=0,n_obs-1 do begin

      keep= where (ch_contr_max(*,iobs) ge cut_contrib)
      if max(keep) ge 0 then begin

         ch_n_contr(iobs) = n_elements(keep)

         ch_wvl(0:ch_n_contr(iobs)-1,iobs)=ch_wvl(keep,iobs)
         ch_wvl(ch_n_contr(iobs):n_ch-1,iobs)=0.

         ch_id(0:ch_n_contr(iobs)-1,iobs)=ch_id(keep,iobs)
         ch_id(ch_n_contr(iobs):n_ch-1,iobs)=' '

         ch_z (0:ch_n_contr(iobs)-1,iobs)   =ch_z(keep,iobs)
         ch_z(ch_n_contr(iobs):n_ch-1,iobs)=0

         ch_ion (0:ch_n_contr(iobs)-1,iobs) =ch_ion(keep,iobs)
         ch_ion(ch_n_contr(iobs):n_ch-1,iobs)=0

         ch_contr_wa(*,0:ch_n_contr(iobs)-1,iobs) = ch_contr_wa(*,keep,iobs)
         ch_contr_wa(*,ch_n_contr(iobs):n_ch-1,iobs)=0.

         ch_term (0:ch_n_contr(iobs)-1,iobs) =ch_term(keep,iobs)
         ch_term (ch_n_contr(iobs):n_ch-1,iobs) =' '
      endif  else begin 

         ch_n_contr(iobs) = 0
      endelse
      
   endfor
endif

;Print the results, store them in the ouput file, and also create an
;-------------------------------------------------------------------
; array with the (eventual) observed wavelengths of the lines that
;------------------------------------------------------------------
;do not have theoretical counterparts. 
;-------------------------------------

i_ex=0
for iobs=0,n_obs-1 do begin
   
   print, obs_id(iobs),obs_wvl(iobs),obs_int(iobs),obs_sig(iobs)
   n_list=ch_n_contr(iobs)
   if n_list ge 1 then begin
      for ilist=0,n_list-1 do begin
         print,ch_wvl(ilist,iobs),ch_id(ilist,iobs),$
               max(ch_contr_wa(*,ilist,iobs)), $
               format='(10x,f10.3,a15,e10.2)'

         format='(f10.3,1x,f10.3,a15,2i5,1x,'+trim(n_dem_temp)+'e11.3,a40)'

         printf,lucontr,obs_wvl(iobs),ch_wvl(ilist,iobs),ch_id(ilist,iobs),$
                ch_z(ilist,iobs),ch_ion(ilist,iobs),ch_contr_wa(*,ilist,iobs),$
                ch_term(ilist,iobs), $
                format=format

      endfor

   endif else begin

      print,'NO CHIANTI LINES FOUND corresponding to ',obs_wvl(iobs),' A'

      exclu_obs_wvl_no_teo(i_ex)=obs_wvl(iobs)
      i_ex=i_ex+1  

   endelse
endfor

free_lun,lucontr

;clean up the zeros:
;-------------------
n=where(exclu_obs_wvl_no_teo ne 0.)
if n(0) eq -1 then exclu_obs_wvl_no_teo=0. else exclu_obs_wvl_no_teo=exclu_obs_wvl_no_teo(n)

;now exclu_obs_wvl_no_teo has the list of wavelength of the lines that do not
;have theoretical counterparts, or is 0.


end
;
;-----------------------------------------------------
