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
;
;		Is the file where all the contribution  functions G(T) are 
;		stored. In the first two lines  the ionization equilibrium 
;		file name, and the constant value of pressure or density 
;		adopted are reported. Then for each line you have reported  
;               the observed wavelength, the theoretical one, the element and
;		ionization stage, then the C(T) values. At the end the 
;		specification for each transition.
;
;               OUT
;               An IDL structure with the data, to be passed back to
;               the main routine.
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
;	N_CH:              
;		In the unlikely event that more than 100 (default value for 
;		N_CH) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_CH equal to a greater number. 
;
;	PRESSURE:     	the value of the pressure (Ne T)
;
;	OUTPUT_NAME:  	the core name for the output
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
;       ADVANCED_MODEL: include density-dependent and CT effects.
;
;       CT: include charge transfer in advanced models
;
;       IONEQ_LOGT: an array of log T [K] values, defining the grid for the
;                   calculation
;
;       ATMOSPHERE: A file with the H,He abundances as a function of temperature.
;                      By default, the file avrett_atmos.dat is read, with data from
;                      Avrett E.H., Loeser R., 2008, ApJ, 175, 229
;
;       HE_ABUND:  The total helium abundance relative to hydrogen. 
;
;       NO_AUTO: If set, then the autoionization rates (contained in
;                the .auto file) are not read. The autoionization states are not
;           included in the calculations, i.e. a single ion rather than the
;           two-ion model  introduced in version 9 is calculated. This speeds
;           up the calculations without affecting the lines from the bound states.
;
;       DR_SUPPRESSION: Switch on DR suppression from Nikolic et al (2018) for all ions 
;              not included in the advanced models. The comparison with Summers (1974) suppression
;              has not been checked for other elements when preparing the models.
;
; Calls       : 
;		ch_synthetic
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
;       V.13,  8-Oct-2020, GDZ
;                  Major revision, removed COMMON blocks and renamed
;                  some variables.
;                  Also added keywords to include photoexcitation in
;                  the input, rphot, radtemp,RADFUNC
;
;       V.14, 4 Nov 2023, GDZ
;           Major rewrite, adding the advanced model option (the default),
;           and also NO_AUTO: If set, then the autoionization states are not
;           included in the calculations, i.e. a single ion rather than the
;           two-ion model I introduced in version 9 is calculated. This speeds
;           up the calculations without affecting the lines from the bound states.
;           Added passing ioneq_name
;
;       v.15, 1-Jul-2024, GDZ, added  dr_suppression
;
; VERSION     :    V.15
;
;
;-        
;---------------------------

pro get_contributions,input, out,  output_name=output_name,density=density, pressure=pressure,$
                      cut_contrib=cut_contrib,n_ch=n_ch,$
                      rphot=rphot, radtemp=radtemp,RADFUNC=RADFUNC, ioneq_name=ioneq_name,$
                no_auto=no_auto,ioneq_logt=ioneq_logt, advanced_model=advanced_model,ct=ct,$
                atmosphere=atmosphere,he_abund=he_abund,dr_suppression=dr_suppression


;ON_ERROR, 2

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
   message, 'system variable !xuvtop must be set  '
xuvtop = !xuvtop

obs_int=input.obs_int
obs_sig= input.obs_sig
obs_wvl= input.obs_wvl
obs_id = input.obs_id
obs_delta_lambda= input.obs_delta_lambda

n_obs=n_elements(obs_int)


;create an array  where eventual observed wavelength that DO NOT have 
;------------------------------------------------------------
;corresponding theoretical lines (at maximum can have an n_obs length)
;------------------------------------------------------------
;will be stored.
;---------------

exclu_obs_wvl_no_teo=fltarr(n_obs)


; inizialise the arrays:                
;----------------------

ch_wvl=fltarr(n_ch,n_obs)
ch_term=strarr(n_ch,n_obs)

ch_z=intarr(n_ch,n_obs)
ch_ion=intarr(n_ch,n_obs)
ch_id=strarr(n_ch,n_obs)


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

print,'% GET_CONTRIBUTIONS:  getting contribution functions G(T) for ALL the observed lines'

; if not set, and if the minimum wavelength is above 50 Angstroms, set NO_AUTO=1

if not keyword_set(NO_AUTO) and w1 gt 50 then begin
print, '% GET_CONTRIBUTIONS: automatically setting NO_AUTO=1 to speed up the calculations '
   NO_AUTO=1
end

if  keyword_set(NO_AUTO) and w1 lt 50 then begin
print, '% GET_CONTRIBUTIONS WARNING: you have set NO_AUTO=1 so no satellite lines will be calculated! '
end

   
ch_synthetic, w1, w2, output=output_ch, err_msg=err_msg, msg=msg, $
              pressure=pressure, density=density, $
              /all, ioneq_name=ioneq_name, $
              noprot=noprot, rphot=rphot, radtemp=radtemp,RADFUNC=RADFUNC, /goft,$
              no_auto=no_auto,ioneq_logt=ioneq_logt, advanced_model=advanced_model,ct=ct,$
              atmosphere=atmosphere,he_abund=he_abund,dr_suppression=dr_suppression


n_dem_temp=n_elements(ioneq_logt)
;; ;this array will have the C(T) 
;; ;-----------------------------
ch_contr_wa=fltarr(n_dem_temp,n_ch,n_obs)


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

savegen, file=output_name+'_gt.genx', st=out
delvarx,out
print, 'savefile '+output_name+'_gt.genx  saved ! '


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;

output_contributions=output_name+'.contributions'
openw,lucontr,output_contributions,/get_lun

; v.11: do not print the ioneq name, but rather the grid of log T values:
; printf,lucontr,ioneq_name
printf, lucontr, arr2str(ioneq_logt,' ',/trim)


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

; get the output data:

out={ioneq_name:ioneq_name, ioneq_logt:ioneq_logt,const_net:const_net,$ 
  exclu_obs_wvl_no_teo:exclu_obs_wvl_no_teo,$
 ch_wvl:ch_wvl,ch_id:ch_id,ch_z:ch_z,ch_ion:ch_ion,ch_contr_wa:ch_contr_wa,$
 ch_term:ch_term,ch_n_contr:ch_n_contr}


end


;
;-----------------------------------------------------
