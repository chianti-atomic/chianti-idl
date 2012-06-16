;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
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
;		Then calculates the G(T) values for the temperature interval
;		log(T)= 4.0 - 8.0  with steps of log(T) = 0.1
;		A constant pressure OR a constant density for all the lines
; 		is used. If you select a constant pressure,
;               for each ion the contribution function is calculated at an 
;               electron density N_e equal to the ratio of the pressure 
;               and the temperature of maximum ionization fraction:  
;               C=C( T, N_e = P/T_ion)
;               The C(T) values are stored by GET_CONTRIBUTIONS in the output 
;		file OUTPUT.CONTRIBUTIONS.
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
;		In the unlikely event that more than 20 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;
;	PRESSURE:     		the value of the pressure (Ne T)
;
;	OUTPUT  :  	-the core name for the output
;
;
; Calls       : 
;		READ_IONEQ 
;               read_masterlist
;		CONVERTNAME
;               ion2spectroscopic
;		ZION2FILENAME
;		READ_WGFA2
;		READ_SPLUPS
;		POP_SOLVER
;		READ_ELVLC
;		READ_IONREC
;		CONVERT_TERMS
;
; Common      : elvlc                          - energy levels
;               wgfa                           - radiative data
;               upsilon                        - upsilon data
;		      
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
;		In the unlikely event that more than 20 (default value for 
;		N_MATCHES) theoretical lines corresponding to an observed
;		line are found, the routine stops; in this case, you have to 
;		start again setting N_MATCHES equal to a greater number. 
;
;		Of course you need to have the enviroment variable CDS_SS_DERE
;		pointing to the CHIANTI database top directory.
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
; VERSION     :    10, 12-Jun-2009
;
;
;-        
;---------------------------

pro get_contributions,output=output,density=density, pressure=pressure,$
       cut_contrib=cut_contrib,n_ch=n_ch

;pick up common definitions:
;---------------------------
COMMON wgfa, wvl,gf,a_value
COMMON upsilon, splstr
COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref


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

;common with pop_solver:

COMMON proton, pstr, pe_ratio


ON_ERROR, 2

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

;
; we assume as default to include the protons.
;
; The proton-to-electron ratio is calculated for all temperatures, using 
; the user-specified ion balance and the default abundance file 
; (!abund_file)

pe_rat_all=proton_dens(ioneq_t)


dlnt=ALOG(10.^(ioneq_t[1]-ioneq_t[0]))      
n_ioneq_t=n_elements(ioneq_t)


;
;pick up the ion masterlist:
;---------------------------

 mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')

read_masterlist,mname,list_ions
nlist=n_elements(list_ions)

;
;   main input and calculation loop  **************
;
print,' getting CHIANTI contribution functions G(T)  for ALL the  observed lines'


FOR  ilist_ion=0,nlist-1 DO  BEGIN 

   gname=list_ions(ilist_ion)

   convertname,gname,iz,ion
   ion2spectroscopic,gname,snote

   locname=strlowcase(gname)
   pos=strpos(locname,'_')
   l=strlen(pos)
   first=strmid(locname,0,pos)
   last=strmid(locname,pos+1,l-pos-1)
;
   if strpos(last,'d') ge 0 then dielectronic=1 else dielectronic=0
;


;convert z and ionisation stage to filename 
;(eg z=26, ion=24 > !xuvtop/fe/fe_24 ) :
;-------------------------------------------

   zion2filename,iz,ion,fname,diel=dielectronic

   wname=fname+'.wgfa'
   elvlcname=fname+'.elvlc'
   upsname=fname+'.splups'
   pname=fname+'.psplups'

 print,'Calculating ...  ', snote


;
; first check if ion fraction  exists 
;

   this_ioneq=ioneq[*,iz-1,ion-1+dielectronic]

; 't_index'  contains the indices of 
; ioneq_t where the   ion fraction is not zero.

   t_index=WHERE(this_ioneq NE 0.)

   IF t_index[0] NE -1 THEN BEGIN

; overlap exists, so now read .wgfa file to see if there are any 
; lines 
;
;  read in level information, wavelengths, gf and A values from .wgfa files:
;---------------------------------------------------------------------------

      read_wgfa2,wname,lvl1,lvl2,wvl1,gf1,a_value1,wgfaref


;DO NOT GET THE LINES WITH UNOBSERVED ENERGY VALUES FOR NOW.


;      check to see if there any wavelength matches for this ion.
; loop over the obs lines.
;---------------------------------------------------------------
      anylines = -1

      FOR iobs=0, n_obs-1 do BEGIN

         dummy=where((abs(obs_wvl(iobs) - wvl1) le obs_delta_lambda(iobs) ) $
                     and (a_value1 NE  0.), nn)
;
         if nn  GE  n_ch then begin
            message,'Error,  increase n_matches to at least '+string(nn)
         ENDIF
         anylines =  anylines > dummy
      ENDFOR
      

      IF  max(anylines) GE  0 THEN  BEGIN 
;
;        there is a wavelength match for this ion so calculate populations
;        -----------------------------------------------------------------
;

         ntrans=n_elements(lvl1)
         nlvls=max([lvl1,lvl2])

; For solving the level balance eqns., need wvl, gf and a_value in 2-D 
; arrays.
;------------------------------

; It's possible that the same transition can have two A-values, 
; one being the standard radiative decay, the other being an 
; autoionisation prob. They need to be added together when 
; constructing the a_value 2-D array.

         ind1 = where(wvl1 EQ 0.)
         ind2 = where(wvl1 NE 0.)

         wvl=FLTARR(nlvls,nlvls)
         gf=FLTARR(nlvls,nlvls)
         a_value=FLTARR(nlvls,nlvls)

         wvl[lvl1-1,lvl2-1]= abs(wvl1)
         gf[lvl1-1,lvl2-1]=gf1

         IF ind1[0] NE -1 THEN BEGIN
            a_value[lvl1[ind1]-1,lvl2[ind1]-1] = a_value1[ind1]
            a_value[lvl1[ind2]-1,lvl2[ind2]-1] = $
              a_value[lvl1[ind2]-1,lvl2[ind2]-1] + a_value1[ind2]
         ENDIF ELSE BEGIN
            a_value[lvl1[ind2]-1,lvl2[ind2]-1] = a_value1[ind2]
         ENDELSE


         read_elvlc,elvlcname,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref

         g=where(ecm EQ 0.)
         IF max(g) GT 0 THEN ecm(g)=ecmth(g)
         mult=2.*jj+1.

;read e collisional data:

      read_splups,upsname,splstr,splref

;read p collisional data:

read_splups, pname, pstr, pref, /prot

;read ion/rec data:

read_ionrec,fname,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status
         
         t_ioneq=ioneq_t             ; Sets the temperature array for ion abundances in variabile ionrec
         IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
         IF status gt 0 THEN BEGIN
              IF ion gt 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
              IF ion eq 1 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-1:ion))  ; No coll.ionization to neutral ions!
              ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
              lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:t_ioneq}
         ENDIF


;  calculate level populations
;------------------------------

         this_ioneq=ioneq[*,iz-1,ion-1+dielectronic]

         temp=10.^ioneq_t[t_index]

         nt = n_elements(t_index)

         IF keyword_set(density) THEN BEGIN 
            dens = fltarr(N_ELEMENTS(temp))
            dens(*)=density 
         ENDIF ELSE  dens=pressure/temp

         pops=DBLARR(nt,nlvls)

         FOR i=0,nt-1 DO BEGIN

            getmin=min(abs(ioneq_t-alog10(temp[i])),i_t)
            pe_ratio=pe_rat_all[i_t]

               input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
                  wvl:wvl, a_value:a_value, splstr:splstr, $
                  pe_ratio:pe_ratio,prot_struc:pstr, ionrec:ionrec}
;, dilute:dilute, radtemp:radt}

         pop_solver,input, temp[i],dens[i],pop

         pop = reform(pop)
         np = N_ELEMENTS(pop)
         pops[i,0:np-1]=pop

         ENDFOR


         tot_n_lines_this_ion = 0

         FOR  iobs=0,n_obs-1 DO  BEGIN 
;   
;exclude the lines that do not have observed energies:
;(otherwise an  ABS(wvl1) should be placed)

            anylines=where((abs(obs_wvl(iobs) - wvl1) le obs_delta_lambda(iobs) ) $
                           and (a_value1 NE  0.), nn)

            IF  nn GT  0 THEN  BEGIN

;;;;;
;print, obs_wvl(iobs), nn

               n_lines_this_iobs = 0

               FOR i=0,nn-1 DO BEGIN ; for each contributing line within iobs

; The check on pops in the line below makes sure that no lines 
; of 0 intensity get into the line list. The check on splups 
; is necessary to stop lines that have already appeared in the 
; standard CHIANTI file turning up in the dielectronic ion 
; spectrum as well.

         IF TOTAL(pops[*,lvl2[anylines[i]]-1]) NE 0. THEN BEGIN
            siz=max([splstr.lvl1,splstr.lvl2])

            IF  lvl2[anylines[i]] LE siz THEN  BEGIN 

;
; With the change to storing the collision data in a structure (splstr), 
; the condition below had to be changed.
;
               i_spl=where(splstr.lvl2 EQ lvl2[anylines[i]])
               IF i_spl[0] NE -1 THEN BEGIN

;
;        there is a wavelength match for this ion and this observed wavelength
;

                           if ch_n_contr(iobs) gt n_ch-1 THEN $
                             message, ' Error, increase n_matches at least to '+string(ch_n_contr(iobs)+1)


                           tot_n_lines_this_ion = tot_n_lines_this_ion + 1

                           n_lines_this_iobs =n_lines_this_iobs +1

;print, snote, wvl1[anylines[i]]
;
                           ch_wvl(ch_n_contr(iobs),iobs)=wvl1[anylines[i]]
                           ch_id(ch_n_contr(iobs),iobs)= snote
                           ch_z(ch_n_contr(iobs),iobs)=iz
                           ch_ion(ch_n_contr(iobs),iobs)=ion


                  l1 = lvl1[anylines[i]]
                  l2 = lvl2[anylines[i]]

                  this_design =  CONVERT_TERMS(l1, l2, result_latex =result_latex)

                ch_term(ch_n_contr(iobs),iobs)=this_design

;                           ch_term(ch_n_contr(iobs),iobs)=strtrim(term(lvl1(anylines[i])-1),2)+$
;                             ' - '+strtrim(term(lvl2(anylines[i])-1),2)


;             hc=6.626d-27*2.998d+10*1.d+8/(4.d*3.14159d*abs(wvl1(anylines[i])))
;           IF KEYWORD_SET(photons) THEN de = 1. ELSE $
                           de = 1.986e-8/ABS(wvl1[anylines[i]])

                  de_nj_Aji=pops[*,lvl2[anylines[i]]-1]*de*a_value1[anylines[i]]

                        tgt = (de_nj_Aji)*this_ioneq[t_index]/dens
                        tgt = tgt/4./!pi
                        get_tmax = MAX(tgt,tmax_ind)
                        tmax = ALOG10(temp[tmax_ind])

                        this_goft = dblarr(n_ioneq_t)

                        this_goft[t_index] = tgt


                           ch_contr_wa(*,ch_n_contr(iobs),iobs)=this_goft

; ch_n_contr is an  intarr(n_obs) defined within CHIANTI_DEM:

                           ch_n_contr(iobs)=ch_n_contr(iobs)+1

                        ENDIF 
                     ENDIF  
                  ENDIF 

               ENDFOR           ; for each contributing line within iobs

; ch_contr_list is an   intarr(n_ch,n_obs) defined within CHIANTI_DEM, but not needed.

            ENDIF               ;nn gt 0
                                ;there were contributing lines  within iobs (else ch_n_contr(iobs)=0)

         ENDFOR                 ;  iobs

         IF tot_n_lines_this_ion GT 0 THEN $
           print,  'Found '+string(tot_n_lines_this_ion)+' lines' ELSE $
           print, 'No lines found'

      ENDIF ELSE  print, 'No lines found'
      ;;        there were  wavelength matches for this ion

   ENDIF  ELSE print, 'No ionization equilibrium values found for this ion...'


;; there were (as they should) ioneq_t values 
   ;;where the   ion fraction is not zero.

ENDFOR                          ; to do this ion from the masterlist
;

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

      print,'NO THEORETHICAL LINES FOUND corresponding to ',obs_wvl(iobs),' A'

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
