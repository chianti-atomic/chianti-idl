;+
;
; PROJECT:  CHIANTI
;
; NAME:
;	bb_rad_loss
;
; PURPOSE:
;
;       calculates energy loss rate by line (bound-bound) radiation
;
; CATEGORY:
;	
;	synthetic spectra
;
; CALLING SEQUENCE:
;
;       BB_RAD_LOSS,Temperature,Loss_rate
;
;
; INPUTS:
;
;	None:  the user will select various parameters such as the 
;              choice of elemental abundances and ionization equilibria
;
;  KEYWORDS:
;
;       Pressure:  pressure in emitting region (cm^-3 K)
;                  density=pressure/temperature(K)
;
;       Density:   density (cm^-3), constant for all temperatures
;                  if neither density or pressure is set, then a 
;                  default constant density of 10^10 cm^-3 is used.
;
;	Sngl_ion:  to calculate the loses spectrum for a single ion
;
;       NOPROT    Switches off inclusion of proton rates.
;
;       RADTEMP   The blackbody radiation field temperature (default 
;                 6000 K).
;
;       RPHOT     Distance from the centre of the star in stellar radius 
;                 units. I.e., RPHOT=1 corresponds to the star's surface. 
;                 (Default is infinity, i.e., no photoexcitation.)
;
; OUTPUTS:
;
;       Temperature:  array of temperatures (K)
;       Loss_rate:  energy loss rate in erg cm^3 s^-1
;
;
; PROCEDURE:
;
;
;  if keyword pressure is set then calculations performed at constant pressure
;  if keyword density is set then calculations performed at constant density
;  otherwise, density is set to 1.e+10
;  
;  pressure = density * temperature  (cm^-3 K)
;
;
;	the user will be asked to select an abundance file. 
;
;
; COMMON BLOCKS:
;
;	ELVLC, WGFA, UPSILON, RADIATIVE, PROTON, ELEMENTS, IONREC
;
; CALLS
;
;       READ_MASTERLIST, 
;       SETUP_ELEMENTS, POP_SOLVER, CONVERTNAME, R2W, PROTON_DENS,
;       ZION2SPECTROSCOPIC, ZION2FILENAME, SETUP_ION
;
; EXAMPLE:
;
;       > bb_rad_loss,t,r
;
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	January 1999:  version 1, adopted from synthetic.pro
;       14-Jul-2000     Peter Young, now calls pop_solver
;
;       Ver.2, 9-Dec-2001, Peter Young
;               Modified for v.4 of CHIANTI.
;
;       V.3, 4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure. Also shorten the whole procedure, by calling
;                  other procedures.
;
;       V.4, 5-Apr-2005, EL
;                  modified the input to POP_SOLVER in order to add ion/rec
;                  in the calculation of level populations. All the rest
;                  is unchanged.
;
;       V 5, 25-May-2005, GDZ 
;                  corrected routine header.
;
;       V.6, 12-Jun-2009, Enrico Landi
;               Changed the definition of the temperature array for ion fractions
;               in the IONREC variable, now taken directly from the output of
;               READ_IONEQ.PRO
;
;       V.7, 25-Jun-2012, Peter Young
;               Changed list() to list[] for compatibility with IDL 8.
;
;
; VERSION     :   7, 25-Jun-2012, Peter Young
;
;-
pro bb_rad_loss,t,loss_rate,pressure=pressure,density=density, $
      sngl_ion=sngl_ion,no_setup=no_setup, noprot=noprot, $
      radtemp=radtemp, rphot=rphot
;
;
;  pressure= electron pressure (Ne x T)  cm-3 K
;
;
;
common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
common wgfa, wvl,gf,a_value
common upsilon, splstr
common elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref
COMMON proton, pstr, pe_ratio
COMMON radiative, radt, dilute
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

;
if n_params() lt 2 then begin
   print,' '
   print,' type> bb_rad_loss,temperature,loss_rate [,pressure= , density= , sngl_ion= ,'
   print,'         /no_setup, radtemp= , rphot= , /noprot ]'
   print,' '
   return
endif
;
;
default_density=1.e+10
if (keyword_set(pressure) ne 1) and (keyword_set(density) ne 1) then begin
   print,' using constant default density = ',default_density
endif
;
;
;   read elemental abundances, ionization equilibrium 
;
if not keyword_set(no_setup) then begin
   setup_elements
;   read_ip,!xuvtop+'/ip/chianti.ip',ionpot,ipref
endif
;
n_ioneq_t=n_elements(ioneq_t)
t=10.^ioneq_t
loss_rate=fltarr(n_ioneq_t)

;
; get proton/electron ratio for all temperatures
;
pe_rat_all=proton_dens(alog10(t))


IF n_elements(rphot) EQ 0 THEN dilute=0d0 ELSE dilute=r2w(rphot)
IF n_elements(radtemp) EQ 0 THEN radt=6d3 ELSE radt=double(radtemp)


;
;  open the file that has the names of the ions
;
;
mname=!xuvtop+'/masterlist/masterlist.ions'
;
;
if keyword_set(sngl_ion) then begin
   list=sngl_ion
endif else begin
   read_masterlist,mname,list
endelse
;
nlist=n_elements(list)
;
;
;   main input and calculation loop  **************
;
;
FOR  ilist=0,nlist-1 DO  BEGIN 

;GDZ-

   gname=list[ilist]
   convertname,gname,iz,ion
   ion2spectroscopic,gname,snote, dielectronic=dielectronic

;   test1=(not dielectronic)
   test2=(abund(iz-1) gt 0.)
;
   if test2 then begin
;
      print,snote
;

      setup_ion,gname,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1, $
        noprot=noprot           ;path=path

; include all lines
;
      wvl1=abs(wvl1)

      ntrans=n_elements(lvl1)

;
;  calculate level populations
;
      this_ioneq=ioneq(*,iz-1,ion-1+dielectronic)
;
;        edensity=pressure/tmax
;
;
      gioneq_t=where(this_ioneq gt 0.,ngt)
      
;
      for ig=0,ngt-1 do begin
         it=gioneq_t(ig)
         tt=t(it)
;
         if keyword_set(pressure) then begin
            edensity=pressure/tt
         endif else if keyword_set(density) then begin
            edensity=density
         endif else edensity=default_density
;
         pe_ratio=pe_rat_all[it]

;EL

         t_ioneq=ioneq_t           ; Sets the temperature array for ion abundance for ionrec structure
 
         IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
         IF status gt 0 THEN BEGIN
              IF ion ge 2 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
              IF ion lt 2 THEN ioneq_ionrec=reform(ioneq(*,iz-1,0:ion))
              ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
                lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:t_ioneq}
         ENDIF


;GDZ

         input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
                  wvl:wvl, a_value:a_value, splstr:splstr, $
                  pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, radtemp:radt, ionrec:ionrec}
;                  pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, radtemp:radt}

         pop_solver,input, tt,edensity,pop

         pop = reform(pop)
         npop=n_elements(pop)
;
         for itrans=0L,ntrans-1 do begin
            if (lvl2(itrans) le npop) and (wvl1(itrans) gt 0.) then begin
               hc=6.626d-27*2.998d+10*1.d+8*pop(lvl2(itrans)-1)*a_value1(itrans)/wvl1(itrans)
               loss_rate(it)=loss_rate(it)+hc*abund(iz-1)*this_ioneq(it)/edensity
            endif
         endfor
      endfor                    ;  ig
;
;printf,luo,lvl1,lvl2,ww,gg,aa,format='$(2i5,f15.5,2e15.3)'
;
   endif                        ;  if block not  dielectronic lines
;
endfor                 ;end of main input and calculation loop   ***************
;
;
end



