
PRO bb_rad_loss,t,loss_rate,pressure=pressure,density=density, $
      sngl_ion=sngl_ion,no_setup=no_setup, noprot=noprot, $
      radtemp=radtemp, rphot=rphot, ioneq_file=ioneq_file, abund_file=abund_file 


;+
; NAME:
;      BB_RAD_LOSS
;
; PURPOSE:
;      Compute the bound-bound radiative losses from a plasma.
;
; CATEGORY:
;      CHIANTI; radiative loss.
;
; CALLING SEQUENCE:
;      BB_RAD_LOSS, T, Loss_Rate
;
; INPUTS:
;      An element abundance file is chose through a widget input
;      (unless abund_file is specified).
;
; OPTIONAL INPUTS:
;       Ioneq_File:  The name of an ionization equilibrium file. If
;                    not set, then the default CHIANTI file is used.
;       Abund_File:  The name of an element abundance file. If not
;                    set, then the user manually selects the file
;                    using a widget.
;       Pressure:    Pressure in emitting region (cm^-3 K)
;                    density=pressure/temperature(K)
;       Density:     Density (cm^-3), constant for all temperatures.
;                    If neither density or pressure is set, then a 
;                    default constant density of 10^10 cm^-3 is used.
;	Sngl_ion:    To calculate the loss spectrum for a single
;                    ion. Ions specified in CHIANTI format (e.g.,
;                    'o_6'). Can be an array of multiple ions.
;       Radtemp:     The blackbody radiation field temperature (default 
;                    6000 K).
;       Rphot:       Distance from the centre of the star in stellar radius 
;                    units. I.e., RPHOT=1 corresponds to the star's surface. 
;                    (Default is infinity, i.e., no photoexcitation.)
;	
; KEYWORD PARAMETERS:
;       NOPROT:   Switches off inclusion of proton rates.
;       NO_SETUP: If set, then the abundance and ioneq common blocks
;                 are assumed to be filled from a previous call.
;
; OUTPUTS:
;       Temperature:  array of temperatures (K)
;       Loss_rate:  energy loss rate in erg cm^3 s^-1
;
; COMMON BLOCKS:
;      ELVLC, WGFA, UPSILON, ELEMENTS, PROTON, RADIATIVE, ION REC.
;
; CALLS
;
;       READ_MASTERLIST, 
;       SETUP_ELEMENTS, POP_SOLVER, CONVERTNAME, R2W, PROTON_DENS,
;       ZION2SPECTROSCOPIC, ZION2FILENAME, SETUP_ION
;
; EXAMPLE:
;      IDL> rad_loss, temp, rad
;      IDL> rad_loss, temp, rad, pressure=1e15
;      IDL> rad_loss, temp, rad, density=1e9
;      IDL> rad_loss, temp, rad, abund_file=!abund_file
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
;       v.8, 8-Aug-2017, Peter Young
;               Added abund_file and ioneq_file optional inputs;
;               updated header; speeded up by using array operations
;               (removed loop over temperatures).
;
; VERSION     :   8, 8-Aug-2017, Peter Young
;-


common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
common wgfa, wvl,gf,a_value
common upsilon, splstr
common elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref
COMMON proton, pstr, pe_ratio
COMMON radiative, radt, dilute
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status


if n_params() lt 2 then begin
   print,' '
   print,' type> bb_rad_loss,temperature,loss_rate [,pressure= , density= , sngl_ion= ,'
   print,'         /no_setup, radtemp= , rphot= , /noprot, abund_file=, ioneq_file= ]'
   print,' '
   return
endif


default_density=1.e+10
if (keyword_set(pressure) ne 1) and (keyword_set(density) ne 1) then begin
   print,' using constant default density = ',default_density
endif


;
; Read elemental abundances, ionization equilibrium.
;  PRY, 8-Aug-2017: if /no_setup is set but abund is empty then I use
;  setup_elements to prevent a crash.
;
chck=n_elements(abund)
if not keyword_set(no_setup) OR chck EQ 0 then begin
  IF keyword_set(no_setup) AND chck EQ 0 THEN BEGIN
    IF NOT keyword_set(quiet) THEN print,'%BB_RAD_LOSS: /no_setup was used, but common block was empty.'
  ENDIF 
   setup_elements, abund_file=abund_file, ioneq_file=ioneq_file
endif


n_ioneq_t=n_elements(ioneq_t)
t=10.^ioneq_t
loss_rate=fltarr(n_ioneq_t)

;
; Get proton/electron ratio for all temperatures.
; Note that ioneq and abund data are passed through the common block,
; so they don't need to be specified.
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

FOR  ilist=0,nlist-1 DO  BEGIN 
  gname=list[ilist]
  convertname,gname,iz,ion
  ion2spectroscopic,gname,snote, dielectronic=dielectronic

  test2=(abund(iz-1) gt 0.)
 ;
  IF test2 THEN BEGIN
   ;
    print,snote
    setup_ion,gname,-1,-1,wvltst,lvl1,lvl2,wvl1,gf1,a_value1, $
              noprot=noprot           ;path=path
    wvl1=abs(wvl1)
   ;
   ; Filter out zero wavelength transitions (2-photon; autoionization)
   ;
    k=where(wvl1 NE 0.,ntrans)
    lvl1=lvl1[k]
    lvl2=lvl2[k]
    wvl1=wvl1[k]
    a_value1=a_value1[k]
    
    ntrans=n_elements(lvl1)

   ;
   ; To help speed up operation I only use temperatures where ioneq is >
   ; 10^-8 of max value. I chose this value based on the losses for H I. 
    this_ioneq=ioneq(*,iz-1,ion-1+dielectronic)
    gioneq_t=where(this_ioneq gt max(this_ioneq)*1e-8,ngt)

    tt=t[gioneq_t]

    if keyword_set(pressure) then begin
      edensity=pressure/tt
    endif else if keyword_set(density) then begin
      edensity=density
    endif else edensity=default_density

    pe_ratio=pe_rat_all[gioneq_t]

    IF status lt 0 THEN ionrec = {rec:0., ci:0., temp:0., lev_up_rec:0., lev_up_ci:0., status:-1, ioneq:0., temp_ioneq:0.}
    IF status gt 0 THEN BEGIN
      IF ion ge 2 THEN ioneq_ionrec=reform(ioneq(*,iz-1,ion-2:ion))
      IF ion lt 2 THEN ioneq_ionrec=reform(ioneq(*,iz-1,0:ion))
      ionrec = {rec:rec_rate, ci:ci_rate, temp:temp_ionrec, lev_up_rec:luprec, $
                lev_up_ci:lupci, status:status, ioneq:ioneq_ionrec, temp_ioneq:ioneq_t}
    ENDIF

    input = {gname:gname, jj:jj, ecm:ecm,ecmth:ecmth, $
             wvl:wvl, a_value:a_value, splstr:splstr, $
             pe_ratio:pe_ratio,prot_struc:pstr, dilute:dilute, radtemp:radt, ionrec:ionrec}
    
    pop_solver,input, tt,edensity,pop, pressure=keyword_set(pressure)

    pop = reform(pop)
    npop=n_elements(pop)

    id_t=make_array(ngt,value=1.)
    id_pop=make_array(ntrans,value=1.)

   ;
   ; This is a 2D array of temperature * transitions (NGT*NTRANS)
   ;
    losses = 6.626d-27*2.998d10*1.d8 * $
             pop[*,lvl2-1]* $
             (id_t#(a_value1/wvl1)) * $
             abund[iz-1]* $
             ((this_ioneq[gioneq_t]/edensity)#id_pop)

   ;
   ; Sum over transitions and add to total loss rate
   ;
    ion_losses=total(losses,2)
    loss_rate[gioneq_t]=loss_rate[gioneq_t]+ion_losses
   ;
  ENDIF 
ENDFOR


END
