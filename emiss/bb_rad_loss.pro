
PRO bb_rad_loss,t,loss_rate,pressure=pressure,density=density, $
                sngl_ion=sngl_ion,no_setup=no_setup, noprot=noprot, $
                radtemp=radtemp, rphot=rphot, ioneq_file=ioneq_file, $
                abund_file=abund_file, element=element, min_abund=min_abund, $
                quiet=quiet 


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
;       Element:     Either the atomic number of an element, or the
;                    symbol name (e.g., 'Fe' for iron). The loss rate
;                    will only be computed for this element.
;       Min_Abund:   A minimum abundance below which elements will not
;                    be included.
;	
; KEYWORD PARAMETERS:
;       NOPROT:   Switches off inclusion of proton rates.
;       NO_SETUP: If set, then the abundance and ioneq common blocks
;                 are assumed to be filled from a previous call. This
;                 keyword is now *obsolete* and no longer does anything
;                 but is retained for backwards compatibility.
;       QUIET:    If set, then no text is printed to the screen.
;
; OUTPUTS:
;       Temperature:  array of temperatures (K)
;       Loss_rate:  energy loss rate in erg cm^3 s^-1
;
; CALLS
;       READ_MASTERLIST, CH_SETUP_ION
;       POP_SOLVER, CONVERTNAME, R2W,
;       ZION2SPECTROSCOPIC, ZION2FILENAME, READ_ABUND, READ_IONEQ
;
; EXAMPLE:
;      IDL> rad_loss, temp, rad
;      IDL> rad_loss, temp, rad, pressure=1e15
;      IDL> rad_loss, temp, rad, density=1e9
;      IDL> rad_loss, temp, rad, abund_file=!abund_file
;      IDL> rad_loss, temp, rad, sngl_ion=['o_6','fe_13']
;      IDL> rad_loss, temp, rad, element='fe'
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
;       v.9, 25-Jan-2018, Peter Young
;            Updated to use ch_setup_ion; removed common blocks; added
;            ELEMENT and MIN_ABUND optional inputs; added keyword
;            QUIET.
;
;       v.10, 26-Jan-2018, Peter Young
;            I forgot to implement rphot and radtemp through
;            ch_setup_ion, so this has been fixed now.
;
; VERSION     :   10, 26-Jan-2018, Peter Young
;-


IF n_params() LT 2 THEN BEGIN 
  print,'Use:  IDL> bb_rad_loss,temperature,loss_rate [,pressure= , density= , sngl_ion= ,'
  print,'           radtemp= , rphot= , /noprot, abund_file=, ioneq_file=, min_abund=, '
  print,'           element= ]'
  print,'Notes:'
  print,'    temperature is an output'
  print,'    sngl_ion allows individual ions to be computed (can be an array)'
  print,"    element allows an individual element to be computed (e.g., 'Fe' or 26)"
  return
ENDIF


default_density=1.e+10
IF (keyword_set(pressure) ne 1) and (keyword_set(density) ne 1) AND NOT keyword_set(quiet) THEN BEGIN
  print,'% BB_RAD_LOSS: keyword DENSITY not specified. Using ',default_density,' cm^-3'
ENDIF 

IF n_elements(min_abund) EQ 0 THEN min_abund=0.


;
; Read abundances. If abund_file has not been set, then just use
; CHIANTI default.
;
IF n_elements(abund_file) EQ 0 THEN BEGIN
  abund_file=!abund_file
  IF NOT keyword_set(quiet) THEN BEGIN 
    print,'% BB_RAD_LOSS: abund_file was not specified so using CHIANTI default (!abund_file).'
    print,'               ',abund_file
  ENDIF 
ENDIF 
read_abund,abund_file,abund,abund_ref


;
; Read ioneq file.
;
IF n_elements(ioneq_file) EQ 0 THEN ioneq_file=!ioneq_file
read_ioneq,ioneq_file,ioneq_t,ioneq,ioneq_ref
;
n_ioneq_t=n_elements(ioneq_t)
t=10.^ioneq_t
loss_rate=fltarr(n_ioneq_t)



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
; Handle the input 'element'.
;
elt_iz=-1
IF n_elements(element) NE 0 THEN BEGIN
  IF datatype(element) EQ 'STR' THEN BEGIN
    z2element,indgen(30)+1,elt,/symbol,/lower_CASE
    k=where(strlowcase(element) EQ elt,nk)
    IF nk NE 0 THEN elt_iz=k[0]+1
  ENDIF ELSE BEGIN 
    elt_iz=element[0]
  ENDELSE
ENDIF 


FOR  ilist=0,nlist-1 DO  BEGIN 
  gname=list[ilist]
  convertname,gname,iz,ion
  ion2spectroscopic,gname,snote, dielectronic=dielectronic

  IF elt_iz NE -1 THEN test1=elt_iz EQ iz ELSE test1=1b
  test2=(abund(iz-1) GT min_abund)
 ;
  IF test2 AND test1 THEN BEGIN
   ;
    IF NOT keyword_set(quiet) THEN print,snote
    input=ch_setup_ion(gname,noprot=noprot,ioneq_file=ioneq_file,abund_file=abund_file, $
                      radtemp=radtemp,rphot=rphot)
   ;
   ; Filter out zero wavelength transitions (2-photon; autoionization)
   ;
    wgfa=input.wgfastr
    k=where(wgfa.aval NE 0.,ntrans)
    lvl2=wgfa[k].lvl2
    a_value1=wgfa[k].aval
    wvl1=abs(wgfa[k].wvl)
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
