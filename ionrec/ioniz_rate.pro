
FUNCTION ioniz_rate,gname,temperature_in,z=z,ion=ion,verbose=verbose, $
                    di_rates=di_rates, ea_rates=ea_rates, data=data

;+
; NAME:
;     IONIZ_RATE
;
; PURPOSE:
;     This routine computes the total ionization rate coefficient of an ion.
;     For more details see CHIANTI Technical Report No. 34.
;
; CATEGORY:
;     CHIANTI; rates; ionization.
;
; CALLING SEQUENCE:
;     Result = IONIZ_RATE( Ion_Name, Temperature )
;
; INPUTS:
;     Gname:  A string, specifying the ion name in Chianti style. For example
;             'o_6' specifies O VI or O+5.
;     Temperature:  Temperature (K). Can be an array.
;
; OPTIONAL INPUTS:
;     Z:  An integer specifying the atomic number of the element. If specified
;         along with ION, then it over-rides GNAME for specifying the ion.
;     Ion:  An integer specifying the spectroscopic number of the element. For
;           example, ion=6 corresponds to VI. If specified along with Z, then
;           it over-rides GNAME for specifying the ion.
;	
; KEYWORD PARAMETERS:
;     VERBOSE:  If specified, then informational messages are printed to the
;               IDL window.
;
; OUTPUTS:
;     The total ionization rate coefficient in units cm^3 s^-1 calculated for
;     the specified temperatures.
;
; OPTIONAL OUTPUTS:
;     Di_Rates:  An array of same size as TEMPERATURE containing the direct
;                ionization rates.
;     Ea_Rates:  An array of same size as TEMPERATURE containing the 
;                excitation-autoionization rates.
;     Data:  An IDL structure containing various parameters from the
;            calculations. The tags are:
;             .ion_name  Name of the ion in CHIANTI format.
;             .temp   Temperature (K).
;             .tot_rate  Total ionization rate coefficient (cm^3 s^-1).
;             .di_rate   Direct ionization rate coefficient (cm^3 s^-1).
;             .ea_rate   EA ionization rate coefficient (cm^3 s^-1).
;             .energy   Energy (eV) at which DI cross section tabulated.
;             .di_cross_section  Direct ionization cross section (cm^2).
;             .ioniz_potential  Ionization potential (eV).
;
;            energy and di_cross_section are float arrays with 12 elements,
;            corresponding to the 12 Gauss-Laguerre points used for the
;            integration. See CHIANTI Technical Report 34 for more details.
;
; CALLS:
;     ZION2FILENAME, SCALE_BT, DESCALE_BT, READ_SPLUPS, DESCALE_ALL
;
; EXAMPLE:
;     IDL> temp=10.^(findgen(11)/10.+5.6)
;     IDL> r=ioniz_rate('fe_11',temp)
;
;     IDL> r=ioniz_rate('',temp,z=26,ion=11)
;
; MODIFICATION HISTORY:
;     Ver.1, 17-Nov-2006, Ken Dere
;     Ver.2, 24-Sep-2015, Peter Young
;       Completed pairs of quotes in string expressions to help with
;       viewing code in certain types of text editor. Otherwise no
;       change to code.
;     Ver.3, 18-Apr-2025, Peter Young
;       Updated header; added informational messages; tidied up code; added
;       ea_rates and di_rates optional outputs.
;     Ver.4, 22-Jul-2025, Peter Young
;       Added data= optional output structure; updated header.
;-



if n_params() lt 2 then begin
   print,' '
   print,'Use: IDL> rate=ioniz_rate(gname,temperature,[z=, ion=, di_rates=, ea_rates='
   print,'                            /verbose, data=data ]) '
   print,'    calculate the ionization rate coefficient in cm^3 s^-1 '
   print,'    as a function of temperature (K) '
   print,' '
   return,-1
endif

alpha=5.287d+13   ; not the fine structure  constant
kb=1.380626d-16   ; Boltzmann constant (erg K^-1)
kbev=8.617343d-5  ; Boltzmann constant (eV)
;
tev=kbev*temperature_in
temperature=temperature_in
ntemps=n_elements(temperature)
;
dirate=dblarr(ntemps)
earate=dblarr(ntemps)
totrate=dblarr(ntemps)
;
;   gauss laguerre coefficients for n=12
;
ngl=12
xgl=[0.115722117358021d,0.611757484515131d,1.512610269776419d,2.833751337743509d     $
    ,4.599227639418353d,6.844525453115181d,9.621316842456871d,13.006054993306350d    $
    ,17.116855187462260d,22.151090379396983d,28.487967250983992d,37.099121044466926]


wgl=[2.647313710554435d-01,3.777592758731382d-01,2.440820113198774d-01,9.044922221168074d-02  $
    ,2.010238115463406d-02,2.663973541865321d-03,2.032315926629993d-04,8.365055856819753d-06  $          
    ,1.668493876540914d-07,1.342391030515027d-09,3.061601635035012d-12,8.148077467426124d-16]
;
;
if not keyword_set(iz) and not keyword_set(ion) then convertname,gname,z,ion
zion2filename,z,ion,fname

difile=fname+'.diparams'
tst=findfile(difile)
IF tst(0) EQ '' THEN BEGIN
  message,/info,/cont,'The direct ionization file (diparams) does not exist. Returning...'
  return,-1
ENDIF 

openr,lur,difile,/get_lun
;
idum=1.  & odum=1.
i5=intarr(5)
f3=fltarr(3)
eastr=''
f1=1.
f11=1.
str1=''

;
; Read first line of diparams file.
;
readf,lur,i5,format='(5i5)'
kz=i5(0)
kon=i5(1)
nspl=i5(2)
nfac=i5(3)

IF keyword_set(verbose) THEN message,/info,/cont,'Number of DI components for '+gname+': '+trim(nfac)+'.'

IF kz NE z AND kon NE ion then BEGIN
  message,/info,/cont,'The ion specification in the diparams file does not match the specified ion. Returning...'
  return,-1
ENDIF 
   
ff1=fltarr(nspl+1)
ff2=fltarr(nspl+1)
x_spline=fltarr(nspl,nfac)
y_spline=fltarr(nspl,nfac)
ev1=fltarr(nfac)
bt_bethe=fltarr(nfac)
btf=fltarr(nfac)

;
; Read the individual "transitions" of the DI rate.
;
FOR ifac=0,nfac-1 DO BEGIN
  readf,lur,ff1
  btf(ifac)=ff1(0)
  x_spline(0,ifac)=ff1(1:*)
  readf,lur,ff2
  ev1(ifac)=ff2(0)
  y_spline(0,ifac)=ff2(1:*)*1.d-14
  bt_bethe(ifac)=y_spline(nspl-1)
ENDFOR


;
; Only if neaev is non-zero do we read the easplups file.
;
neaev=i5[4]
if neaev gt 0 then begin
  readf,lur,eastr  ; lur,str1   ;   ,format='(e12.3)'
  eastra=strsplit(eastr,' ',/extract)
  f1=float(eastra)
  ;
  ; read in EA collision strengths
  ;
  eaname=fname+'.easplups'
  read_splups,eaname,ea_splups,ea_splupsref
ENDIF ELSE BEGIN
  IF keyword_set(verbose) THEN message,/info,/cont,'There are no excitation-autoionization rates for '+gname+'.'
ENDELSE 
;
free_lun,lur

;
; Compute the direct ionization rate coefficients
;
dum=1.   ; dummy variable
FOR it=0,ntemps-1 DO BEGIN 
  ;
  FOR ifac=0,nfac-1 DO BEGIN
    x0=double(ev1(ifac)/tev(it))
    beta=sqrt(kb*temperature(it))
    ;
    ; Energies corresponding to the Legendre polynomial roots (xgl).
    ; It is necessary to add the ionization potential (ev1) as the
    ; cross-section energies correspond to incident electron energies.
    ;
    egl=ev1(ifac)+xgl*tev(it)
    ;
    ; Convert the Legendre polynomial energies to scaled energies (btgl)
    ;
    scale_bt,egl,dum,dum,btf(ifac),ev1(ifac),btgl,dum,dum
    ;
    ; Run a spline through the diparams spline points to determine the
    ; scaled cross-sections corresponding to the Legendre polynomial roots.
    ; 
    y2=spl_init(x_spline(*,ifac),y_spline(*,ifac))
    btcross=spl_interp(x_spline(*,ifac),y_spline(*,ifac),y2,btgl)
    ;
    ; Convert the scaled cross-sections to regular cross-sections.
    ;
    btcross=double(btcross)
    descale_bt,btgl,btcross,dum,btf(ifac),ev1(ifac),evd2,crossgl,odum
    ;
    ; Integrate the cross-section using the Gauss-Laguerre technique to
    ; yield the direct ionization rate coefficient (newcross)
    ;
    newcross=alpha*beta*exp(-x0)*(total(wgl*xgl*crossgl)+x0*total(wgl*crossgl))
    dirate(it)=dirate(it)+newcross
  ENDFOR 
ENDFOR 

;
; Add the excitation-autoionization rates (if available).
;
IF neaev GT 0 THEN BEGIN
  nups=n_elements(ea_splups.lvl1)
  if n_elements(f1) lt nups then f1=replicate(f1,nups)
  ;
  IF keyword_set(verbose) THEN message,/info,/cont,'Number of EA components for '+gname+': '+trim(nups)+'.'
  ;
  for iups=0,nups-1 do begin
    descale_all,temperature,ea_splups,iups,ups
    xt=tev/(ea_splups(iups).de*13.6056981)
    earate=earate+f1(iups)*8.63d-6*ups*exp(-1./xt)/(sqrt(temperature))
  endfor
  totrate=dirate + earate
endif else totrate=dirate

di_rates=dirate
IF n_elements(earate) EQ 0 THEN ea_rates=dblarr(ntemps) ELSE ea_rates=earate

;
; Define the DATA output structure.
;
data={ion_name: gname, $
      temp: temperature, $
      tot_rate: totrate, $
      di_rate: di_rates, $
      ea_rate: ea_rates, $
      energy: evd2, $
      di_cross_section: crossgl, $
      ioniz_potential: ev1[0]}

return,totrate

END
