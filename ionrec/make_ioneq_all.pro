;+
; PROJECT:
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;               Florence (Italy), and the University of Cambridge .
;
; NAME
;     MAKE_IONEQ_ALL
;
;  PURPOSE:
;     Calculates the ionization equilibrium for all ions with Z=1->30
;     using ionization and recombination rates stored in CHIANTI.
;
; INPUTS
;    Temperature:  The temperature (in K) for which the ion fraction
;                  curves are needed. Temperature (K).
;
; OUTPUTS
;    Creates a new ion fraction file called 'new.ioneq' in the
;    user's working directory. An alternative name can be
;    specified with OUTNAME=.  
;
; OPTIONAL INPUTS:
;    Outname:  A string specifying the name of the output ioneq
;              file. It is recommended to append ".ioneq" to the
;              name. 
;
; OPTIONAL OUTPUTS:
;    Ion_rate: The array of ionization rates used to compute the ion
;              fractions. 
;    Rec_rate: The array of recombination rates used to compute the ion
;              fractions. 
;
; KEYWORDS
;    None.
;
; CALLS:
;    ZION2NAME, RECOMB_RATE, IONIZ_RATE
;
; EXAMPLE:
;    IDL> ltemp=findgen(51)/10.+4.0
;    IDL> make_ioneq_all,10.^ltemp
;    
; MODIFICATION HISTORY:
;    Ver.1, 3-March-2009, Ken Dere
;    Ver.2, 19-May-2010, Ken Dere
;         correct various errors
;    Ver.3, 19-Jul-2016, Peter Young
;         modified implementation of OUTNAME; introduced ion_rate= and
;         rec_rate= outputs; plots are no longer made; added a warning
;         if ion fractions below 10^4 K are computed.
;
; VERSION     :  3, 19-Jul-2016
;
;-
pro make_ioneq_all,temp,outname=outname, ion_rate=ion_rate, rec_rate=rec_rate
;
if n_params() lt 1 then begin
   print,' > make_ioneq_all, temp, [outname = ]'
   print, ' temp can be a single temperature or an array'
   print, ' if the name of the output file is not specified by outname'
   print, ' a file new.ioneq will be created that can be read by '
   print, ' read_ioneq'
   return
endif
;
ioneqmin=1.e-20
;
;   get ionization file list
;
zlabl=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
       'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', $
       'Mn','Fe','Co','Ni','Cu','Zn']
;
;  make ioneq file
;
;
ntemp=n_elements(temp)
tformat='('+trim(ntemp)+'f6.2)'
iformat='('+trim(ntemp)+'e10.3)'
;
;home = getenv('HOME')
IF n_elements(outname) EQ 0 THEN BEGIN
  outname='new.ioneq'
ENDIF 
;; if not keyword_set(outname) then begin
;;     outname = 'new.ioneq'
;;     print,' creating output file = ',outname
;; endif
;name=concat_dir(home,outname)
openw,luw,outname,/get_lun
printf,luw,ntemp,30,format='(2i3)'
printf,luw,alog10(temp),format=tformat

FOR z=1,30 DO BEGIN 
  el=strlowcase(zlabl(z-1))
  ion_rate=fltarr(ntemp,z+2)
  rec_rate=fltarr(ntemp,z+2)
 ;
 ; Get ionization rates
 ;
  for ion=0,z-1 do begin
    zion2name,z,ion+1,gname
    ion_rate(0,ion) = ioniz_rate(gname,temp)
  endfor
 ;
;;   if ntemp gt 1 then begin
;;     yr=[1.e-10,1.]*max(ion_rate)
;;     plot_oo,fltarr(2),fltarr(2),xr=[1.e+3,1.e+9],yr=yr,/nodata,title=zlabl(z-1)
;; ;
;;     for ion=0,z-1 do begin
;;       oplot,temp,ion_rate(*,ion)
;;     endfor
;;   endif
 ;
 ; Get recombination rates
 ;
  for ion=1,z do begin
    zion2name,z,ion+1,gname
    rec_rate(0,ion)=recomb_rate(gname,temp)
  endfor
 ;
 ;
  ;; if ntemp gt 1 then begin
  ;;   for ion=1,z do begin
  ;;     oplot,temp,rec_rate(*,ion)
  ;;   endfor
  ;;   ;
  ;;   print,' hit any key to continue'
  ;;   x=get_kbrd(1)
  ;; endif
 ;
 ;
  ioneq=dblarr(ntemp,z+2)
;
  for it=0,ntemp-1 do begin
    factor=fltarr(z+1)
   ;
    for ion=1,z-1 do begin
      rat=ion_rate(it,ion)/rec_rate(it,ion)
      factor(ion)=rat^2+rat^(-2)
    endfor
   ;
    factor(0)=max(factor)
    factor(z)=max(factor)
   ;
    idx=where(factor eq min(factor))
   ;
    most=idx(0)
    ioneq(it,most)=1.d
   ;
   ; Get ions above most
   ;
    for ion=most+1,z+1 do begin
      if rec_rate(it,ion) gt 0. then begin
        ioneq(it,ion)=ion_rate(it,ion-1)*ioneq(it,ion-1)/rec_rate(it,ion)
      endif else ioneq(it,ion)=0.
      if ioneq(it,ion) lt ioneqmin then ioneq(it,ion)=0.    
    endfor
   ;
    for ion=most-1,0,-1 do begin
      ioneq(it,ion)=rec_rate(it,ion+1)*ioneq(it,ion+1)/ion_rate(it,ion)
      if ioneq(it,ion) lt ioneqmin then ioneq(it,ion)=0.
    endfor
   ;
    ioneq(it,0)=ioneq(it,*)/total(ioneq(it,*))
   ;
  endfor   ; it
 ;
  ;; if ntemp gt 1 then begin
  ;;   xr=minmax(temp)
  ;;   yr=[1.e-5,1.]
  ;;   plot_oo,fltarr(2),fltarr(2),xr=xr,yr=yr,/nodata
  ;;   for ion=0,z do begin
  ;;     oplot,temp,ioneq(*,ion)
  ;;   endfor
  ;; endif
 ;
  for ion=0,z do begin
    printf,luw,z,ion+1,ioneq(*,ion),format='(2i3,'+iformat+')'
  endfor
;
endfor  ; loop over z
;
printf,luw,' -1'
printf,luw,'%filename:  ',file_basename(outname)
printf,luw,'%ionization equilibrium:  CHIANTI'
printf,luw,' produced as part of the CHIANTI atomic data base collaboration'
datetime=systime()
printf,luw,'  Created on ',datetime
printf,luw,' -1'
free_lun,luw

print,'% MAKE_IONEQ_ALL: the ion fractions have been written to the file'
print,'                  '+outname


k=where(temp LT 1e4,nk)
IF nk GT 0 THEN BEGIN
  print,'**WARNING**'
  print,'Ion fractions have been computed below 10^4 K. These values may not be accurate as charge'
  print,'exchange can be significant at low temperatures and this process is not included in CHIANTI.'
ENDIF 

END
