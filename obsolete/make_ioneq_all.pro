;+
; NAME:
;	MAKE_IONEQ_ALL
;
; PURPOSE:
;       Compute an equilibrium ion fraction table in a format suitable
;       for input to the CHIANTI database.
;
; CATEGORY:
;       CHIANTI; ionization; recombination, ion fractions.
;
; CALLING SEQUENCE:
;	MAKE_IONEQ_ALL, Temp
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       Temp:  1D array specifying the temperatures (in K) for which
;              the ion fraction curves are needed. If not
;              specified, then the default CHIANTI range of logT=4.0
;              to logT=9.0 in 0.05 dex intervals is 
;              used. 
;       Outname:  The name of the output file. If not specified, then
;                 new ion fraction file is sent to 'new.ioneq'.
;       Density:  Specifies the electron number density (units: cm^-3)
;                 for the which the ion fractions are calculated.
;       Pressure: Specifies the pressure (N_e*T; units cm^-3 K) at
;                 which the ion fractions are calculated.
;
; OUTPUTS:
;       Creates the file 'new.ioneq' in the current working directory
;       (use OUTNAME= to specify a different filename). The file
;       contains ion fractions of all ions of all elements up to and
;       including zince.
;
; EXAMPLE:
;       
;	Please provide a simple example here. An example from the
;	DIALOG_PICKFILE documentation is shown below. Please try to
;	include examples that do not rely on variables or data files
;	that are not defined in the example code. Your example should
;	execute properly if typed in at the IDL command line with no
;	other preparation. 
;
; MODIFICATION HISTORY:
;    Ver.1, 3-March-2009, Ken Dere
;    Ver.2, 19-May-2010, Ken Dere
;         correct various errors
;    Ver.3, 19-Jul-2016, Peter Young
;         modified implementation of OUTNAME; introduced ion_rate= and
;         rec_rate= outputs; plots are no longer made; added a warning
;         if ion fractions below 10^4 K are computed.
;    Ver.4, 27-Jul-2017, Peter Young
;         added PRESSURE and DENSITY optional inputs; removed ION_RATE
;         and REC_RATE as they were'nt implemented properly;
;         updated header.
;    Ver.5, 05-Apr-2023, Peter Young
;         Previously the routine required TEMP to be specified. This is
;         no longer needed.
;-


PRO make_ioneq_all,temp,outname=outname, $
                   density=density, pressure=pressure

;
ioneqmin=1.e-20
;
;   get ionization file list
;
zlabl=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
       'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr', $
       'Mn','Fe','Co','Ni','Cu','Zn']
;


;
; If TEMP not set, then use the default CHIANTI values.
;
IF n_elements(temp) EQ 0 THEN BEGIN
  log_temp=findgen(101)/20.+4.0
  temp=10.^log_temp
ENDIF 
ntemp=n_elements(temp)
tformat='('+trim(ntemp)+'f6.2)'
iformat='('+trim(ntemp)+'e10.3)'


IF n_elements(outname) EQ 0 THEN outname='new.ioneq'
;
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
    ion_rate(*,ion) = ioniz_rate(gname,temp)
  endfor
 ;
 ;
 ; Get recombination rates
 ;
  for ion=1,z do begin
    zion2name,z,ion+1,gname
    rec_rate(*,ion)=recomb_rate(gname,temp,density=density, pressure=pressure)
  endfor
 ;
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
