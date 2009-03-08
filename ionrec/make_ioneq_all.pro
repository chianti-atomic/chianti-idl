
;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;               Florence (Italy), and the University of Cambridge .
;
;
; NAME
;
;     MAKE_IONEQ_ALL
;
;  PURPOSE:
;
;     Calculates the ionization equilibrium for all ions with Z=1->30
;     Stores the result in the users home direction.  the file name is specified
;     by the variable 'outname', with a default value of 'new.ioneq'
;
; INPUTS
;
;
;    TEMPERATURE    Temperature (K).
;
;
; OUTPUTS
;
;    None
;
; OPTIONAL INPUTS
;
;    NONE
;
; KEYWORDS
;
;    outname
;
; CALLS
;
;
; COMMON BLOCKS
;
;    NONE
;
; PROGRAMMING NOTES
;
;    NONE
;
; MODIFICATION HISTORY
;
;    Ver.1, 3-March-2009, Ken Dere
;
; VERSION     :  1, 3-March-2009
;
;-
pro make_ioneq_all,temp,outname='new.ioneq'
;
if n_params() lt 1 then begin
   print,' > make_ioneq_all, temp'
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
;
home = getenv('HOME')
name=concat_dir(home,outname)
openw,luw,name,/get_lun
printf,luw,ntemp,30,format='(2i3)'
printf,luw,alog10(temp),format='(1000f6.2)'
;
for z=1,30 do begin
;       
el=strlowcase(zlabl(z-1))
;
ion_rate=fltarr(ntemp,z+2)
;
rec_rate=fltarr(ntemp,z+2)
;
;
;  get ionization rates
;
for ion=0,z-1 do begin
  zion2name,z,ion+1,gname
;   print,' ionizdat gname = ',gname,iz,ion
   ion_rate(0,ion) = ionizrate(gname,temp)
endfor
;
;
yr=[1.e-10,1.]*max(ion_rate)
plot_oo,fltarr(2),fltarr(2),xr=[1.e+3,1.e+9],yr=yr,/nodata,title=zlabl(z-1)
;
for ion=0,z-1 do begin
      oplot,temp,ion_rate(*,ion)
endfor
;
;  get recombination rates
;
for ion=1,z do begin
	zion2name,z,ion+1,gname
;   print,' recombdat gname = ',gname,iz,ion
	rec_rate(0,ion)=recombrate(gname,temp)
endfor
;
;
for ion=1,z do begin
       oplot,temp,rec_rate(*,ion)
endfor
;
print,' hit any key to continue'
x=get_kbrd(1)
;
;
;
;
;
ioneq=dblarr(ntemp,z+2)
;
for it=0,ntemp-1 do begin
;  print,' t = ',t(it)
;
	factor=fltarr(z+1)
;
	for ion=1,z-1 do begin
 	  rat=ion_rate(it,ion)/rec_rate(it,ion)
 	  factor(ion)=rat^2+rat^(-2)
	;   print,' ion, rat = ',ion,rat,factor(ion),
	endfor
;
	factor(0)=max(factor)
	factor(z)=max(factor)
;
	idx=where(factor eq min(factor))
	; print,' min factor = ', idx(0),factor(idx(0))
;
	most=idx(0)
	ioneq(it,most)=1.d
;
	;  get ions above most
;
	for ion=most+1,z+1 do begin
		if rec_rate(it,ion) gt 0. then begin
			ioneq(it,ion)=ion_rate(it,ion-1)*ioneq(it,ion-1)/rec_rate(it,ion)
		endif else ioneq(it,ion)=0.
   
	;   if ioneq(it,ion) lt ioneqmin then ioneq(it,ion)=0.
	;   print,it,ion,ioneq(it,ion),ion_r(it,ion-1),rec_r(it,ion)
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
xr=minmax(temp)
yr=[1.e-5,1.]
plot_oo,fltarr(2),fltarr(2),xr=xr,yr=yr,/nodata
for ion=0,z do begin
   oplot,temp,ioneq(*,ion)
endfor
;
;
for ion=0,z do begin
   printf,luw,z,ion+1,ioneq(*,ion),format='(2i3,1000e10.3)'
endfor
;
endfor  ; loop over z
;
printf,luw,' -1'
printf,luw,'%filename:  ',name
printf,luw,'%ionization equilibrium:  CHIANTI'
printf,luw,' produced as part of the CHIANTI atomic data base collaboration'
datetime=systime()
printf,luw,'  K.P. Dere (GMU)  ',datetime
printf,luw,' -1'
free_lun,luw
end
