;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It
;       is a collaborative project involving George Mason University USA),
;       the University of Michigan (USA), and Cambridge University (UK).
;
;
; NAME:
;	READ_IP
;
; PURPOSE:
;
;	to read values of ionization potentials 
;
;
; CALLING SEQUENCE:
;
;       READ_IP, File, IP, Ref
;
;
; INPUTS:
;
;	File:  the name of the file containing the IP values, usually
;               !xuvtop/ip/chianti.ip	
;
;
; OUTPUTS:
;
;	IP:  Array values of ionization potential (cm^-1)
;       Ref:  the reference to the IP values in the scientific literature
;
;
;
;
;
; EXAMPLE:
;
;             > read_ip,!xuvtop+'/ip/chianti.ip',ip,ref
;                  ip(2-1,2-1) give the ionization potential of He II  (Z=2, Ion=2)
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1998:     Version 1.0
;
;-
pro read_ip,ipname,ip,ref
;
;
if n_params() lt 2 then begin
   print,' '
   print,'  IDL> read_ip,!xuvtop+''/ip/chianti.ip'',ip,ref
   print,'         ip(2-1,2-1) give the ionization potential (cm^-1) of He II  (Z=2, Ion=2)'
   print,' '
   return
endif
   
;
ip=dblarr(30,30)
;
openr,lu,filename,/get_lun
;
;
iz=1
ion=1
;
string1=' '
ip1=1.d
;
index=0
;
while strpos(string1, '-1') lt 0 do begin
readf,lu,string1
if(strpos(string1,'-1') lt 0) then begin

  reads,string1,iz,ion,ip1,format='(2i5,f20.5)
  ip(iz-1,ion-1)=ip1
endif
endwhile
;
;
;  get references
refstring=strarr(100)
nref=0
;
string1=' '
while strpos(string1,'-1') lt 0  do begin
readf,lu,string1
if(strpos(string1,'-1') lt 0) and (strpos(string1,'%file') lt 0) then begin
  refstring(nref)=string1
  nref=nref+1
endif
endwhile
;
ref=refstring(0:nref-1)
;
;
;
;
free_lun,lu
;
;
end
