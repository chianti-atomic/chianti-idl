;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
;
; NAME:
;	READ_ELVLC_DIRECT
;
; PURPOSE:
;
;	to read files containing observed energy levels
;       does not reformat some variables like READ_ELVL
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_ELVLC_DIRECT, File, L1, Term, Conf, ss, ll, spd, jj, Ecm, Eryd, Ecmth, Erydth, Ref
;
;
; INPUTS:
;
;	File:	the name of the file 
;               i.e., !xuvtop+'/si/si_12/si_12.elvl' for Si XII
;
; OPTIONAL INPUTS:
;
;	None:
;	
;
; OUTPUTS:       L1    - level index
;                Term  - configuration index
;                Conf  - configuration description
;                ss    - 2S+1
;                ll    - L
;                spd   - 'S', 'P', etc to denote L value
;                jj    - J
;                Mult  - multiplicity  2J+1
;                Ecm   - energy (cm^-1)
;                Eryd  - energy (Rydbergs)
;                Ecmth - energy (cm^-1)
;                Erydth- energy (Rydbergs)
;                Ref   - reference
;               
;
;
;
;
; EXAMPLE:
;
;             > file = !xuvtop+'/si/si_12/si_12.elvl'
;             > read_elvlc_direct,file,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref
;             > 
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;-
pro read_elvlc_direct,filename,l1,conf,term,ss,ll,spd,jj,mult,ecm,eryd,ecmth,erydth,ref
;
;
if n_params(0) lt 14 then begin
     print,' '
     print,' IDL> read_elvlc_direct,filename,l1,conf,term,ss,ll,spd,jj,mult,ecm,eryd,ecmth,erydth,ref'
     print,' '
     return
endif
;
;
openr,lue,filename,/get_lun
;
l1=intarr(1000)
conf=intarr(1000)
desig=strarr(1000)
ss=intarr(1000)
ll=intarr(1000)
spd=strarr(1000)
jj=fltarr(1000)
mult=fltarr(1000)
ecm=dblarr(1000)
eryd=dblarr(1000)
ecmth=dblarr(1000)
erydth=dblarr(1000)
;
;
string1=' '
l11=1  & desig1=' ' & ll1=1. & conf1=1 & ss1=1. & jj1=1. & mult1=1 & ecm1=1.d
spd1=' ' & eryd1=1.d & ecmth1=1.d  & erydth1=1.d
;
;
while (strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 10) do begin
readf,lue,string1
if(strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 10) then begin

  reads,string1,l11,conf1,desig1,ss1,ll1,spd1,jj1,mult1,ecm1,eryd1,ecmth1,erydth1,  $
   format='$(i3,i6,a15,2i3,a3,f4.1,i3,f15.3,f15.6,f15.3,f15.6,f15.3,f15.6,f15.3,f15.6)'
  l=l11-1
  l1(l)=l11
  conf(l)=conf1
  desig(l)=desig1
  ss(l)=ss1
  ll(l)=ll1
  spd(l)=spd1
  jj(l)=jj1
  mult(l)=float(mult1)
  ecm(l)=ecm1
  eryd(l)=eryd1
  ecmth(l)=ecmth1
  erydth(l)=erydth1
endif
endwhile
;
;  get references
refstring=strarr(100)
nref=0
;
string1=' '
while (strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 10)  do begin
readf,lue,string1
if (strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 10) and (strpos(string1,'%file') lt 0) then begin
  refstring(nref)=string1
  nref=nref+1
endif
endwhile
;
ref=refstring(0:nref-1)
;
;
nlvls=max(l1)
l1=l1(0:nlvls-1)
conf=conf(0:nlvls-1)
term=desig(0:nlvls-1)
ss=ss(0:nlvls-1)
ll=ll(0:nlvls-1)
spd=spd(0:nlvls-1)
jj=jj(0:nlvls-1)
mult=mult(0:nlvls-1)
ecm=ecm(0:nlvls-1)
eryd=eryd(0:nlvls-1)
ecmth=ecmth(0:nlvls-1)
erydth=erydth(0:nlvls-1)
;
;term=strarr(nlvls)
;;
;for lvl1=0,nlvls-1 do begin
;  term(lvl1)=strtrim(desig(lvl1),2)
;endfor
;
;
free_lun,lue
;
;
end
