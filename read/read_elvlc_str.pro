;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), 
;       Cambridge University (United Kingdom), George Mason University (USA), and
;       the University of Michigan (USA).
;
;
; NAME:
;	READ_ELVLC_STR
;
; PURPOSE:
;
;	to read files containing observed and theoretical energy levels and return an IDL
;   structure with the energy level information
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_ELVLC, File, ElvlcStr, Ref
;
;
; INPUTS:
;
;	File:	the name of the file 
;               i.e., !xuvtop+'/si/si_12/si_12.elvlc' for Si XII
;
; OPTIONAL INPUTS:
;
;	None:
;	
;
; OUTPUTS: 
;       ElvlcStr, an IDL structure with the following tags:
;                filename - the filename
;                l1      - level index
;                term    - configuration index
;                conf    - configuration description
;                ss      - 2S+1
;                ll      - L
;                jj      - J
;                ecm     - observed energy (cm^-1)
;                eryd    - observed energy (Rydbergs)
;                ecmth   - theoretical energy (cm^-1)
;                erydth  - theoretical energy (Rydbergs)
;   and
;       Ref     1D string array of references to the data in the scientific literature
;               
;   note:  the theoretical energies are usually those used in the scattering
;          calculation and are only useful for predicting approximate 
;          wavelengths
;
;
; EXAMPLE:
;
;             > file = !xuvtop+'/si/si_12/si_12.elvlc'
;             > read_elvlc,file,elvlcstr,ref
;             > 
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	November 2010:     Version 1.0
;
;
;   this is a modified version of read_elvlc to return all variables in a structure
;-
pro read_elvlc_str,filename,elvlcstr,ref
;
;
;
if n_params(0) lt 3 then begin
   print,''
   print,' type> read_elvlc,filename,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref'
   print,''
   return
endif
;
;
;
openr,lue,filename,/get_lun
;
l1a=intarr(1000)
confa=intarr(1000)
desiga=strarr(1000)
ssa=intarr(1000)
lla=intarr(1000)
spda=strarr(1000)
jja=fltarr(1000)
multa=intarr(1000)
ecma=dblarr(1000)
eryda=dblarr(1000)
ecmtha=dblarr(1000)
erydth=dblarr(1000)
;
;
string1=' '
l1=1  & desig=' ' & ll=1. & conf=1 & ss=1. & jj=1. & mult=1 & ecm=1.d
spd=' ' & eryd=1.d & ecmzh=1.d  & erydzh=1.d
;
;
while (strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 10) do begin
readf,lue,string1
if (strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 10) then begin

  reads,string1,l1,conf,desig,ss,ll,spd,jj,mult,ecm,eryd,ecmzh,erydzh,  $
   format='$(i3,i6,a15,2i3,a3,f4.1,i3,f15.3,f15.6,f15.3,f15.6,f15.3,f15.6)'
l=l1-1
  l1a(l)=l1
  confa(l)=conf
  desiga(l)=desig
  ssa(l)=ss
  lla(l)=ll
  spda(l)=spd
  jja(l)=jj
  multa(l)=float(mult)
  ecma(l)=ecm
  eryda(l)=eryd
  ecmtha(l)=ecmzh
  erydth(l)=erydzh
endif
endwhile
;
nlvls=max(l1a)
l1a=l1a(0:nlvls-1)
confa=confa(0:nlvls-1)
ssa=ssa(0:nlvls-1)
lla=lla(0:nlvls-1)
spda=spda(0:nlvls-1)
jja=jja(0:nlvls-1)
ecma=ecma(0:nlvls-1)
eryda=eryda(0:nlvls-1)
ecmtha=ecmtha(0:nlvls-1)
erydth=erydth(0:nlvls-1)
;
term=strarr(nlvls)


;
; The section below has been modified to give the J-value in fractional 
; form. PRY.
;
FOR lvl1=0,nlvls-1 DO BEGIN
  IF (fix(2*jja[lvl1]/2) NE 2*jja[lvl1]/2) THEN BEGIN
    jstring=strtrim(string(format='(i4)',2*jja[lvl1]),2)+'/2'
  ENDIF ELSE BEGIN
    jstring=strtrim(string(format='(i4)',2*jja[lvl1]/2),2)
  ENDELSE
;  term(lvl1)=strtrim(desiga(lvl1),2)+' '+string(ssa(lvl1),'(i1)') $
;     +strtrim(spda(lvl1),2) +string(jja(lvl1),'(f3.1)')
  term(lvl1)=strtrim(desiga(lvl1),2)+' '+string(ssa(lvl1),'(i1)') $
       +strtrim(spda(lvl1),2) +jstring
ENDFOR
;
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
elvlcstr = {  $
    filename:filename, $
    l1:l1a, $
    term:term, $
    conf:confa, $
    ss:ssa,$
    ll:lla, $
    spd:spda, $
    jj:jja, $
    ecm:ecma,  $
    eryd:eryda, $
    ecmth:ecmtha, $
    erydth:erydth }
;
free_lun,lue
;
;

end
