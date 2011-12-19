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
;	READ_WGFA
;
; PURPOSE:
;
;	:
;
; CATEGORY:
;
;	science
;
; CALLING SEQUENCE:
;
;       READ_WGFA, File, Lvl1, Lvl2, Wvl, Gf, A_value, Ref
;
;
; INPUTS:
;
;	File:  name of the file containing the radiative data
;                i.e. !xuvtop/c/c_4/c_4.wgfa
;
;
; OUTPUTS:
;
;	Lvl1:  1D array of indices of the lower level (starting at 1)
;       Lvl2:  1D array of indices of the upper level (starting at 1)
;       Wvl:   2D array of transition wavelengths in Angstroms
;       Gf:    2D array of weighted oscillator strength gf
;       A_value:  2D array of radiative transition probability (s^-1)
;       Ref:   1D string array of references to the data in the scientific literature
;
;
;
; EXAMPLE:
;
;             > read_wgfa,!xuvtop+'/c/c_4/c_4.wgfa',lvl1,lvl2,wvl,gf,a,ref
;             
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       v.3, 6-Mar-02, Peter Young
;           Changed the way the reference string is read in order to 
;           prevent '-1' problems.
;
;       v.4, 12-Mar-02, Peter Young
;           Corrected bug following above change.
;
;-
pro read_wgfa,filename,lvl1,lvl2,wvl,gf,a_value,ref
;
;
;   to read data in *.wgfa (gf, A value) files
;
if n_params() lt 7 then begin
   print,' '
   print,'   IDL> read_wgfa,filename,lvl1,lvl2,wvl,gf,a_value,ref'
   print,'  i.e.> read_wgfa,!xuvtop+''/c/c_4/c_4.wgfa'',lvl1,lvl2,wvl,gf,a,ref'
   print,'         to read the gf and A values for C IV'
   print,' '
   return
endif
;
ename=filename
;
openr,lue,ename,/get_lun
;
lvl1=intarr(10000)
lvl2=intarr(10000)
wvl0=fltarr(10000)
gf0=fltarr(10000)
a_value0=fltarr(10000)
;
;
;
index=0
string1=' '
while strpos(string1,'-1') lt 0 or strpos(string1,'-1') gt 10 do begin
readf,lue,string1
;  print,string1
if(strpos(string1,'-1') lt 0 or strpos(string1,'-1') gt 10) then begin

  reads,string1,l1,l2,wvl1,gf1,a_value1,  $
   format='$(2i5,f15.3,2e15.3)'
     lvl1(index)=l1
     lvl2(index)=l2
     wvl0(index)=wvl1
     gf0(index)=gf1
     a_value0(index)=a_value1
     index=index+1
endif
endwhile


refstring=''
tst1=1
WHILE tst1 EQ 1 DO BEGIN
  readf,lue,string1
  IF strcompress(string1,/rem) NE '-1' THEN refstring=[refstring,string1] $
       ELSE tst1=2
ENDWHILE

free_lun,lue

ref=refstring[1:*]

nindex=index   ;nindex=index-1
nlvls=max([lvl1,lvl2])
;
lvl1=lvl1(0:nindex-1)           ;  lower level #
lvl2=lvl2(0:nindex-1)           ;  upper level #
wvl=fltarr(nlvls,nlvls)         ;  wavelength in Angstroms
gf=fltarr(nlvls,nlvls)          ;  gf
a_value=fltarr(nlvls,nlvls)     ;  Einstein A coefficient
;
for index=0,nindex-1 do begin
    l1=lvl1(index)-1
    l2=lvl2(index)-1
    wvl(l1,l2)=wvl0(index)
    gf(l1,l2)=gf0(index)
    a_value(l1,l2)=a_value0(index)
endfor
;
end
