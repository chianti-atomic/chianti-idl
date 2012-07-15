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
;	READ_WGFA2
;
; PURPOSE:
;
;	read radiative data files
;         a modified version of read_wgfa
;         needed to take account of two types of transitions between the same 2 levels:
;         for example, the M1 and the 2 photon E1 transition 1s-2s in hydrogenic ions
;
;
; CATEGORY:
;
;	science
;
; CALLING SEQUENCE:
;
;       READ_WGFA2, File, Lvl1, Lvl2, Wvl, Gf, A_value, Ref
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
;       Wvl:   1D array of transition wavelengths in Angstroms
;       Gf:    1D array of weighted oscillator strength gf
;       A_value:  1D array of the total radiative transition probability (s^-1)
;       Ref:   1D string array of references to the data in the scientific literature
;
;
;
; EXAMPLE:
;
;             > read_wgfa2,!xuvtop+'/c/c_4/c_4.wgfa',lvl1,lvl2,wvl,gf,a,ref
;             
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       V.3, 3-Nov-03, Giulio Del Zanna (GDZ)
;        Modified size of input arrays. Now the routine can read a file up to
;        30 000  lines.
;
;       V.4, 26-May-09, Enrico Landi (EL)
;        Modified size of input arrays. Now the routine can read a file up to
;        100,000 lines.
;        Modified definition of variable index (now index=0L).
;
;       V.5, 26-May-09, Peter Young
;        Modified the way references are read.
;
;          
; VERSION     :  5, 26-May-09
;
;-
pro read_wgfa2,filename,lvl1,lvl2,wvl,gf,a_value,ref
;
;
;   to read data in *.wgfa (gf, A value) files
;
;  a modified version of read_wgfa
;     needed to take account of two types of transitions between the same 2 levels:
;         for example, the M1 and the 2 photon E1 transition 1s-2s in hydrogenic ions
;
;
if n_params() lt 7 then begin
   print,' '
   print,'   IDL> read_wgfa2,filename,lvl1,lvl2,wvl,gf,a_value,ref'
   print,'  i.e.> read_wgfa2,!xuvtop+''/c/c_4/c_4.wgfa'',lvl1,lvl2,wvl,gf,a,ref'
   print,'         to read the gf and A values for C IV'
   print,'         the arrays are structured differently than with read_wgfa'
   print,' '
   return
endif
;
;
ename=filename
;
openr,lue,ename,/get_lun
;
lvl1=intarr(100000)             
lvl2=intarr(100000)            
wvl=fltarr(100000)            
gf=fltarr(100000)            
a_value=fltarr(100000)      
;
;
;
index=0L                   
string1=' '
while strpos(string1,'-1') lt 0 or strpos(string1,'-1') gt 10 do begin
readf,lue,string1
;  print,string1
if(strpos(string1,'-1') lt 0 or strpos(string1,'-1') gt 10) then begin

  reads,string1,l1,l2,wvl1,gf1,a_value1,  $
   format='$(2i5,f15.3,2e15.3)'
    lvl1(index)=l1
    lvl2(index)=l2
    wvl(index)=wvl1
    gf(index)=gf1
    a_value(index)=a_value1
    index=index+1
endif
endwhile
;
;
; Get references
;
tst1=0
ref=''
str1=''
WHILE tst1 EQ 0 DO BEGIN
  readf,lue,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    tst1=1
  ENDIF ELSE BEGIN
    ref=[ref,str1]
  ENDELSE
ENDWHILE
ref=ref[1:*]

free_lun,lue
;
nindex=index   ;nindex=index-1
nlvls=max([lvl1,lvl2])
;
lvl1=lvl1(0:nindex-1)           ;  lower level #
lvl2=lvl2(0:nindex-1)           ;  upper level #
wvl=wvl(0:nindex-1)         ;  wavelength in Angstroms
gf=gf(0:nindex-1)          ;  gf
a_value=a_value(0:nindex-1)     ;  Einstein A coefficient
;
;
end
