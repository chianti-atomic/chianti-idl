

PRO which_line, ionname, wvl, narrow=narrow, all=all

;+
; NAME
;
;     WHICH_LINE
;
; PROJECT
;
;     CHIANTI
;
; PURPOSE:
;
;     Upon given an ion name and wavelength, this routine prints out a list
;     of possible line IDs for the wavelength. Wavelengths within 1% of the
;     input wavelength are searched for.
;
; INPUTS
;
;     IONNAME  Name of an ion in the CHIANTI format. E.g., 'fe_13' for Fe XIII.
;
;     WVL      A wavelength in angstroms.
;
; OUTPUTS (to screen)
;
;     Prints a list of atomic transitions and wavelengths for lines close to
;     the input wavelength. A '*' is placed next to the closest wavelength
;     match.
;
; KEYWORDS
;
;     NARROW   Narrows the search range to 0.02% of the specified
;     wavelength.
;
;     ALL      If set, then lines with theoretical wavelengths are
;              included in the check.
;
; EXAMPLE
;
;     IDL> which_line,'o_6',1032
;     Wavelength   i   j Lower level           Upper level             A-value
;       1037.615   1   2 1s2.2s 2S1/2        - 1s2.2p 2P1/2          4.21e+008
;       1031.914   1   3 1s2.2s 2S1/2        - 1s2.2p 2P3/2          4.28e+008
;
; CALLS
;
;     CONVERTNAME, ZION2FILENAME, READ_WGFA2, READ_ELVLC
;
; HISTORY
;
;     Ver.1, 22-Jun-2004, Peter Young
;     Ver.2, 16-Jul-2013, Peter Young
;        now works for the dielectronic ions.
;     Ver.3, 20-Feb-2014, Peter Young
;        modified to only print lines with observed wavelengths (use
;        /all to show all wavelengths).
;-

IF n_params() LT 2 THEN BEGIN
  print,'Use: IDL> which_line, ionname, wvl [, /narrow, /all] '
  return
ENDIF

convertname,ionname,iz,ion,diel=diel
zion2filename,iz,ion,filename,diel=diel

wgfaname=filename+'.wgfa'
elvlcname=filename+'.elvlc'

read_wgfa2,wgfaname,lvl1,lvl2,wavel,gf,a_value,ref
read_elvlc,elvlcname,l1,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,ref

IF keyword_set(narrow) THEN factor=0.002 ELSE factor=0.01

chck=wvl*factor

;
; Find lines within wavelength range
;
ind=where(abs(abs(wavel)-wvl) LE chck,n_ind)

;
; Filter out those lines with negative wavelengths.
;
IF NOT keyword_set(all) AND n_ind NE 0 THEN BEGIN
  k=where(wavel[ind] GT 0.,nk)
  IF nk NE 0 THEN ind=ind[k] ELSE n_ind=0
ENDIF 

IF n_ind EQ 0 THEN BEGIN
  print,'% WHICH_LINE: no lines found within 1% of the wavelength '+string(wvl)
  return
ENDIF

k=sort(abs(wavel[ind]))
ind=ind[k]
getmin=min(abs(abs(wavel[ind])-wvl),imin)

n=n_elements(ind)

print,'     Wavelength   i   j  Lower level           Upper level             A-value'
FOR i=0,n-1 DO BEGIN
  j=ind[i]
  term1=strpad(term[lvl1[j]-1],20,/after,fill=' ')
  term2=strpad(term[lvl2[j]-1],20,/after,fill=' ')
  IF i EQ imin THEN str1='*' ELSE str1=''
  print,format='(a1,f14.3,2i4,"  ",a20,"- ",a20,e11.2)', $
       str1,wavel[j],lvl1[j],lvl2[j], $
       term1,term2,a_value[j]
ENDFOR

IF NOT keyword_set(all) THEN BEGIN
  print,''
  print,'Use keyword /all to include lines with theoretical wavelengths.'
ENDIF 

END
