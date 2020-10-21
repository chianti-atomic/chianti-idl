

PRO which_line, ionname, wvl, narrow=narrow, all=all, path=path

;+
; NAME:
;     WHICH_LINE
;
; PURPOSE:
;     Upon given an ion name and wavelength, this routine prints out a list
;     of possible line IDs for the wavelength. Wavelengths within 1% of the
;     input wavelength are searched for.
;
; CATEGORY:
;     CHIANTI; search.
;
; CALLING SEQUENCE:
;     WHICH_LINE, IonName, Wvl
;
; INPUTS:
;     IonName: Name of an ion in the CHIANTI format. E.g., 'fe_13' for
;              Fe XIII. 
;     Wvl:     A wavelength in angstroms.
;
; OPTIONAL INPUTS:
;	Parm2:	Describe optional inputs here. If you don't have any, just
;		delete this section.
;	
; KEYWORD PARAMETERS:
;     NARROW:  Narrows the search range to 0.02% of the specified
;              wavelength.
;     ALL:     If set, then lines with theoretical wavelengths are
;              included in the check.
;
; OUTPUTS:
;     Prints a list of atomic transitions and wavelengths for lines close to
;     the input wavelength. A '*' is placed next to the closest wavelength
;     match.
;
; CALLS:
;     CONVERTNAME, ZION2FILENAME, READ_WGFA2, READ_ELVLC
;
; EXAMPLE:
;     IDL> which_line,'o_6',1032
;     Wavelength   i   j Lower level           Upper level             A-value
;       1037.615   1   2 1s2.2s 2S1/2        - 1s2.2p 2P1/2          4.21e+008
;       1031.914   1   3 1s2.2s 2S1/2        - 1s2.2p 2P3/2          4.28e+008
;
; MODIFICATION HISTORY:
;     Ver.1, 22-Jun-2004, Peter Young
;     Ver.2, 16-Jul-2013, Peter Young
;        now works for the dielectronic ions.
;     Ver.3, 20-Feb-2014, Peter Young
;        modified to only print lines with observed wavelengths (use
;        /all to show all wavelengths).
;     Ver.4, 25-Jun-2020, Peter Young
;        added PATH= optional input; updated header format.
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use: IDL> which_line, ionname, wvl [, /narrow, /all, path=] '
  return
ENDIF

convertname,ionname,iz,ion,diel=diel

IF n_elements(path) NE 0 THEN BEGIN
  chck=file_info(path)
  IF chck.exists EQ 0 OR chck.directory EQ 0 THEN BEGIN
    print,'% WHICH_LINE: the specified PATH does not exist. Returning...'
    return
  ENDIF
  basename=concat_dir(path,ionname)
ENDIF ELSE BEGIN
  zion2filename,iz,ion,basename,diel=diel
ENDELSE 
;
wgfaname=basename+'.wgfa'
elvlcname=basename+'.elvlc'

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
  IF keyword_set(narrow) THEN perc_string='0.2%' ELSE perc_string='1%'
  print,'% WHICH_LINE: no lines found within '+perc_string+' of the wavelength '+trim(string(wvl))+' angstroms.'
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
