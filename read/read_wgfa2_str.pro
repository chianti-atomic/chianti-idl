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
;	READ_WGFA2_STR
;
; PURPOSE:
;
;	read radiative data files
;         a modified version of read_wgfa
;         needed to take account of two types of transitions between the same 2 levels:
;         for example, the M1 and the 2 photon E1 transition 1s-2s in hydrogenic ions.
;
;
; CATEGORY:
;
;	science
;
; CALLING SEQUENCE:
;
;       READ_WGFA2, File, Wgfastr  Ref
;
;
; INPUTS:
;
;	File:  name of the file containing the radiative data
;                i.e. !xuvtop/c/c_4/c_4.wgfa
;
;
; OUTPUTS:
;   Wgfastr, a structure with the following tags
;
;       lvl1:  1D array of indices of the lower level (starting at 1)
;       lvl2:  1D array of indices of the upper level (starting at 1)
;       wvl:   1D array of transition wavelengths in Angstroms
;       gf:    1D array of weighted oscillator strength gf
;       a_value:  1D array of the total radiative transition probability (s^-1)
; and
;   Ref:   1D string array of references to the data in the scientific literature
;
;
;
; EXAMPLE:
;
;             > read_wgfa2_str, wgfastr, ref
;             
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere (GMU)
;	17 November 2010:     Version 1.0
;
;
;-
pro read_wgfa2_str,filename, wgfastr, ref
;
;
;   to read data in *.wgfa (gf, A value) files
;
;  a modified version of read_wgfa
;     needed to take account of two types of transitions between the same 2 levels:
;         for example, the M1 and the 2 photon E1 transition 1s-2s in hydrogenic ions
;      and return a structure
;
;
if n_params() lt 3 then begin
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
if file_test(filename) eq 1 then begin
  status = 1
  
  openr,lue,ename,/get_lun
  ;
  lvl1=intarr(100000)
  lvl2=intarr(100000)
  wvl1=fltarr(100000)
  gf1=fltarr(100000)
  a_value1=fltarr(100000)
  ;
  ;
  ;
  index=0L
  string1=' '
  while strpos(string1,'-1') lt 0 or strpos(string1,'-1') gt 10 do begin
    readf,lue,string1
    ;  print,string1
    if(strpos(string1,'-1') lt 0 or strpos(string1,'-1') gt 10) then begin
    
      reads,string1,l1,l2,wvl,gf,a_value,  $
        format='$(2i5,f15.3,2e15.3)'
      lvl1(index)=l1
      lvl2(index)=l2
      wvl1(index)=wvl
      gf1(index)=gf
      a_value1(index)=a_value
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
  wvl1=wvl1(0:nindex-1)         ;  wavelength in Angstroms
  gf1=gf1(0:nindex-1)          ;  gf
  a_value=a_value1(0:nindex-1)     ;  Einstein A coefficient
  ;
  nlvls=max([lvl1,lvl2])
  wvl=fltarr(nlvls,nlvls)
  gf=fltarr(nlvls,nlvls)
  a_value=fltarr(nlvls,nlvls)
  
  ind1 = where(wvl1 EQ 0.)
  ind2 = where(wvl1 NE 0.)
  
  wvl[lvl1-1,lvl2-1]= abs(wvl1)
  gf[lvl1-1,lvl2-1]=gf1
  
  IF ind1[0] NE -1 THEN BEGIN
    a_value[lvl1[ind1]-1,lvl2[ind1]-1] = a_value1[ind1]
    a_value[lvl1[ind2]-1,lvl2[ind2]-1] = $
      a_value[lvl1[ind2]-1,lvl2[ind2]-1] + a_value1[ind2]
  ENDIF ELSE BEGIN
    a_value[lvl1[ind2]-1,lvl2[ind2]-1] = a_value1[ind2]
  ENDELSE
  wgfastr = {lvl1:lvl1, $
    lvl2:lvl2, $
    wvl:wvl, $
    gf:gf, $
    a_value:a_value, $
    status:status }
endif else begin
  status = 0
  wgfastr = {status:status}
endelse 
    
;
end
