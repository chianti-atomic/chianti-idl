
PRO read_wgfa_str, wgfaname, str, ref

;+
; NAME
;
;    READ_WGFA_STR
;
; PROJECT
;
;    CHIANTI
;
; EXPLANATION
;
;    Reads the CHIANTI .wgfa file into a structure. For the most part
;    CHIANTI wgfa files contain radiative decay rates. I.e., rates for
;    transitions that yield photons of a definite wavelength. However
;    for some ions the .wgfa file is also used to contain
;    autoionization rates and/or two photon transitions. These are
;    needed for accurately modelling the level populations within
;    ions, but do not yield photons of a definite wavelength
;    (autoionization yields no photons; two photon transitions yield a
;    continuum which is separately modelled in CHIANTI).
;
;    These "radiationless" transitions are denoted in the .wgfa file
;    with a zero wavelength. 
;
;    In the output structure all radiative decay rates are stored in
;    the .AVAL tag, while autoionization and two photon rates are
;    stored in the .AUTO tag.
;
; INPUTS
;
;    WGFANAME  The name of the file to be read.
;
; OUTPUTS
;
;    STR   A structure with the following tags:
;          .lvl1  The index of the lower level of the transition. 
;          .lvl2  The index of the upper level of the transition.
;          .wvl   The wavelength of the transition (angstroms).
;          .gf    The weighted oscillator strength.
;          .aval  The radiative decay rate (s^-1).
;          .auto  The autoionization rate (s^-1).
;
;    REF   A string array containing the file references.
;
; HISTORY
;
;    Ver.1, 13-Feb-2009, Peter Young
;-

chck=file_search(wgfaname)
IF chck EQ '' THEN BEGIN
  print,'%READ_WGFA_STR: file not found. Returning...'
  return
ENDIF

openr,lin,wgfaname,/get_lun

ss={lvl1: 0, lvl2: 0, wvl: 0., gf: 0., aval: 0., auto: 0.}
str=0

str1=''
tst1=0
WHILE tst1 EQ 0 DO BEGIN
  readf,lin,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    tst1=1
  ENDIF ELSE BEGIN
    reads,str1,format='(2i5,f15.0,2e15.0)',lvl1,lvl2,wvl,gf,aval
    IF wvl EQ 0. THEN BEGIN
      auto=aval 
      aval=0.
    ENDIF ELSE BEGIN
      auto=0.
    ENDELSE      
    ss.lvl1=lvl1
    ss.lvl2=lvl2
    ss.wvl=wvl
    ss.gf=gf
    ss.aval=aval
    ss.auto=auto
   ;
    IF n_tags(str) EQ 0 THEN BEGIN
      str=ss
    ENDIF ELSE BEGIN
      k=where(ss.lvl1 EQ str.lvl1 AND ss.lvl2 EQ str.lvl2,nk)
     ;
     ; If a transition already exists in STR, then it means that the previous
     ; A-value was a genuine radiative decay and the new one is an
     ; autoionization rate, or vice versa. This can be checked by looking
     ; at WVL - a zero wavelength indicates the autoionization rate. If
     ; both wavelengths are non-zero there's an error!
     ;
      IF nk NE 0 THEN BEGIN
        IF str[k].wvl NE 0. AND ss.wvl NE 0. THEN BEGIN
          print,'%READ_WGFA_STR: a transition is duplicated. Returning...'
          print,' Transition: ',trim(ss.lvl1),' - ',trim(ss.lvl2)
          free_lun,lin
          return
        ENDIF
        IF ss.wvl EQ 0. THEN str[k].auto=ss.auto
        IF ss.wvl NE 0. THEN BEGIN
          str[k].wvl=ss.wvl
          str[k].gf=ss.gf
          str[k].aval=ss.aval
        ENDIF
      ENDIF ELSE BEGIN
        str=[str,ss]
      ENDELSE
    ENDELSE
  ENDELSE
ENDWHILE

tst1=0
ref=''
WHILE tst1 EQ 0 DO BEGIN
  readf,lin,str1
  IF trim(str1) EQ '-1' THEN BEGIN
    tst1=1
  ENDIF ELSE BEGIN
    ref=[ref,str1]
  ENDELSE
ENDWHILE
ref=ref[1:*]

free_lun,lin

END
