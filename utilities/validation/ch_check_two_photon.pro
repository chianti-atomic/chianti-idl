
PRO ch_check_two_photon

;+
; NAME:
;     CH_CHECK_TWO_PHOTON
;
; PURPOSE:
;     This checks the two photon transitions in the H-like and He-like wgfa
;     files to make sure they are present.
;
; CATEGORY:
;     CHIANTI; validity.
;
; CALLING SEQUENCE:
;     CH_CHECK_TWO_PHOTON
;
; INPUTS:
;     None.
;
; OUTPUTS:
;     Text messages are printed to the IDL window. If no problems are found,
;     then a message stating this will be printed.
;
; EXAMPLE:
;     CH_CHECK_TWO_PHOTON
;
; MODIFICATION HISTORY:
;     Ver.1, 07-Feb-2025, Peter Young
;-


mlist=ch_read_list_ions(count=n)
mlist=mlist.list_ions

problems=0b

FOR i=0,n-1 DO BEGIN
  convertname,mlist[i],iz,ion,diel=diel
  ;
  ; H-like ions
  ;
  IF iz-ion EQ 0 AND diel EQ 0 THEN BEGIN
    zion2filename,iz,ion,fname
    read_elvlc,fname+'.elvlc',elvlc=elvlc
    read_wgfa_str,fname+'.wgfa',wgfa,two_photon=two_photon
    ;
    ; If the two_photon structure exists then we're OK. Otherwise check if
    ; the transition exists, but has a non-zero wavelength.
    ; In case the level index of the 2S1/2 levels changes with ion, then I
    ; check the elvlc file to get the index.
    ;
    IF n_tags(two_photon) EQ 0 THEN BEGIN
      problems=1b
      message,/info,/cont,'H-like: no two photon transition for '+mlist[i]+'.'
      d=elvlc.data
      k=where(d.level EQ '2S1/2' AND d.conf EQ '2s',nk)
      IF nk EQ 1 THEN BEGIN
        upp_lev=d[k].index
        k=where(wgfa.lvl1 EQ 1 AND wgfa.lvl2 EQ upp_lev,nk)
        CASE nk OF
          0: message,/info,/cont,'H-like: no transitions from 2s 2S1/2 level for '+mlist[i]+'.'
          1: BEGIN
            IF wgfa[k[0]].wvl NE 0. THEN message,/info,/cont,'H-like: wavelength for 2s 2S1/2 transition for '+mlist[i]+' is not zero!'
          END
          ELSE: message,/info,/cont,'H-like: more than one transition from 2s 2S1/2 level for '+mlist[i]+'.'
        ENDCASE 
      ENDIF ELSE BEGIN
        message,/info,/cont,'H-like: the 2s 2S1/2 level is missing for '+mlist[i]+'.'
      ENDELSE
    ENDIF 
  ENDIF
  ;
  ; H-like ions
  ;
  IF iz-ion EQ 1 AND diel EQ 0 THEN BEGIN
    zion2filename,iz,ion,fname
    read_elvlc,fname+'.elvlc',elvlc=elvlc
    read_wgfa_str,fname+'.wgfa',wgfa,two_photon=two_photon
    ;
    ; If the two_photon structure exists then we're OK. Otherwise check if
    ; the transition exists, but has a non-zero wavelength.
    ;
    IF n_tags(two_photon) EQ 0 THEN BEGIN
      problems=1b
      message,/info,/cont,'He-like: no two photon transition for '+mlist[i]+'.'
      d=elvlc.data
      k=where(d.level EQ '1S0' AND d.index LE 7,nk)
      IF nk EQ 1 THEN BEGIN
        upp_lev=d[k].index
        k=where(wgfa.lvl1 EQ 1 AND wgfa.lvl2 EQ upp_lev,nk)
        CASE nk OF
          0: message,/info,/cont,'He-like: no transitions from 2s 1S0 level for '+mlist[i]+'.'
          1: BEGIN
            IF wgfa[k[0]].wvl NE 0. THEN message,/info,/cont,'He-like: wavelength for 2s 1S0 transition for '+mlist[i]+' is not zero!'
          END
          ELSE: message,/info,/cont,'He-like: more than one transition from 2s 1S0 level for '+mlist[i]+'.'
        ENDCASE 
      ENDIF ELSE BEGIN
        message,/info,/cont,'He-like: the 2s 1S0 level is missing for '+mlist[i]+'.'
      ENDELSE
    ENDIF 
  ENDIF
  
ENDFOR


IF ~ problems THEN message,/info,/cont,'The two photon transitions are all present and correct!'

END
