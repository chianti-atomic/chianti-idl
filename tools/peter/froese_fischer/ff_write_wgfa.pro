
PRO ff_write_wgfa, trans, levstr, outfile, maxlev=maxlev

;+
; NAME:
;    FF_WRITE_WGFA
;
; PURPOSE:
;    Takes the complete set of Froese Fischer radiative data for an
;    ion and level map structure, and writes out a CHIANTI .wgfa
;    file.
;
; INPUTS:
;    TRANS:  A structure containing the FF data (see the routine
;            ff_read_data.pro).
;    LEVSTR: A structure that maps the FF level notation to a CHIANTI
;            level index (see ff_read_levels.pro).
;    OUTFILE:Name of the CHIANTI .wgfa file to write to.
;
; OPTIONAL INPUTS:
;    Maxlev: Only transitions with level indices less than this
;            number will be printed to the wgfa file.
;
; OUTPUTS:
;    Sends the radiative data to OUTFILE in the CHIANTI .wgfa format.
;
; MODIFICATION HISTORY:
;    Ver.1, 4-Oct-2009, Peter Young
;    Ver.2, 9-Jun-2017, Peter Young
;       Added check for multiple matches to same level.
;    Ver.3, 12-Jun-2017, Peter Young
;       Now accepts new format for levstr structure.
;    Ver.4, 26-Jun-2025, Peter Young
;       Added maxlev= optional input.
;-

IF n_params() LT 3 THEN BEGIN
  print,'Use:  IDL> ff_write_wgfa, trans, levstr, outfile [ maxlev= ]'
  return
ENDIF 

IF n_elements(maxlev) EQ 0 THEN maxlev=max(levstr.ind)

openw,lout,outfile,/get_lun

n=n_elements(trans)

FOR i=0,n-1 DO BEGIN
  k1=where(trim(trans[i].conf1) EQ trim(levstr.conf) AND $
           trim(trans[i].lev1) EQ trim(levstr.lev) AND $
           round(trans[i].e1_cm) EQ levstr.energy,nk1)
  k2=where(trim(trans[i].conf2) EQ trim(levstr.conf) AND $
           trim(trans[i].lev2) EQ trim(levstr.lev) AND $
           round(trans[i].e2_cm) EQ levstr.energy,nk2)
 ;
  IF nk1 GT 1 THEN BEGIN
    print,'There are multiple matches for level '+trans[i].conf1+' '+trans[i].lev1
  ENDIF
 ;
  IF nk2 GT 1 THEN BEGIN
    print,'There are multiple matches for level '+trans[i].conf2+' '+trans[i].lev2
  ENDIF 
 ;
  IF nk1 EQ 0 AND nk2 EQ 0 THEN BEGIN
    print,'Transition '+trim(i+1)+'/'+trim(n)
    print,' - both levels not found in level map'
  ENDIF
 ;
  IF nk1 EQ 0 AND nk2 EQ 1 THEN BEGIN
    print,'Transtion '+trim(i+1)+'/'+trim(n)
    print,' - lower level not found: '+trim(trans[i].conf1)+' '+trim(trans[i].lev1)
  ENDIF
 ;
  IF nk1 EQ 1 AND nk2 eq 0 THEN BEGIN
    print,'Transtion '+trim(i+1)+'/'+trim(n)
    print,' - upper level not found: '+trim(trans[i].conf2)+' '+trim(trans[i].lev2)
  ENDIF
 ;
  IF nk1 EQ 1 AND nk2 EQ 1 THEN BEGIN
    lev1=levstr[k1].ind
    lev2=levstr[k2].ind
    wvl=1d8/abs(trans[i].e1_cm-trans[i].e2_cm)
    IF lev1 LE maxlev AND lev2 LE maxlev THEN BEGIN 
      printf,lout,format='(2i5,f15.3,2e15.3)',lev1,lev2,wvl,trans[i].gf, $
             trans[i].aval
    ENDIF 
  ENDIF
ENDFOR

printf,lout,' -1'
printf,lout,'%file: '
printf,lout,'%Prepared for the CHIANTI database by Peter Young'
printf,lout,' -1'

free_lun,lout

END
