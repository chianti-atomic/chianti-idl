
PRO write_all_ff_data, trans, levstr, outfile

;+
; NAME
;
;    WRITE_ALL_FF_DATA
;
; EXPLANATION
;
;    Takes the complete set of Froese Fischer radiative data for an
;    ion and level map structure, and writes out a CHIANTI .wgfa
;    file.
;
; INPUTS
;
;    TRANS   A structure containing the FF data (see the routine
;            process_ff_data.pro).
;
;    LEVSTR  A structure that maps the FF level notation to a CHIANTI
;            level index (see read_ff_map.pro).
;
;    OUTFILE Name of the CHIANTI .wgfa file to write to.
;
; OUTPUTS
;
;    Sends the radiative data to OUTFILE in the CHIANTI .wgfa format.
;
; HISTORY
;
;    Ver.1, 4-Oct-2009, Peter Young
;-

openw,lout,outfile,/get_lun

n=n_elements(trans)

FOR i=0,n-1 DO BEGIN
  k1=where(trim(trans[i].conf1) EQ trim(levstr.conf) AND trim(trans[i].lev1) EQ trim(levstr.lev),nk1)
  k2=where(trim(trans[i].conf2) EQ trim(levstr.conf) AND trim(trans[i].lev2) EQ trim(levstr.lev),nk2)
 ;
  IF nk1 EQ 0 AND nk2 EQ 0 THEN BEGIN
    print,'Transtion '+trim(i+1)+'/'+trim(n)
    print,' - both levels not found in level map'
  ENDIF
 ;
  IF nk1 EQ 0 AND nk2 NE 0 THEN BEGIN
    print,'Transtion '+trim(i+1)+'/'+trim(n)
    print,' - lower level not found: '+trim(trans[i].conf1)+' '+trim(trans[i].lev1)
  ENDIF
 ;
  IF nk1 NE 0 AND nk2 eq 0 THEN BEGIN
    print,'Transtion '+trim(i+1)+'/'+trim(n)
    print,' - upper level not found: '+trim(trans[i].conf2)+' '+trim(trans[i].lev2)
  ENDIF
 ;
  IF nk1 NE 0 AND nk2 NE 0 THEN BEGIN
    lev1=levstr[k1].ind
    lev2=levstr[k2].ind
    printf,lout,format='(2i5,f15.3,2e15.3)',lev1,lev2,0.,trans[i].gf, $
           trans[i].aval 
  ENDIF
ENDFOR

printf,lout,' -1'
printf,lout,'%file: '
printf,lout,'%Prepared for the CHIANTI database by Peter Young'
printf,lout,' -1'

free_lun,lout

END
