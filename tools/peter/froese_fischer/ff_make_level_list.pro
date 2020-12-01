
PRO ff_make_level_list, trans, outfile=outfile

;+
; NAME:
;      FF_MAKE_LEVEL_LIST
;
; PURPOSE:
;      Given the levels identified in the TRANS structure (produced by
;      process_ff_data.pro), this routine prints an energy ordered list of
;      levels. 
;
; CATEGORY:
;      CHIANTI; data formatting.
;
; CALLING SEQUENCE:
;      FF_MAKE_LEVEL_LIST, TRANS
;
; INPUTS:
;      Trans:   A structure in the format produced by
;               process_ff_data.pro. 
;
; OPTIONAL INPUTS:
;      Outfile: The name of a file to send the output to. If not
;               specified, then it is sent to ff_lev_map.txt.
;
; OUTPUTS:
;      A list of levels is sent to OUTFILE. There are four columns:
;      configuration, level descriptor, level index and level energy
;      in cm^-1.
;
; EXAMPLE:
;      IDL> ff_read_data,'ar_6_ff.txt',trans
;      IDL> ff_make_level_list,trans,outfile='ar_6_ff_levs.txt'
;       -> may need to edit the ar_6_ff_levs.txt file at this point.
;      IDL> ff_write_wgfa, trans, levstr, 'ar_6.wgfa_ff'
;
; MODIFICATION HISTORY:
;      Ver.1, 11-Jun-2017, Peter Young
;-



str={conf: '', lev: '', index: 0, energy: 0l}
lvlstr=0

n=n_elements(trans)

count=0
FOR i=0,n-1 DO BEGIN
  str.conf=trans[i].conf1
  str.lev=trans[i].lev1
  str.energy=round(trans[i].e1_cm)
  IF n_tags(lvlstr) NE 0 THEN BEGIN
    k=where(str.conf EQ lvlstr.conf AND str.lev EQ lvlstr.lev AND str.energy EQ lvlstr.energy,nk)
    IF nk EQ 0 THEN BEGIN
      count=count+1
      str.index=count
      lvlstr=[lvlstr,str]
    ENDIF 
  ENDIF ELSE BEGIN
    count=count+1
    str.index=count
    lvlstr=str
  ENDELSE
 ;
  str.conf=trans[i].conf2
  str.lev=trans[i].lev2
  str.energy=round(trans[i].e2_cm)
  IF n_tags(lvlstr) NE 0 THEN BEGIN
    k=where(str.conf EQ lvlstr.conf AND str.lev EQ lvlstr.lev AND str.energy EQ lvlstr.energy,nk)
    IF nk EQ 0 THEN BEGIN
      count=count+1
      str.index=count
      lvlstr=[lvlstr,str]
    ENDIF 
  ENDIF ELSE BEGIN
    count=count+1
    str.index=count
    lvlstr=str
  ENDELSE
  
ENDFOR 

n=n_elements(lvlstr)
print,'***There are '+trim(n)+' unique levels'

IF n_elements(outfile) EQ 0 THEN outfile='ff_lev_map.txt'

openw,lout,outfile,/get_lun

k=sort(lvlstr.energy)
lvlstr=lvlstr[k]
lvlstr.index=indgen(n)+1

FOR i=0,n-1 DO BEGIN
  printf,lout,format='(a20,a10,i5,i12)',strpad(lvlstr[i].conf,20,/after), $
         strpad(lvlstr[i].lev,10,/after),lvlstr[i].index,lvlstr[i].energy
ENDFOR 

free_lun,lout


END
