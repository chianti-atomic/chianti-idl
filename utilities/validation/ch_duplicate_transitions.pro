
FUNCTION ch_duplicate_transitions, filename, status=status, quiet=quiet

;+
; NAME:
;     CH_DUPLICATE_TRANSITIONS
;
; PURPOSE:
;     Checks the specified CHIANTI file to see if there are any transitions
;     that are duplicated in the file.
;
; CATEGORY:
;     CHIANTI; validity check.
;
; CALLING SEQUENCE:
;     Result = CH_DUPLICATE_TRANSITIONS( Filename )
;
; INPUTS:
;     Filename:  The name of a CHIANTI file to check. 
;
; KEYWORD PARAMETERS:
;     QUIET:  If set, then no information messages are printed to IDL window.
;
; OUTPUTS:
;     If status=1, then the output is a structure array with the following
;     tags:
;      .lvl1  Lower level of duplicate transition.
;      .lvl2  Upper level of duplicate transition.
;      .n_trans  No. of transitions in the file.
;
;     For other values of status, the output is the number -1.
;
; OPTIONAL OUTPUTS:
;     Status:  A byte scalar with the value
;              0 - file exists, and there are no duplicates.
;              1 - file exists, and there are duplicates.
;              2 - file does not exist
;              3 - file exists, but the file is not processed by this routine.
;
; EXAMPLE:
;     IDL> output=ch_duplicate_transitions('c_5.wgfa')
;
; MODIFICATION HISTORY:
;     Ver.1, 06-Feb-2025, Peter Young
;     Ver.2, 22-Aug-2025, Peter Young
;      Now checks if the filename contains the filename extension rather than
;      if it ends with the filename extension.
;-


status=0b

chck=file_info(filename)
IF chck.exists EQ 0 THEN BEGIN
  IF NOT keyword_set(quiet) THEN message,/info,/cont,'The specified file does not exist. Returning...'
  status=2b
  return,-1
ENDIF 

CASE 1 OF
  filename.contains('.wgfa'): BEGIN
    read_wgfa_str,filename,wgfa
    lvl1=wgfa.lvl1
    lvl2=wgfa.lvl2
  END
  filename.contains('.scups'): BEGIN
    read_scups,filename,scups
    lvl1=scups.data.lvl1
    lvl2=scups.data.lvl2
  END
  filename.contains('.rrlvl'): BEGIN
    rrdata=read_rrlvl(filename.remove(-6),status_rr)
    lvl1=rrdata.initial_level
    lvl2=rrdata.final_level
  END
  ELSE: BEGIN
    message,/info,/cont,'No matches to the CHIANTI file types. Returning...'
    status=3b
    return,-1
  END 
ENDCASE


lvl1_str=trim(lvl1)
lvl1_str=lvl1_str.reverse()
lvl1_str=lvl1_str.insert('0',6,fill='0')
lvl1_str=lvl1_str.substring(0,3)
lvl1_str=lvl1_str.reverse()

lvl2_str=trim(lvl2)
lvl2_str=lvl2_str.reverse()
lvl2_str=lvl2_str.insert('0',6,fill='0')
lvl2_str=lvl2_str.substring(0,3)
lvl2_str=lvl2_str.reverse()

lvl_str=lvl1_str+lvl2_str
n=n_elements(lvl_str)
IF NOT keyword_set(quiet) THEN message,/info,/cont,'There are '+trim(n)+' transitions in the file.'

a=lvl_str[uniq(lvl_str,sort(lvl_str))]
n2=n_elements(a)
IF NOT keyword_set(quiet) THEN message,/info,/cont,'There are '+trim(n2)+' unique transitions in the file.'

IF n NE n2 THEN BEGIN
  status=1b
  str={ lvl1: 0, lvl2: 0, n_trans: 0}
  FOR j=0,n2-1 DO BEGIN
    k=where(lvl_str EQ a[j],nk)
    IF nk GT 1 THEN BEGIN
      str.lvl1=fix(a[j].substring(0,3))
      str.lvl2=fix(a[j].substring(4,7))
      str.n_trans=nk
      IF n_tags(output) EQ 0 THEN output=str ELSE output=[output,str]
    ENDIF 
  ENDFOR 
ENDIF ELSE BEGIN
  output=-1
ENDELSE 

return,output

END
