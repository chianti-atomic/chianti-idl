
PRO compare_ups_scups, ion_name, missing_trans=missing_trans,  $
                       prob_trans=prob_trans, auto=auto, $
                       perc_check=perc_check

;+
; NAME:
;     COMPARE_UPS_SCUPS
;
; PURPOSE:
;     Compare the upsilons de-scaled from the CHIANTI SCUPS file with
;     the original upsilons in the UPS file. If one of the descaled
;     upsilons is discrepant by 0.5% or more than the transition is
;     flagged as a problem transition.
;
; CATEGORY:
;     CHIANTI; data checking.
;
; CALLING SEQUENCE:
;     COMPARE_UPS_SCUPS, ion_name
;
; INPUTS:
;     Ion_Name:  The name of an ion in CHIANTI format (e.g., 'fe_13'
;                for Fe XIII). The routine will looks for the .ups and
;                .scups file in the current working directory.
;	
; OPTIONAL INPUTS:
;     Perc_Check: The routine defines a problem transition by a
;                 de-scaled upsilons being 0.5% discrepant with the
;                 original upsilon. This keyword allows you to change
;                 the percentage. For example, perc_check=1 use 1%
;                 instead. 
;	
; KEYWORD PARAMETERS:
;     AUTO:   If set, then the routine will check the .scups_auto file
;             instead of the .scups file.
;
; OUTPUTS:
;     Some text will be printed to the IDL window to state if there
;     are any problem or missing transitions.
;
; OPTIONAL OUTPUTS:
;     Missing_Trans:  A structure containing a list of the transitions
;                     that are missing in the SCUPS file. The tags are:
;                     .lvl1   The lower level
;                     .lvl2   The upper level
;     Prob_Trans:     A structure containing a list of the problem
;                     transitions. The tags are:
;                     .lvl1   The lower level
;                     .lvl2   The upper level
;                     .scup_ind  The index within scup structure.
;                     .frac_diff The max fractional difference between
;                                descaled and actual upsilon
;                     .ttype  Transition type
;                     .c_val  Scaling parameter
;
; EXAMPLES:
;      IDL> compare_ups_scups, 'mg_5'
;      IDL> compare_ups_scups, 'mg_5', /auto
;      IDL> compare_ups_scups, 'mg_5', missing_trans=missing_trans
;      IDL> compare_ups_scups, 'mg_5', prob_trans=prob_trans
;
; MODIFICATION HISTORY:
;     Ver.1, 3-Feb-2017, Peter Young
;     Ver.2, 23-Jun-2020, Peter Young
;       Added PERC_CHECK optional input.
;     Ver.3, 28-Jul-2020, Peter Young
;       Prints a message saying what perc_check is.
;     Ver.4, 04-Feb-2021, Peter Young
;       Added checks as to whether the files exist.
;-


IF n_params() LT 1 THEN BEGIN
  print,'Use: IDL> compare_ups_scups, ionname [, /auto, missing_trans=, prob_trans=, perc_check=]'
  return
ENDIF 


;
; Delete pre-existing structures
;
junk=temporary(missing_trans)
junk=temporary(prob_trans)

ups_file=ion_name+'.ups'
IF keyword_set(auto) THEN scups_file=ion_name+'.scups_auto' ELSE scups_file=ion_name+'.scups'
chck=file_info(scups_file)
IF chck.exists EQ 0 THEN BEGIN
   print,'% COMPARE_UPS_SCUPS: The SCUPS file was not found. Returning...'
   return
ENDIF 

chck=file_info(ups_file)
IF chck.exists EQ 0 THEN BEGIN
   print,'% COMPARE_UPS_SCUPS: The UPS file was not found. Returning...'
   return
ENDIF 


IF n_elements(perc_check) EQ 0 THEN perc_check=0.5


read_ups,ups_file,upsstr
read_scups,scups_file,scupstr

nt=upsstr.info.ntrans

str={lvl1: 0, lvl2: 0}
str2={lvl1: 0, lvl2: 0, scup_ind: 0l, frac_diff: 0., ttype: 0, c_val: 0.}

FOR i=0,nt-1 DO BEGIN
  t=upsstr.data[i].temp
  ups=upsstr.data[i].ups
  j=upsstr.data[i].lvl1
  k=upsstr.data[i].lvl2
 ;
  ii=where(scupstr.data.lvl1 EQ j AND scupstr.data.lvl2 EQ k,nii)
  IF nii EQ 0 THEN BEGIN
    str.lvl1=j
    str.lvl2=k
    IF n_tags(missing_trans) EQ 0 THEN missing_trans=str ELSE missing_trans=[missing_trans,str]
  ENDIF ELSE BEGIN
    descale_scups,t,scupstr,ii[0],ups_descale
    ichck=where((ups_descale-ups)/ups GE perc_check/100.,nchck)
    IF nchck GT 0 THEN BEGIN
      str2.lvl1=j
      str2.lvl2=k
      str2.scup_ind=ii[0]
      str2.ttype=scupstr.data[ii[0]].t_type
      str2.c_val=scupstr.data[ii[0]].c_ups
      str2.frac_diff=max(abs((ups_descale-ups)/ups))
      IF n_tags(prob_trans) EQ 0 THEN prob_trans=str2 ELSE prob_trans=[prob_trans,str2]
    ENDIF 
  ENDELSE 
ENDFOR

IF n_tags(missing_trans) GT 0 THEN BEGIN
  n=n_elements(missing_trans)
  print,'% COMPARE_UPS_SCUPS: There are '+trim(n)+' missing transitions in the SCUPS file.'
  print,'                     Give optional output MISSING_TRANS= to retrieve the list.'
ENDIF ELSE BEGIN
  print,'% COMPARE_UPS_SCUPS: There are no missing transitions in the SCUPS file.'
ENDELSE 

print,'% COMPARE_UPS_SCUPS: Checking for upsilons discrepant by '+trim(string(format='(f8.2)',perc_check))+'%...'

IF n_tags(prob_trans) GT 0 THEN BEGIN
  n=n_elements(prob_trans)
  print,'% COMPARE_UPS_SCUPS: ...there are '+trim(n)+' problem transitions in the SCUPS file.'
  print,'                     Give optional output PROB_TRANS= to retrieve the list.'
ENDIF ELSE BEGIN
  print,'% COMPARE_UPS_SCUPS: ...there are no problem transitions in the SCUPS file.'
ENDELSE 



END
