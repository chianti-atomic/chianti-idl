
PRO convert_splups_to_scups, infile, outfile

;+
; NAME:
;      CONVERT_SPLUPS_TO_SCUPS
;
; PURPOSE:
;      Converts the old style SPLUPS file to the new SCUPS format.
;
; CATEGORY:
;      CHIANTI; file conversion.
;
; CALLING SEQUENCE:
;      CONVERT_SPLUPS_TO_SCUPS, Infile, Outfile
;
; INPUTS:
;      Infile:  The name of a CHIANTI file in the ".SPLUPS" format.
;      Outfile: The name of the file to be written.
;
; OUTPUTS:
;      The file OUTFILE is written and has the format of the CHIANTI
;      .SCUPS files. For more information about the format, please
;      check CHIANTI Technical Report No. 13.
;
; EXAMPLE:
;      IDL> convert_splups_to_scups, 'o_6.splups', 'o_6.scups'
;
; CALLS:
;      READ_SPLUPS
;
; MODIFICATION HISTORY:
;      Ver.1, 20-May-2013, Peter Young
;      Ver.2, 31-May-2013, Peter Young
;      Ver.3, 5-Mar-2018, Peter Young
;         Updated to the final SCUPS format; updated header.
;-


IF n_params() LT 2 THEN BEGIN
  print,'Use: IDL> convert_splups_to_scups, infile, outfile'
  return
ENDIF 

read_splups,infile,splstr,ref

n=n_elements(splstr)

openw,lout,outfile,/get_lun


FOR i=0,n-1 DO BEGIN
  l1=splstr[i].lvl1
  l2=splstr[i].lvl2
  ttype=splstr[i].t_type
  gf=splstr[i].gf
  de=splstr[i].de
  c_val=splstr[i].c_ups
  t_type=splstr[i].t_type
  nt=splstr[i].nspl
  sups=splstr[i].spl[0:nt-1]

  st=findgen(nt)/float(nt-1)

  IF gf NE 0. THEN lim=4.*gf/de ELSE lim=-1
  
  IF c_val GE 1e4 OR c_val LT 0.1 THEN cform='e12.4' ELSE cform='f12.5'
  IF lim GE 0 THEN limform='(e12.3)' ELSE limform='(i12)'
  format_str='(2i7,2e12.3,'+limform+',2i5,'+cform+')'

  printf,lout,format=format_str,l1,l2,de,gf,lim,nt,t_type,c_val
        
  format_str='('+trim(nt)+'e12.3)'
  printf,lout,format=format_str,st
  printf,lout,format=format_str,sups

ENDFOR 

printf,lout,' -1'

nref=n_elements(ref)
FOR i=0,nref-1 DO BEGIN
  IF trim(ref[i]) NE '' THEN printf,lout,ref[i]
ENDFOR
printf,lout,' -1'

free_lun,lout

END
