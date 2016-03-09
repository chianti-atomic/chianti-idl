
PRO convert_splups_to_scups, infile, outfile

;+
; NAME
;
;    CONVERT_SPLUPS_TO_SCUPS
;
; PROJECT
;
;    CHIANTI
; 
; EXPLANATION
;
;    Converts the old style SPLUPS file to the new format.
;
; INPUTS
;
;    INFILE   Name of the .splups file to convert.
;
;    OUTFILE  Name of the output .scups file.
;
; OUTPUTS
;
;    Creates the new .scups file in the current working directory..
;
; HISTORY
;
;    Ver.1, 20-May-2013, Peter Young
;    Ver.2, 31-May-2013, Peter Young
;       changed name to reflect new name for the output file.
;-


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

  IF gf NE 0. THEN BEGIN
    lim=4.*gf/de 
    format_str='(2i7,3e12.3,i3,e12.3,i5)'
  ENDIF ELSE BEGIN
    lim=-1.
    format_str='(2i7,2e12.3,i12,i3,e12.3,i5)'
  ENDELSE 

  printf,lout,format=format_str,l1,l2,de,gf,lim,t_type,c_val,nt
        
  format_str='('+trim(nt)+'e12.3)'
  printf,lout,format=format_str,st
  printf,lout,format=format_str,sups

ENDFOR 

printf,lout,' -1'

nref=n_elements(ref)
FOR i=0,nref-1 DO BEGIN
  IF trim(ref[i]) NE '' THEN printf,lout,ref[i]
ENDFOR 

free_lun,lout

END
