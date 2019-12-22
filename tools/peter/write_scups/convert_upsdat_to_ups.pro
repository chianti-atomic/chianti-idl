

PRO convert_upsdat_to_ups, infile, outfile

;+
; NAME
;
;      CONVERT_UPSDAT_TO_UPS
;
; PROJECT
;
;      CHIANTI
;
; EXPLANATION
;
;      Takes an old-style .upsdat file and writes it a new-style
;      .ups file.
;
; INPUTS
;
;      INFILE   The name of the .upsdat file to be read.
;
;      OUTFILE  The name of the output .ups file.
;
; OUTPUTS
;
;      The .ups file (OUTFILE) is written to the user's
;      working directory. 
;
; HISTORY
;
;      Ver.1, 20-May-2013, Peter Young
;-


read_upsdat_str, infile, upsdat, ref

n=n_elements(upsdat)

openw,lout,outfile,/get_lun

nt=n_elements(upsdat[0].temp)

FOR i=0,n-1 DO BEGIN
 ;
  IF upsdat[i].gf GT 0 THEN BEGIN
    lim=4.*upsdat[i].gf/upsdat[i].de
    format_str='(2i7,3e12.3,i5)'
  ENDIF ELSE BEGIN
    lim=-1
    format_str='(2i7,2e12.3,i12,i5)'
  ENDELSE 
 ;
  printf,lout,format=format_str,upsdat[i].lvl1,upsdat[i].lvl2, $
         upsdat[i].de,upsdat[i].gf,lim,nt
 ;
  format_str='('+trim(nt)+'e12.3)'
  printf,lout,format=format_str,upsdat[i].temp
 ;
  printf,lout,format=format_str,upsdat[i].ups
ENDFOR 

printf,lout,' -1'

nref=n_elements(ref)
FOR i=0,nref-1 DO BEGIN
  printf,lout,ref[i]
ENDFOR 

free_lun,lout

END
