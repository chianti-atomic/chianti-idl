;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas.  See www.chiantidatabase.org 
;
;     		          
; Purpose     : converts a CHIANTI ascii SCUPS file to a FITS table.
;
; Explanation : the  CHIANTI ascii SCUPS file is read into an IDL
;               structure, which is then written to a FITS table.
;    
; Inputs      : the file names
;
; Outputs     : none
;
; Calls:
;             read_scups and  wrt_fits_bin_exten
;
; 
; Written     : 
;       Version 1, Giulio Del Zanna (GDZ)  13 Oct 2020
;
; Modified    :
;
; VERSION     :    V.1, 13 Oct 2020
;
;-        
;---------------------------

pro ch_scups2fits, file_scups, file_fits

; read the structure with the data:
  read_scups,  file_scups,  st, /verb
  
  wrt_fits_bin_exten,st.info, file_fits

  wrt_fits_bin_exten,st.data, file_fits,/append
  

end
