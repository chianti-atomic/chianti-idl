

pro read_klgfb,pe,gfb,n

;+
; NAME:
;     READ_KLGB
;
; PURPOSE:
;     Return the Karzas & Latter (1961) freebound Gaunt factors for
;     the specified n-levels.
;
; CATEGORY:
;     CHIANTI; continuum; freebound.
;
; CALLING SEQUENCE:
;     READ_KLGFB, Pe, Gfb, N
;
; INPUTS:
;     N:   Principal quantum number. An integer between 1 and 6.
;
; KEYWORD PARAMETERS:
;	KEY1:	Document keyword parameters like this. Note that the keyword
;		is shown in ALL CAPS!
;
; OUTPUTS:
;     Pe:  The outgoing photon energy in Rydbergs (1D array).
;     Gfb:  Free bound Gaunt factors for the orbital sub-shells of the
;           specified principal quantum number. The array has
;           dimensions [NL,NE] where NL is the number of sub-shells
;           (=N-1), and NE is the number of energies.
;
; EXAMPLE:
;     Read the Gaunt factors for the 3s, 3p and 3d levels:
;
;     IDL> read_klgfb, pe, klgfb, 3
;
; MODIFICATION HISTORY:
;     Ver.1, 12-Feb-2021, Peter Young
;       A complete rewrite of the original read_klgfb to read new data
;       files. 
;-


IF n_params() LT 3 THEN BEGIN
   print,'Use:  IDL> read_klgfb, pe, klgfb, n'
   print,''
   print,'   n is an input with value between 1 and 6'
   return
ENDIF 

dir=concat_dir(!xuvtop,'continuum')
file=concat_dir(dir,'klgfb_'+trim(n)+'.dat')

chck=file_info(file)
IF chck.exists EQ 0 THEN BEGIN
   print,'% READ_KLGFB: The Karzas data file does not exist. Returning...'
   return
ENDIF 

;
; Get number of lines in file.
;
r=query_ascii(file,info)
nlines=info.lines

;
; Read data into array.
;
openr,lin,file,/get_lun
data=fltarr(n+1,nlines)
readf,lin,data
free_lun,lin

;
; Reformat data.
;
pe=reform(data[0,*])
gfb=reform(data[1:n,*])

IF n GT 1 THEN gfb=transpose(gfb)

END
