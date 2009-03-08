;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an atomic database package for the calculation of
;       continuum and emission line spectra from astrophysical plasmas. It is a 
;       collaborative project involving the Naval Research Laboratory
;       (Washington DC, USA), the Arcetri Observatory (Firenze, Italy), and the
;       Cambridge University (United Kingdom).
;
;
; NAME:
;	READ_MASTERLIST
;
; PURPOSE:
;
;	read the masterlist.ions types of file and output a list of ions
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_MASTERLIST,filename,list
;
;
; INPUTS:
;
;	filename:   name of the masterlist file
;	
;
;	
; KEYWORD PARAMETERS:
;
;	none
;
; OUTPUTS:
;
;	list:  list of ions
;
;
;
; COMMON BLOCKS:
;
;	none;
;
;
; EXAMPLE:
;
;             > read_masterlist,'!xuvtop/masterlist.masterlist.ions',list
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	December 1998:  first version  
;
;-
pro read_masterlist,filename,list
;
;
if n_params(0) lt 2 then begin
  print,' '
  print,'  IDL>  read_masterlist,filename,list
  print,'           or'
  print,'  IDL>  read_masterlist,''',''',list  for a choice
  print,' '
  return
endif
;
if filename eq '' then begin
   path=!xuvtop+'/masterlist/'
   file='masterlist.ions'
   filename=dialog_pickfile(path=path,file=file,title='Select a masterlist file')
   if filename eq '' then begin
       print,' Did not select a masterlist.ions file
       return
   endif
endif
;
list=''
;
openr,lum,filename,/get_lun
;
gname=''
elstage=''
;
;   main input and calculation loop  **************
;
while not eof(lum) do begin
;
;   read the name of the ions unless a single ion (sngl_ion) has been specified
;
if (not keyword_set(sngl_ion)) then begin
   readf,lum,gname
   index=strpos(gname,';') ;  to sort out comments
;
   if index ge 0 then gname=strmid(gname,0,index-1)
   gname=strtrim(gname,2)
   list=[list,gname]
endif 
;
endwhile
;
free_lun,lum
;
list=list(1:*)
end
;


