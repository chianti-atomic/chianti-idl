;+
;
; PROJECT:  CHIANTI
;
;        This program was developed as part of CHIANTI-VIP. 
;        CHIANTI-VIP (Virtual IDL and Python) is  a member of the CHIANTI family
;        mantained by Giulio Del Zanna, to develop additional features and
;        provide them to the astrophysics community.
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the
;       University of Cambridge, Goddard Space Flight Center, and University of Michigan.
;
;
; NAME:
;
;	CH_READ_LIST_IONS
;
;
; PURPOSE:
;
;	Read information about ions stored in a file. Adaption of read_masterlist
;       that also allows for comments at the bottom of the file and will read
;       extra information about each ion, in this case the number of levels to be
;       included in the advanced models when solving the ion balance for each ion
;
;
; CALLING SEQUENCE:
;
;       mlist=CH_READ_LIST_ions(filename)
;
;
; INPUTS:
;
;	None.
;
; OPTIONAL INPUTS:
;
;       filename:   Name of the file to be read. If not specified then the
;                   default masterlist.ions file will be read.
;	
; KEYWORD PARAMETERS:
;
;	ADVANCED:  If set, then reads the default advanced model masterlist file.
;
;
; OUTPUTS:
;
;	out:   a structure with
;              mlist:  string with the ion name
;              nlevels: the number of levels to include in the advanced models. For
;                       the regular masterlist, the values will all be -1.
;
;
; OPTIONAL OUTPUTS:
;
;       Count:  An integer giving the number of ions in the masterlist.
;
;
; PREVIOUS HISTORY:
;
;       Modified from read_masterlist
;
;
; WRITTEN:
;
;       Giulio Del Zanna (GDZ)
;       DAMTP, University of Cambridge, 11 Oct 2023
;
;
; MODIFICATION HISTORY:
;
;       Ver.2, Peter Young, 11-Oct-2024
;         If filename is not input, then the routine reads the default masterlist
;         file ; the /advanced keyword has now been implemented; added count=
;         optional output.
;
;
; VERSION     : 2
;
;- 

function ch_read_list_ions, filename, advanced=advanced, count=count
;

count=0

IF n_elements(filename) EQ 0 THEN BEGIN
  IF keyword_set(advanced) THEN BEGIN
    dir=concat_dir(!xuvtop,'ancillary_data')
    dir=concat_dir(dir,'advanced_models')
    filename=concat_dir(dir,'advmodel_list.ions')
  ENDIF ELSE BEGIN 
    dir=concat_dir(!xuvtop,'masterlist')
    filename=concat_dir(dir,'masterlist.ions')
  ENDELSE 
ENDIF 


if not file_exist(filename) then begin 
   print,' ERROR, input file does not exist'
       return,-1
    endif


;
mlist=''
nlevels=-1

;
openr,lum,filename,/get_lun
;
str=''
elstage=''
read_data=1
;
;   main input and calculation loop  **************
;
while not eof(lum) do begin
;
;   read the ion name and the number of levels 
;

   readf,lum, str
   
   if trim(str) eq  '-1' then read_data=0

   if read_data then begin
      
     index=strpos(str,';') ;  to remove comments
  ;
     if index ge 0 then str=trim(strmid(str,0,index-1))

     pp=str_sep(str,' ',/trim)
     g=where(pp ne '')
     pp=pp[g]
     
     mlist=[mlist, pp[0]]

     if n_elements(pp) gt 1 then nlevels=[nlevels, fix(pp[1])] else nlevels=[nlevels, -1]
     
   endif
;
endwhile
;
free_lun,lum
;
mlist=mlist[1:*]
nlevels=nlevels[1:*]

out={list_ions:mlist, nlevels:nlevels}

count=n_elements(mlist)

return, out 


end
;


