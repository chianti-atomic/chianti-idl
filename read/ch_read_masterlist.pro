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
;	CH_READ_MASTERLIST
;
; PURPOSE:
;
;	read information about ions stored in a file. 
;
;
; CALLING SEQUENCE:
;
;       mlist=CH_READ_MASTERLIST(filename)
;
;
; INPUTS:
;
;	filename:   name of the file
;	
; KEYWORD PARAMETERS:
;
;	none
;
; OUTPUTS:
;
;	mlist:  string with the information about ions
;
; Prev. Hist. : modified from read_masterlist
;
; Written     : Roger Dufresne (RPD) and Giulio Del Zanna (GDZ)
;               DAMTP, University of Cambridge
;
; Modified    :
;
; VERSION     : 1, 16 Sept 2023 
;
;- 

function ch_read_masterlist,filename
;

;
if filename eq '' then begin
   print,' ERROR, file not defined'
       return,-1
    endif

if not file_exist(filename) then begin 
   print,' ERROR, input file does not exist'
       return,-1
    endif


;
mlist=''
;
openr,lum,filename,/get_lun
;
gname=''
elstage=''
read_data=1
;
;   main input and calculation loop  **************
;
while not eof(lum) do begin
;
;   read the name of the ions unless a single ion (sngl_ion) has been specified
;
if (not keyword_set(sngl_ion)) then begin
   readf,lum,gname
   
   ; changed from here for tr_models
   if trim(gname) eq '-1' then read_data=0
   
   if read_data eq 1 then begin
     index=strpos(gname,';') ;  to sort out comments
  ;
     if index ge 0 then gname=strmid(gname,0,index-1)
     gname=strtrim(gname,2)
     mlist=[mlist,gname]
   endif
   ; end of changes
endif 
;
endwhile
;
free_lun,lum
;
mlist=mlist[1:*]

return, mlist


end
;


