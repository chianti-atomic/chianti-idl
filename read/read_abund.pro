;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;
; NAME:
;	READ_ABUND
;
; PURPOSE:
;
;	to read CHIANTI abundance files
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_ABUND,File,Abundance,Reference
;
;
; INPUTS:
;
;	File:  the name (string) of the file containing the abundance values
;                (relative to hydrogen) usually of the form 
;                 '!xuvtop/abundance/*.abund'
;
;
;	
;
; OUTPUTS:
;
;	Abundance:  an array of abuncance values
;       Reference:  a string containing the reference to the chosen set
;                   of abundances in the scientific literature
;
;
;
; PROCEDURE:
;
;	You can describe the foobar superfloatation method being used here.
;
; EXAMPLE:
;
;             > read_abund,'allen.abund',abundance,ref
;               abundance(26) = abundance of iron relative to hydrogen
;               quoted by C.W. Allen in Astrophysical Quantities
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       V.   3, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
; VERSION     :   3, 21-May-2002 
;
;-

PRO  read_abund,filename,abund,ref

;
;
if n_params(0) lt 3 then begin
   print,''
   print,' > read_abund, filename, abund, ref '
   print,'      or'
   print,' > read_abund,''',''' , abund, ref'
   print,''
   return
endif
;
if strtrim(filename,2) eq '' then begin
;
    dir=concat_dir(!xuvtop,'abundance')
    filename=dialog_pickfile(path=dir,filter='*.abund',title='Select Abundance File')
    print,' selected:  ',filename
endif
;
;
;
openr,lua,filename,/get_lun
;
abund=fltarr(50)
;
;
string1=' '
z1=1 & a1=1.
;
;
while strpos(string1,'-1') lt 0 or (strpos(string1,'-1') gt 3) do begin
readf,lua,string1
if (strpos(string1,'-1') lt 0) or (strpos(string1,'-1') gt 3) then begin
  reads,string1,z1,a1
  abund(z1-1)=a1
endif
endwhile
;
;  get references
refstring=strarr(100)
nref=0
;
string1=' '
while strpos(string1,'-1') lt 0  do begin
readf,lua,string1
if(strpos(string1,'-1') lt 0) and (strpos(string1,'%file') lt 0) then begin
  refstring(nref)=string1
  nref=nref+1
endif
endwhile
;
ref=refstring(0:nref-1)
;
;
g=where(abs(abund) ne 0.)
abund(g)=10.^(abund(g)-abund(0))
;
;
free_lun,lua
;
;
end


































