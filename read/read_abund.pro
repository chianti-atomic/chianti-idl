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
; OPTIONAL INPUTS:
;	Element: Either the atomic number of an element, or the
;                symbol name (e.g., 'Fe' for iron). If set, then the
;                output array will only contain non-zero values for
;                the specified element.
;
; OUTPUTS:
;
;	Abundance:  an array of abuncance values
;       Reference:  a string containing the reference to the chosen set
;                   of abundances in the scientific literature
;
; EXAMPLE:
;
;             > read_abund,'allen.abund',abundance,ref
;               abundance(26) = abundance of iron relative to hydrogen
;               quoted by C.W. Allen in Astrophysical Quantities
;
;             > read_abund,!abund_file,abundance,ref,element='fe'
;             > read_abund,!abund_file,abundance,ref,element=26
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	March 1996:     Version 2.0
;
;       V.   3, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS.
;
;      Ver.4, 26-Apr-2019, Peter Young
;          Added ELEMENT optional input.
;-

PRO  read_abund,filename,abund,ref, element=element

;
;
if n_params(0) lt 3 then begin
   print,''
   print,' > read_abund, filename, abund, ref [, element= ]'
   print,'      or'
   print,' > read_abund,''',''' , abund, ref [, element= ]'
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
if file_test(filename) then begin
	openr,lua,filename,/get_lun
endif else begin
	print,' file not found - ',filename
	return
endelse
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
; Handle the input 'element'.
;
elt_iz=-1
IF n_elements(element) NE 0 THEN BEGIN
  IF datatype(element) EQ 'STR' THEN BEGIN
    z2element,indgen(30)+1,elt,/symbol,/lower_CASE
    k=where(strlowcase(element) EQ elt,nk)
    IF nk NE 0 THEN elt_iz=k[0]+1
  ENDIF ELSE BEGIN 
    elt_iz=element[0]
  ENDELSE
  ab_save=abund
  abund=abund*0.
  abund[elt_iz-1]=ab_save[elt_iz-1]
ENDIF 

end


































