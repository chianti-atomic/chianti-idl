
pro read_ioneq,filename,t,ioneq,ref, sngl_ion=sngl_ion, element=element


;+
; NAME:
;     READ_IONEQ
;
; PURPOSE:
;     Reads a CHIANTI format ionization equilibrium (ioneq) file. 
;
; CATEGORY:
;     CHIANTI; read.
;
; CALLING SEQUENCE:
;     READ_IONEQ, Filename, Temperature, Ioneq, Ref
;
; INPUTS:
;     Filename:  The name of the file to be read.
;
; OPTIONAL INPUTS:
;     Sngl_Ion:  A string or string array containing a CHIANTI ion
;                name (or names). If specified, then the output array
;                will only contain non-zero values for those ions
;                specified by SNGL_ION.
;     Element:   Either the atomic number of an element, or the
;                symbol name (e.g., 'Fe' for iron). If set, then the
;                output array will only contain non-zero values for
;                the specified element.
;	
; OUTPUTS:
;     T:     Log10 temperature array (units: K) at which ionization
;            fractions are tabulated.
;     Ioneq: A 3D array (T,element,ion) containing ionization fractions.
;
; EXAMPLE:
;     IDL> read_ioneq, !ioneq_file, temp, ioneq
;     IDL> read_ioneq, !ioneq_file, temp, ioneq, sngl_ion=['fe_11','fe_12']
;     IDL> read_ioneq, !ioneq_file, temp, ioneq, element='fe'
;     IDL> read_ioneq, !ioneq_file, temp, ioneq, element=26
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere (KPD)
;	March 1996:     Version 2.0
;       March 1999:     KPD to read both number of temperature and number 
;                       of elements
;
;       25-Oct-2000     V. 4, Giulio Del Zanna (GDZ).
;                       Corrected to interpret the '-1' as a reference only
;                       if within the first 3 columns
;
;       V.  5, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       v.6, 25-Oct-2004, Peter Young
;            modified format statement so that it will read any number of
;            temperatures.
;
;       V 7, 25-May-2005, GDZ 
;                  corrected routine header.
;       V 7.1, 14-Apr-2013, Ken Dere
;                  minor change to name 'chianti for ioneq file
;
;       Ver.8, 26-Apr-2019, Peter Young
;           Added SNGL_ION and ELEMENT optional inputs; tidied up
;           header.
;
;       Ver.9, 10-Mar-2020, Peter Young
;           Fixed bug whereby dielectronic ions (specified by
;           sngl_ion) were not assigned the correct ion fractions. 
;
;       Ver.10, 27-Mar-2024, RPD
;           Increased length of comment array for advanced models
;-




if n_params(0) lt 3 then begin
   print,''
   print,' >  read_ioneq,filename,t,ioneq,ref [, sngl_ion=, element= ] '
   print,'      or'
   print,' > read_ioneq,''',''' , t,ioneq, ref [, sngl_ion=, element= ] '
   print,''
   return
endif
;
if strtrim(filename,2) eq '' then begin
;
    dir=concat_dir(!xuvtop,'ioneq')
    filename=dialog_pickfile(path=dir,filter='*.ioneq',title='Select Ionization Equilibrium File')
    print,' selected:  ',filename
endif
;
;
;
;
if file_test(filename) then begin
	openr,lu,filename,/get_lun
endif else begin
	print, ' file not found - ',filename
	return
endelse
;
string1=' '
;
str=''
nt=1
nz=1
;
readf,lu,str  ; read number of temperatures and elements
str=strtrim(str,2)
if strlen(str) le 3 then begin
   reads,str,nt    ;  v1.0 style
   ioneq=fltarr(nt,30,31)
endif else begin
    reads,str,nt,nz   ; v2.0 style
    ioneq=fltarr(nt,nz,nz+1)
endelse
   ioneq=fltarr(nt,30,31)   ;  hard-wired for now
;
z1=1 & ion1=1. & f1=fltarr(nt) 
;
readf,lu,f1
;
t=f1
;
;
ioneqform='(2i3,'+trim(nt)+'e10.2)'
;
;
while strpos(string1,'-1') EQ  -1 or strpos(string1,'-1') GT 2  do begin
readf,lu,string1
if(strpos(string1,'-1')   EQ  -1 or strpos(string1,'-1') GT 2) then begin

  reads,string1,z1,ion1,f1,format=ioneqform
  ioneq(0,z1-1,ion1-1)=f1(*)
endif
endwhile
;
;  get references
refstring=strarr(500)
nref=0
;
string1=' '
while strpos(string1,'-1') EQ  -1  do begin
readf,lu,string1


if(strpos(string1,'-1') EQ -1) and (strpos(string1,'%file') EQ -1) then begin
  refstring(nref)=string1
  nref=nref+1
endif
endwhile
;
ref=refstring(0:nref-1)
;
;
;g=where(ioneq gt 0.)
;ioneq(g)=10.^(-ioneq(g))
;
;
free_lun,lu

;
; If SNGL_ION is set, then set all elements of ioneq to zero apart
; from the ions in SNGL_ION
; 
IF n_elements(sngl_ion) NE 0 THEN BEGIN
  n=n_elements(sngl_ion)
  ioneq_save=ioneq
  ioneq=ioneq*0.
  FOR i=0,n-1 DO BEGIN 
    convertname,sngl_ion[i],iz,ion,dielectronic=dielectronic
    IF dielectronic EQ 1 THEN BEGIN
      ioneq[*,iz-1,ion]=ioneq_save[*,iz-1,ion]
    ENDIF ELSE BEGIN 
      ioneq[*,iz-1,ion-1]=ioneq_save[*,iz-1,ion-1]
    ENDELSE 
  ENDFOR 
  ioneq_save=0.
ENDIF


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
  ioneq_save=ioneq
  ioneq=ioneq*0.
  ioneq[*,elt_iz-1,*]=ioneq_save[*,elt_iz-1,*]
  ioneq_save=0.
ENDIF 


end
