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
;	READ_IONEQ
;
; PURPOSE:
;
;	to read files containing the ionization equilibrium values
;
; CATEGORY:
;
;	science.
;
; CALLING SEQUENCE:
;
;       READ_IONEQ, File, T, Ioneq, Ref
;
;
; INPUTS:
;
;	File:	for example, !xuvtop+'/ioneq/arnaud_rothenflug.ioneq'
;
; OPTIONAL INPUTS:
;
;	None:
;	
;
; OUTPUTS:
;
;	T:  array of log10 temperatures
;       Ioneq:  3-d array (T,element,ion) 
;               of the fractional abundance of the ion in
;              ionization equilibrium.
;       Ref:  reference in the scientific literature
;
;
; EXAMPLE:
;
;             > read_ioneq,!xuvtop+'/ioneq/arnaud_rothenflug.ioneq'
;             > help,t,ioneq
;             > T               FLOAT     = Array(41)
;             > IONEQ           FLOAT     = Array(41, 28, 29)
;             > print, minmax(t)
;             >  4.00000      8.00000
;             > print,t(20)
;             >  6.0000
;             > print,ioneq(20,25,9)
;             >   0.269                  = fractional abundance of  Fe X in
;                                          ionization equilibrium 
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
;
; VERSION     :   7, 25-May-2005
;
;-
pro read_ioneq,filename,t,ioneq,ref
;
;
if n_params(0) lt 4 then begin
   print,''
   print,' >  read_ioneq,filename,t,ioneq,ref '
   print,'      or'
   print,' > read_ioneq,''',''' , t,ioneq, ref'
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
refstring=strarr(100)
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
;
end
