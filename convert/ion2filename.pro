;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME:
;	ion2filename
;
; PURPOSE:
;	Ion names as character strings are converted to
;	provide their complete file name (with out suffix)
;
; CATEGORY:
;	
;	naming utility
;
; CALLING SEQUENCE:
;
;       ION@FILENAME,Ion,Filename
;
;
; INPUTS:
;	Name:   such as 'c_2'
;
;
; OUTPUTS:
;
;	Filename:  !xuvtop/c/c_2/c_2
;
;
; EXAMPLE:
;
;                     > ion2filename,'c_2d',filename
;                     > print,filename
;                     > !xuvtop/c/c_2d/c_2d
;
; MODIFICATION HISTORY:
;
; 	Written by:	Ken Dere
;
;	V.2,  September 1999:  added for use with Version 3
;
;       V.3, 21-May-2002, Giulio Del Zanna (GDZ) 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
; VERSION     :   3, 21-May-2002 
;
;
;-
pro ion2filename,ion,filename
;
if n_params(0) lt 2 then begin
   print,'    use> ion2filename,ion,filename'
   print,'    use> ion2filename, "c_6",filename'
   print,'           giving filename = /data/c/c_6 , for example'
   return
endif
;
;
name=strtrim(ion,2)
;
;
locname=strlowcase(name)
pos=strpos(locname,'_')
l=strlen(pos)
first=strmid(locname,0,pos)
last=strmid(locname,pos+1,l-pos-1)
;
if pos ge 0 then begin
;
;
   name=ion

dir=concat_dir(!xuvtop, first)

filename=concat_dir(concat_dir(dir, name), name)

;
endif else print,' incorrect format for ion name (i.e., c_2) ',ion
;
return
end
