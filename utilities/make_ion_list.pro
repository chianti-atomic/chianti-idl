;+
; Project     : SOHO - CDS     
;                   
; Name        : MAKE_ION_LIST
;               
; Purpose     : Reads masters ions list and interprets.
;               
; Explanation : Reads the CHIANTI master ion list and returns the list of 
;               available elements and ions in suitable arrays.
;               
; Use         : IDL> make_ion_list, elements, ions
;    
; Inputs      : None
;               
; Opt. Inputs : None
;               
; Outputs     : elements - string array of all elements available
;               ions     - 2-d string array of roman ionization stages
;                          available for each element. (N,M) where N is
;                          number of elements and M is large enough to 
;                          contain all ionization stages (Null if not used)
;               
; Opt. Outputs: None
;               
; Keywords    : None
;
; Calls       : None
;
; Common      : None
;               
; Restrictions: None
;               
; Side effects: None
;               
; Category    : Spectrum
;               
; Prev. Hist. : None
;
; Written     : C D Pike, RAL, 19-Jan-96
;               
; Modified    : V.2 Increased array sizes.   CDP, 25-Jun-96
;               V.3. Further increase of the array sizes, and  removed the 'd'
;               directories from the masterlist, to be compatible with CHIANTI
;                removed call to concat_dir    Giulio Del Zanna 16-Oct-2000. 
;               V.4 21-May-2002, GDZ: 
;                  generalize directory concatenation to work for Unix, Windows
;                  and VMS.
; 
; Version     : Version 4 21-May-2002
;-            

pro make_ion_list, list_el,list_ions

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
message, 'Error, system variable !xuvtop must be set !'
xuvtop = !xuvtop


element=['','H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na',$
     'Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',$
     'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn']


ionlist=['','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
         'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII',$
         'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII']


;
;  read the master ions list (should be in the masterion subdirectory
;  of $CDS_SS_DERE
;
file = concat_dir(!xuvtop, 'masterlist')
file = concat_dir(file, 'masterlist.ions')
openr,lun, (file), /get_lun,error=err

if err ne 0 then begin
   print,'Unable to read masterlist file.'
   return
endif

;
;  read the file
;   
en = intarr(1000)
in = intarr(1000)
i = -1
while not eof(lun) do begin
   text = ''
   readf,lun,text
   text = str2arr(text,' ')
   text = text(0)

;exclude the 'd' directories
index=strpos(text,'d')

IF  index EQ -1 THEN BEGIN 
   i = i + 1
   convertname, text(0),a,b
   en(i) = a
   in(i) = b
ENDIF
endwhile
free_lun,lun

;
;  now we are all in numbers
;

n = where(en ne 0)
en = en(n)
in = in(n)

;
;  so they can be sorted
;
ns = sort(en)
en = en(ns)
in = in(ns)


;
;  get a list of elements used
;
n =rem_dup(en)
el = en(n)


;
;  set up storage for ions associated with each element
;
numel = n_elements(el)
ions = intarr(60,30)

list_el = strarr(numel)
for i=0,numel-1 do begin
   n = where(en eq el(i))
   ion = in(n)
   ion = ion(sort(ion))
   list_el(i) = element(el(i))
   ions(el(i),0:n_elements(ion)-1) = ion
endfor
;
;  trim the ions array to max needed
;
n = where(ions(*,0) ne 0)
ions = ions(n,*)

nmax = 0
for i=0,numel-1 do begin
   n = max(where(ions(i,*) gt 0))
   if n gt nmax then nmax = n
endfor

;
;  and convert back to Roman
;
ions =ions(*,0:nmax)
list_ions = ionlist(ions)

end
