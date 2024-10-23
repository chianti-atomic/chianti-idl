
function max_temp, ion, all=all, ioneq_name=ioneq_name, quiet=quiet

;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; Name        : MAX_TEMP
;               
; Purpose     : Calculates temperature at max ionisation ratio for an ion
;
; Explanation : 
;               
; Use         : IDL> print,max_temp(ion)
;    
; Inputs      :  ion - the specific ion in the form eg 'Fe_XII' or 'Fe XII'
;               
; Opt. Inputs :  None
;               
; Outputs     :  Function returns log of max temperature
;               
; Opt. Outputs:  None
;               
; Keywords    :  ALL - if set produces a plot of all temperatures
;                IONEQ_NAME:  	Name of the ionization equilization name to use.
;                    	    	If not passed, then the CHIANTI
;                               default, chianti.ioneq, is used.
;   	    	 QUIET:     	quiet mode.
;
; Calls       :  Other CHIANTI routines
;
; Common      :
;               
; Restrictions:  None
;               
; Side effects:  None
;               
; Category    :  Spectral
;               
; Prev. Hist. :  None
;
; Written     :  C D Pike, RAL, 9-Jun-97
;               
; Modified    :  V.2 Update element list.  CDP, 18-Jun-99 
;                V.3 modified definition of XUVTOP, and allowed selection of
;                ionization eq. file. Giulio Del Zanna (DAMTP) 10-Oct-2000
;
;                V.4, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;                V.5, I. Ugarte-Urra (IUU) 5-Jan-2009
;   	    	    added IONEQ_NAME and QUIET keywords.
;
;                V.6, Peter Young
;                   converted from DOS to UNIX format; now uses
;                   chianti.ioneq by default (rather than select with
;                   widget); got rid of check on xuvtop; updated
;                   format of table (if /all set)
;
; VERSION     :   6, 11-Mar-2020, Peter Young 
;
;-                             

element=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
         'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti',$
         'V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII',$
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI']




IF n_elements(ioneq_name) EQ 0 THEN BEGIN
  ioneq_name=!ioneq_file
  IF NOT keyword_set(quiet) THEN print,'% MAX_TEMP: using '+ioneq_name
ENDIF ELSE BEGIN
  chck=file_search(ioneq_name,count=count)
  IF count EQ 0 THEN BEGIN
    print,'% MAX_TEMP: the specified IONEQ_NAME could not be found. Returning...'
    return,-1
  ENDIF
ENDELSE
read_ioneq,ioneq_name,ioneq_t,ioneq,ioneq_ref



if not keyword_set(quiet) then begin
    print, ''
    FOR i=0, n_elements(ioneq_ref)-1 DO print, ioneq_ref(i)
    print, ''
endif
;
;  overall plot or individual?
;
if not keyword_set(all) then begin

   if n_params() lt 1 then begin
      print,'Use: IDL> print,max_temp(ion)'
      return,-1
   endif

;
;  parse input
;
   p = str2arr(ion,'_')
   if n_elements(p) eq 2 then begin
      ciz = p(0)
      cion = p(1)
   endif else begin
      p = str2arr(ion,' ')
      if n_elements(p) eq 2 then begin
         ciz = p(0)
         cion = p(1)
      endif else begin
         print,'Unparsable input - try eg Fe_XII or Fe XII.'
         return,-1
      endelse
   endelse
   
   iz = where(strlowcase(element) eq strlowcase(ciz))
   if iz(0) eq -1 then begin
      bell
      print,'Unrecognised element' 
      return,-1
   endif else iz = iz(0)+1

   ion = where(ionstage eq strupcase(cion))
   if ion(0) eq -1 then begin
      bell
      print,'Unrecognised ionization stage'   
      return,-1
   endif else ion = ion(0)+1
endif

;
;  work out the temperature
;
if keyword_set(all) then BEGIN

   if !d.name eq 'X' then window,(!d.window > 0),xs=1200,ys=650
   plot,[0,32],[-1,33],/nodata,xst=5,yst=5,/noclip
   s = size(ioneq)
   xyouts,9,33,'CHIANTI - Tmax using '+ioneq_name,chars=1

   for i=0,n_elements(element)-1 do begin
      oplot_string,[i+1,i+1],[0,0],element(i),chars=1.3
      oplot_string,[i+1,i+1],[32,32],element(i),chars=1.3
   endfor
   for j=0,n_elements(ionstage)-1 do begin
      oplot_string,[0,0],[j+1,j+1],ionstage(j),chars=1.3
      oplot_string,[31,31],[j+1,j+1],ionstage(j),chars=1.3
   endfor
   for i=0,s(2)-1 do begin
      for j=0,i+1 do begin
         this_ioneq=ioneq(*,i,j)
         hit=where(this_ioneq eq max(this_ioneq))
         if n_elements(hit) gt 1 then BEGIN
           temp=string(format='(f4.2)',ioneq_t[hit[0]])
          endif else BEGIN
            temp=string(format='(f4.2)',ioneq_t[hit[0]])
;            temp = trim(ioneq_t(hit(0)))
;            if strlen(temp) eq 1 then temp = temp + '.0'
         endelse
         oplot_string,[i+1,i+1],[j+1,j+1],temp,chars=1.3
      endfor
   endfor
endif else BEGIN

   this_ioneq=ioneq(*,iz-1,ion-1)
   hit=where(this_ioneq eq max(this_ioneq))

   if n_elements(hit) gt 1 then begin
      print,'multiple values found... using the first'
      bell
      hit = hit(0) 
   endif
   return,ioneq_t(hit)
endelse

end
