;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;                   
; Name        : PLOT_IONEQ
;               
; Purpose     : Plots the ionisation equilibrium values for an element.
;
; Explanation : 
;               
; Use         : IDL> plot_ioneq, element [ ion=ion]
;    
; Inputs      : element - the element name
;               
; Opt. Inputs : Ion (as keyword)
;               
; Outputs     : None
;               
; Opt. Outputs: a postscript file.
;               
; Keywords    : ION_RANGE - specify range of ions. E.g., ion_range=[5,8] 
;                    means V to VIII inclusive.
;               IONEQ_NAME:  Name of the ionization equilization name to use.
;                    If not passed, then the user is prompted for it.
;               NOT_INTERACTIVE Avoid interactive use.
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
; Modified    :  V.2. Update element list. modified definition of XUVTOP, and
;               allowed selection of ionization eq. file and creation of
;               postscript file. 
;               Giulio Del Zanna  (DAMTP) 10-Oct-2000 
;
;       V.3, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V 4, 15-July-2002, GDZ 
;         Added keywords ioneq_name, not_interactive
;
;       V.5, 9-Feb-2005, Peter Young
;            Changed ion= keyword to ion_range=
;
; VERSION     :    5, 9-Feb-2005
;
;
;-                             

pro plot_ioneq, elem,  ion_range=ion_range, ioneq_name=ioneq_name, $
                not_interactive=not_interactive

element=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
         'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti',$
         'V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII',$
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII','XXIX']

;
;  enough input?
;
if n_params() eq 0 then begin
   print,'Use:  IDL> plot_ioneq, element [, ion_range=, ioneq_name=, '
   print,'                /not_interactive ]'
   return
endif

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
message, 'system variable !xuvtop must be set '
xuvtop = !xuvtop

IF n_elements(ioneq_name) EQ 0 THEN BEGIN 
dir=concat_dir(!xuvtop,'ioneq')
ioneq_name=pickfile(path=dir,filter='*.ioneq', $
                       title='Select Ionization Equilibrium File')
END 

ff = findfile(ioneq_name)
IF  ff(0)  NE ''  THEN $
read_ioneq,ioneq_name,ioneq_t,ioneq,ioneq_ref ELSE BEGIN 
   print, 'Error,  no ioneq file found !'
   return
END 


print, ''
FOR i=0, n_elements(ioneq_ref)-1 DO print, ioneq_ref(i)
print, ''

iz = where(strlowcase(element) eq strlowcase(elem))
if iz(0) eq -1 then begin
   bell
   print,'Unrecognised element' 
   return
endif else iz = iz(0)

;
;  which ions to plot
;
if n_elements(ion_range) gt 0 then begin
   min_ion = min(ion_range)-1
   max_ion = max(ion_range)-1
endif else begin
   min_ion = 0
   max_ion = n_elements(ionstage)-1
endelse

;print,min_ion,max_ion

 break_file,ioneq_name, disk,dir,f,ext

;tit = '!6Ionization equilibrium !3for !17'+element(iz)+'!3'
tit = 'Ionization equilibrium for '+element(iz)+' -  '+f


window, 0
circle_sym,/fill

x_min=min(ioneq_t) 
x_max=max(ioneq_t)
y_min=0.0
y_max=1.2

pr=''
begin_post:

!p.multi=0


plot,ioneq_t,ioneq(*,iz,min_ion),psym=-8,$
   xrange=[x_min,x_max], yrange=[y_min,y_max],$
    xstyle=1,ystyle=1,$
     tit=tit,$
     chars=1.5,xtit='log ( T / K )', ytit='Ion fraction'

hit = where(ioneq(*,iz,min_ion) eq max(ioneq(*,iz,min_ion)))
t = ioneq_t(hit(0))
xyouts,t+0.02,max(ioneq(*,iz,min_ion))+0.02,'!17'+ionstage(min_ion)+'!3',$
        chars=1.5

if max_ion ne min_ion then begin
   for i = min_ion+1, max_ion do begin
      if i mod 2 eq 1 then circle_sym else circle_sym,/fill
      oplot,ioneq_t,ioneq(*,iz,i),psym=-8
      hit = where(ioneq(*,iz,i) eq max(ioneq(*,iz,i)))
      if max(ioneq(*,iz,i)) gt 0 then begin
         t = ioneq_t(last_item(hit)) + 0.02
         if t lt x_max and t gt x_min then begin
            xyouts,t,max(ioneq(*,iz,i))+0.02,'!17'+ionstage(i)+'!3',chars=1.5
         endif
      endif
   endfor
ENDIF 

IF NOT keyword_set(not_interactive) THEN BEGIN 
print2d_plot, pr=pr , x_min=x_min,x_max=x_max,y_min=y_min,y_max=y_max,$
          go_to_line=go_to_line, out_name=out_name , /ask_name

if go_to_line eq 'y' then goto,begin_post
END 

end
