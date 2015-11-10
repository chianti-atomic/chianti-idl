;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
;                   
; Name        : CALC_DMM_DR
;               
; Purpose     : Calculates CHIANTI density sensitive line ratios.
;               
; Explanation : 
;               
; Use         : Called from DMM_NE
;    
; Inputs      :  ciz  - element number
;                cion - ion number
;                denmin, denmax - log density range
;                mess_win - message widget ID for DMM_NE use.
;               
; Opt. Inputs :  None
;               
; Outputs     :  None
;               
; Opt. Outputs:  None
;               
; Keywords    :  Temperature - if set use this temperature rather than 
;                              temp at max of ioneq.
;
; Calls       :  EMISS_CALC
;
; Common      :
;               
; Restrictions:  None
;               
; Side effects:  None
;               
; Category    :  Spectral
;               
; Prev. Hist. :  Based on density_ratios by K Dere.
;
; Written     :  C D Pike, RAL, 22-Jan-96
;               
; Modified    :  Send warning message to widget.  CDP, 27-Jan-96
;                Update generally.  CDP, 7-Jun-97
;                Corrected typo in reading elvlc file.  CDP, 14-Jul-97
;                Added minimum line ratio factor.  CDP, 17-Jul-97
;                Added multiple wavelength ranges. CDP, 18-Jul-97
;                Fix typo for XXII stage.          CDP, 05-Mar-98
;                Update list of elements.          CDP, 18-Jun-99
;                v. 9. Correct hiccup in above.          CDP, 13-Jul-99
;
;                Update for compatibility with CHIANTI v.3. Replaced populate
;                with pop_solver. Left out 'dielectronic' lines. ouput
;                temperature. plus various small additions.
;                Now ratios are calculated with  0.1 increment in log Ne
;                Giulio Del Zanna (GDZ) 10-Oct-2000 
;
;                Ver.11, 7-Dec-01, Peter Young (PRY)
;                modified for v.4 of CHIANTI 
;
;                Ver. 12, 1-May-02, GDZ
;                Fixed a bug:  Added nlist = 0 when no lines are present. 
;                V. 13, 21-May-2002, GDZ 
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;                V.14, 06-Aug-02 GDZ
;                   Changed the use of CHIANTI system variables. 
;
;               V.15, 16-Sep-2004, GDZ
;                    added /NO_DE in call to emiss_calc as chianti_te expects
;                    emissivities in photon units
;                 
; Version     :  Version 15, 16-Sep-2004
;-                             

pro calc_dmm_dr,ciz,cion,denmin,denmax,mess_win,$
        temperature=temperature

;
;  denmin= log 10 of minimum density to be considered
;  denmax= log 10 of maximum density to be considered
;

;these COMMONS are used by pop_solver:

;COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth ,eref
;COMMON wgfa, wvl,gf,a_value
;COMMON upsilon,t_type,deu,c_ups,splups
;

;common with the caller routine, chianti_ne.pro :

common dmm_lines, list_wvl, list_int, list_descr1, list_descr2, species,$
  density, nden, ratio, description, nlist, savetext,$
  temper, intens

;commmon with the read_* routines:
common dmm_refs, ioneq_ref, wgfaref, copy_eref, upsref

;common with the caller routine, chianti_ne.pro :

common wavmm, wminw1, wminw2, wminw3, wminw4, $
  wmaxw1, wmaxw2, wmaxw3, wmaxw4, $
  minw1, minw2, minw3, minw4, $
  maxw1, maxw2, maxw3, maxw4

defsysv,'!xuvtop', EXISTS = EXISTS 
IF NOT EXISTS THEN $
message, 'system variable !xuvtop must be set using  '
xuvtop = !xuvtop


element=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
         'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti',$
         'V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII',$
          'XXIII','XXIV','XXV','XXVI','XXVII','XXVIII']


iz = where(strlowcase(element) eq strlowcase(ciz))
if iz(0) eq -1 then begin
   bell
   if n_elements(mess_win) gt 0 then begin
      widget_control, mess_win, /append, set_v='Unrecognised element'
   endif else begin 
      print,'Unrecognised element' 
   endelse
   return
endif else iz = iz(0)+1

ion = where(ionstage eq strupcase(cion))
if ion(0) eq -1 then begin
   bell
   if n_elements(mess_win) gt 0 then begin
      widget_control, mess_win, /append, set_v='Unrecognised ionization stage'   
   endif else begin 
      print,'Unrecognised ionization stage'   
   endelse
   return
endif else ion = ion(0)+1

;
;  Use T at max ioneq or user supplied
;
if temperature gt 0 then begin
   tmax = temperature
endif else begin

ioneqdir=concat_dir(!xuvtop, 'ioneq')
ioneq_name=!ioneq_file

   read_ioneq,ioneq_name,ioneq_t,ioneq,ioneq_ref

;dlnt=alog(10.^(ioneq_t(1)-ioneq_t(0)))
;

;   this_ioneq=ioneq(*,iz-1,ion-1+dielectronic)

   this_ioneq=ioneq(*,iz-1,ion-1)
   hit=where(this_ioneq eq max(this_ioneq))
   tmax=ioneq_t(hit)
;
   if n_elements(mess_win) gt 0 then begin
      widget_control, mess_win, $
        set_val='  deriving maximum T from the ionization eq. file: '+$
        !ioneq_file
   END 
   print,  ' deriving maximum T from the ionization eq. file: '+$
     !ioneq_file

   tmax=10.^(tmax(0))

;define temperature to pass on to caller
   temperature = tmax

endelse

;tmax=10.^(tmax(0))

if n_elements(mess_win) gt 0 then begin
   widget_control, mess_win, $
     set_val= string('Using Tmax= ',tmax,format='(a12,e10.3)'), /append
ENDIF
print,'Using Tmax=',tmax,format='(a12,e10.3)'


species=element(iz-1)+' '+ionstage(ion-1)


dlogden=0.1                    ;  increment in log Ne for calculating ratios

nden=fix((denmax-denmin)/dlogden) +1
density = 10.^(denmin+findgen(nden)*dlogden)

em=emiss_calc(iz,ion,diel=dielectronic,temp=alog10(tmax), $
              dens=alog10(density), noprot=noprot, $
              /quiet,rphot=rphot,radtemp=radtemp,/no_de)

index=where( ((em.lambda GE minw1) AND (em.lambda LE maxw1)) OR $
             ((em.lambda GE minw2) AND (em.lambda LE maxw2)) OR $
             ((em.lambda GE minw3) AND (em.lambda LE maxw3)) OR $
             ((em.lambda GE minw4) AND (em.lambda LE maxw4)) AND $
             (em.lambda NE 0.))


;
;  If no lines found then exit
;
if index[0] EQ -1 then begin
   bell
nlist = 0
   if n_elements(mess_win) gt 0 then begin
      widget_control, mess_win, set_val='No lines found in that region'
   endif else begin
      print,'No lines found in that region'
   endelse
   return
endif

nlist=n_elements(index)


list_descr1=em[index].lvl1_desc
list_descr2=em[index].lvl2_desc
list_int=fltarr(nlist,nden)
list_wvl=em[index].lambda

FOR ilist=0,nlist-1 DO $
     list_int[ilist,*]=reform(em[index[ilist]].em)/density

dd = list_int(*,nden-1)
ddmax = max(dd)
n = where(dd ge intens*ddmax)
list_wvl = list_wvl(n)
list_int = list_int(n,*)
list_descr1 = list_descr1(n)
list_descr2 = list_descr2(n)
nlist = n_elements(list_wvl)



;  print results

dash=' - '
;
widget_control, mess_win, /append,$
  set_val='----------------------------------------'
print,'------------------------------------------'

widget_control, mess_win, set_val= '#    Wavelength  Intensity    Transition ',/append
widget_control, mess_win, set_val= '        (A)   ',/append
widget_control, mess_win, set_val= ' ',/append

print, '#    Wavelength  Intensity    Transition '
print, '        (A)   '
print, ''

savetext = strarr(nlist)

for i=0,nlist-1 do begin
   text = string(i,list_wvl(i),list_int(i,nden-1),list_descr1(i),dash,$
                 strpad(list_descr2(i),15,/after), $
                 format='(i4,f10.4,e10.2,a15,a3,a15)')
   savetext(i) = text
   widget_control, mess_win, set_val=text,/append
   print,text
endfor
end
