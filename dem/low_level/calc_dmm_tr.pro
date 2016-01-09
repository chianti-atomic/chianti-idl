;+
; Project     : SOHO - CDS     
;                   
; Name        : CALC_DMM_TR
;               
; Purpose     : Calculates CHIANTI temperature sensitive line ratios.
;               
; Explanation : 
;               
; Use         : Called from DMM_TE
;    
; Inputs      :  ciz  - element number
;                cion - ion number
;                wmin, wmax - wavelength range
;                temmin, temmax - log temperature range
;                mess_win - message widget ID for DMM_TE use.
;               
; Opt. Inputs :  None
;               
; Outputs     :  None
;               
; Opt. Outputs:  None
;               
; Keywords    :  None
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
; Prev. Hist. :  Based on temperature_ratios by K Dere.
;
; Written     :  C D Pike, RAL, 22-Jan-96,
; Modified    :  H E Mason, 03-Oct-96 (density to temperature)
;               
; Modified    :  Send warning message to widget.  CDP, 27-Jan-96
;                Update list of elements and cut call to 
;                       read_elvl.                CDP, 18-Jun-99
;                V.3. Fix hiccupp in above.            CDP, 13-Jul-99
;
;                V.4. Rewritten completely, adding possibility to pass on the
;                density at which the intensities are calculated, and making
;                this routine compatible ith CHIANTI v.3.
;                Giulio Del Zanna (DAMTP), 10 Oct-2000
;
;                Ver.5, 7-Dec-2001, Peter Young
;                Modified for v.4 of CHIANTI.
;
;                Ver. 6, 1-May-02, GDZ
;                Fixed a bug:  Added nlist = 0 when no lines are present. 
;
;               V.7, 16-Sep-2004, Peter Young
;                    added /NO_DE in call to emiss_calc as chianti_te expects
;                    emissivities in photon units
;
; Version     :  Version 7, 16-Sep-04
;-                             

pro calc_dmm_tr,ciz,cion,wmin,wmax,temmin,temmax,mess_win, $
    density=density

;

;these COMMONS are used by pop_solver:

;COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;COMMON wgfa, wvl,gf,a_value
;COMMON upsilon,t_type,deu,c_ups,splups
;

;common with the caller routine, chianti_te.pro :

common dmm_lines_te, list_wvl, list_int, list_descr1, list_descr2, species,$
                  temperature, ntem, ratio, description, nlist, savetext, $
                  dens, intens

common wavmm, wminw1, wminw2, wminw3, wminw4, $
              wmaxw1, wmaxw2, wmaxw3, wmaxw4, $
              minw1, minw2, minw3, minw4, $
              maxw1, maxw2, maxw3, maxw4

;commmon with the read_* routines:
common dmm_refs, ioneq_ref, wgfaref, copy_eref, upsref



element=['H','He','Li','Be','B','C','N','O','F','Ne','Na',$
         'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti',$
         'V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

;
ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII',$
          'XIII','XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI',' XXII',$
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


IF  density EQ  0 THEN BEGIN 
density=10.^8
ENDIF 

if n_elements(mess_win) gt 0 then begin
   widget_control, mess_win, $
      set_val='calculating intensities at constant Ne = '+$
      string(density,format='(e10.3)')
END 
;
print,'calculating intensities at constant Ne = '+string(density,format='(e10.3)')

species=element(iz-1)+' '+ionstage(ion-1)

;
dlogtem=0.1       ;  increment in log T for calculating ratios
;
ntem=fix((temmax-temmin)/dlogtem) +1
;

temperature=10.^(temmin+findgen(ntem) *dlogtem)

em=emiss_calc(iz,ion,diel=dielectronic,temp=alog10(temperature), $
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
 nlist=0
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
list_int=fltarr(nlist,ntem)
list_wvl=em[index].lambda

FOR ilist=0,nlist-1 DO $
     list_int[ilist,*]=reform(em[index[ilist]].em)/density

dd = list_int(*,ntem-1)
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
   text = string(i,list_wvl(i),list_int(i,ntem-1),list_descr1(i),dash,$
                 strpad(list_descr2(i),15,/after), $
                 format='(i4,f10.4,e10.2,a15,a3,a15)')
   savetext(i) = text
   widget_control, mess_win, set_val=text,/append
   print,text
endfor


end
