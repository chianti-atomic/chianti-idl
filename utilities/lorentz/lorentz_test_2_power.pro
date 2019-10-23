
PRO lorentz_test_2_power, loss_rate

;+
; NAME:
;     LORENTZ_TEST_2_POWER
;
; PURPOSE:
;     Computes the CHIANTI output for Lorentz test No. 2 (Power). 
;
; CATEGORY:
;     CHIANTI; Lorentz test.
;
; CALLING SEQUENCE:
;     LORENTZ_TEST_2_POWER
;
; INPUTS:
;     None.
;
; OUTPUTS:
;     Creates the file 'lorentz_test_2_chianti.txt' in the current
;     working directory.
;
; EXAMPLE:
;     IDL> lorentz_test_2_power
;
; MODIFICATION HISTORY:
;     Ver.1, 1-May-2019, Peter Young
;     Ver.2, 8-May-2019, Peter Young
;        Call to proton_dens was missing the /hydrogen keyword, so
;        fixed this now.
;-

basename='lorentz_test_2_chianti.txt'

nt=51
logt=findgen(nt)/10.+4.0

abundfile=concat_dir(!xuvtop,'abundance')
abundfile=concat_dir(abundfile,'proto_solar_2009_lodders.abund')

;
; Power is computed for an electron number density of 1 cm^-3.
;
n_e=1.0

;
; Set plasma volume in cm^3
;
v=1e6

;
; Compute NH/Ne.
;
nh_NE=proton_dens(logt,abund_file=abundfile,ioneq_file=!ioneq_file,/hydrogen)

openw,lout,basename,/get_lun

;
; Write the header.
;
printf,lout,'#Column 1: Atomic number of element'
printf,lout,'#Columns 2-52: Total radiative power in erg s^-1 from a 1 m^3 volume with electron density 10^6 m^-3, tabulated for log temperature = 4.0 to 9.0 at 0.1 dex intervals'
printf,lout,'#Derived using CHIANTI version '+ch_get_version()
printf,lout,'#Abundance file: '+file_basename(abundfile)
printf,lout,'#File created: '+systime()

;
; Write first row.
;
printf,lout,format='(a3,51f10.1)','Z',logt

;
; Now write the radiative losses for each element up to zinc.
;
FOR i=0,29 DO BEGIN
  bb_rad_loss,t,r_bb,element=i+1,abund_file=abundfile,density=n_e
  fb_rad_loss,t,r_fb,element=i+1,abund_file=abundfile
  ff_rad_loss,t,r_ff,element=i+1,abund_file=abundfile
  r=r_bb+r_fb+r_ff
  r_out=dblarr(nt)
 ;
 ; r_out is tabulated for more temperatures than we need, so here we
 ; pull out the required values.
 ;
  FOR j=0,nt-1 DO BEGIN
    getmin=min(abs(logt[j]-alog10(t)),imin)
    r_out[j]=r[imin]
  ENDFOR
 ;
 ; Multiply by V*NH*Ne
 ;
  r_out=r_out*v*n_e^2*nh_NE
 ;
  z2element,i+1,elt,/symbol
  elt=strpad(elt,3,/after,fill=' ')
  printf,lout,format='(i3,51e10.3)',i+1,r_out
ENDFOR 


free_lun,lout

END
