
;+
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics 
;       of Astrophysical Plasmas. It is a collaborative project involving 
;       the Naval Research Laboratory (USA), the University of Florence 
;       (Italy), the University of Cambridge and the Rutherford Appleton 
;       Laboratory (UK). 
;
; NAME:
;
;       FB_RAD_LOSS
;
; PURPOSE:
;
;       Calculate the total radiative losses of a plasma due to the 
;       free-bound (radiative recombination) continuum.
;
; EXPLANATION
;
;       This routine does not use the same method of calculating the ion 
;       continuum emissivities as the FREEBOUND routine. This is because a 
;       modified version of FREEBOUND would be very slow for this purpose. 
;       Instead, we use the method of Mewe et al. (A&AS 65, 511, 1986) which 
;       is outlined in their Sect.2.2. The integration over the quantity 
;       P_c(lambda,T) is very simple as the gaunt factor, G_c, has no 
;       intrinsic lambda dependence other than through the limits.
;
;       Comparisons between the wavelength-resolved continuum emission 
;       derived from the Mewe et al. method with the more sophisticated 
;       method employed in FREEBOUND show excellent agreement with at most 
;       10% differences at specific wavelengths.
;
; CALLING SEQUENCE:
;
;      fb_rad_loss, temp, int, min_abund=min_abund, /no_setup
;
; OUTPUTS
;
;	TEMP    Temperatures (K). These are the temperatures at which the 
;               ion fractions are defined (typically 10^4 to 10^8 in 0.1 
;               dex intervals).
;
;       INT     The emissivity in units of erg cm^3 s^-1.
;
;
; OPTIONAL INPUTS:
;
;       MIN_ABUND Exclude elements whose abundances are less than MIN_ABUND. 
;                 (Note that Ab(H)=1.)
;	
; KEYWORD PARAMETERS:
;
;	NO_SETUP   If the procedure setup_elements has already been called 
;                  then the keyword /nosetup should be set to avoid 
;                  repeating this step
;
; PROGRAMMING NOTES
;
;       This routine computes the free-bound gaunt factor following the 
;       prescription set out in Sect.2.2 of Mewe et al. (1986). The 
;       expression for f_2 (Eq. 16 of Mewe et al.) contains several 
;       quantities that need to be computed. Zeta_0 is computed internally 
;       by the function ZETA_0 through a prescription evident from browsing 
;       Table I of Mewe et al. Z_0 is computed from the ionization potential 
;       of the recombined ion which is contained in the .ip file within the 
;       database. The quantity n_0 (also used in deriving Z_0) is derived 
;       using the routine CONF2N which extracts the highest n value from 
;       the configuration description of the ground term of each ion.
;
;       The quantity Z is just the charge on the recombined ion; E_0 is the 
;       ionization potential of the recombined ion, while E_n_0+1 is derived 
;       from Mewe's Eq.7 with the prescription that z_n=Z and n=n_0+1.
;
;       The ions that are considered for the continuum are those for which 
;       .fblvl files exist.
;
; COMMON BLOCKS:
;
;	common elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref
;
; INTERNAL FUNCTIONS
;
;       ZETA_0
;
; CALLS
;
;       READ_IP,  ZION2FILENAME, SETUP_ELEMENTS, FILE_EXIST,
;       READ_FBLVL, CONCAT_DIR, GET_IEQ
;
; EXAMPLES
;
;       IDL> fb_rad_loss,temp,int
;       IDL> plot,temp,int,/xlog,/ylog
;
; MODIFICATION HISTORY:
;
;     Ver.1, 1-Aug-2002, Peter Young
;          Completely new version of fb_rad_loss. Incorporates code from 
;          a version of freebound not available in CHIANTI.
;
;       V 2, 25-May-2005, GDZ 
;                  corrected routine header.
;
; VERSION:   2, 25-May-2005
;
;-

FUNCTION zeta_0, iz, ion

;+
; NAME 
;
;     ZETA_0
;
; EXPLANATION
;
;     Returns the value of zeta_0 (the number of vacancies in the ion given 
;     by IZ and ION). See Sect. 2.2 of Mewe et al. (1986, A&AS 65, 511).
;
; INPUTS
;
;     IZ    Atomic number of ion (e.g., 26 -> Fe)
;
;     ION   Spectroscopic number of ion (e.g., 13 -> XIII)
;
; OUTPUTS
;
;     Value of zeta_0.
;-

CASE 1 OF
  (iz-ion) GT 22: return,double(ion-iz+55)
  ((iz-ion) LE 22) AND ((iz-ion) GT 8): return,double(ion-iz+27)
  ((iz-ion) LE 8) AND ((iz-ion) GT 0): return,double(ion-iz+9)
  (iz-ion) LE 0: return,double(ion-iz+1)
ENDCASE

END


PRO fb_rad_loss, temp, int, min_abund=min_abund, $
                    no_setup=no_setup

COMMON elements,abund,abund_ref,ioneq,ioneq_t,ioneq_ref


IF n_params() LT 2 THEN BEGIN
  print,'Use:  IDL> fb_rad_loss, temp, int [, /no_setup, min_abund= ]'
  return
ENDIF

z_lbl=['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si',$
       'p','s','cl','ar','k','ca','sc','ti','v','cr','mn','fe','co','ni',$
      'cu','zn']

IF NOT keyword_set(no_setup) THEN setup_elements

read_ip,concat_dir(concat_dir(!xuvtop, 'ip'), 'chianti.ip'),ionpot,ipref

temp=double(10.^ioneq_t)

nt=n_elements(temp)

g_fb=dblarr(nt)
int=dblarr(nt)

temp6=temp/1d6

ehcm=ionpot[0,0]       ; H ionization potential

IF n_elements(min_abund) EQ 0 THEN min_abund=0.

FOR iz=1,30 DO BEGIN
  IF abund[iz-1] GE min_abund THEN BEGIN
  FOR ion=1,iz DO BEGIN
    ieq=get_ieq(temp,iz,ion+1,ioneq_logt=ioneq_t,ioneq_frac=ioneq)
   ;
   ; need to read in .elvlc file to get n for ground state of recombined ion 
   ;
    zion2filename,iz,ion,name
    name=name+'.fblvl'
    IF file_exist(name) EQ 0 THEN GOTO,lbl1
    read_fblvl,name,l1,conf,pqn,ll,spd,mult,ecm,ecmth,ref

   ;
   ; compute gaunt factor
   ;
    n_0=pqn[0]
    z_0=sqrt(ionpot[iz-1,ion-1]/ehcm)*n_0

    e_0=ionpot[iz-1,ion-1]/8065.54d0/1d3 ; convert to keV
    e_n1=ehcm*(ion)^2/(n_0+1)^2/8065.54d0/1d3

    l_0=1d8/(ionpot[iz-1,ion-1]-ecm[0])
    l_n1=1d8/ehcm/(ion)^2*(n_0+1)^2

    zeta0=zeta_0(iz,ion+1)

    ind_t=where(ieq NE 0.)
    IF ind_t[0] NE -1 THEN BEGIN 
      f2=dblarr(nt)
      f2[ind_t]=0.9d0*zeta0*z_0^4/n_0^5*exp(0.1578*z_0^2/n_0^2/temp6[ind_t])
      g_fb=(0.1578/temp6*abund[iz-1]*ieq*f2)

      int=int+2.051d-22*temp6/143.9*g_fb/sqrt(temp6)*exp(-143.9/l_0/temp6)
      
      f2[ind_t]=0.42d0*n_0^(-1.5)*double(ion)^4* $
           exp(0.1578*double(ion)^2/(n_0+1)^2/temp6[ind_t])
      g_fb=(0.1578/temp6*abund[iz-1]*ieq*f2)
      
      int=int+2.051d-22*temp6/143.9*g_fb/sqrt(temp6)*exp(-143.9/l_n1/temp6)
    ENDIF

    lbl1:
  ENDFOR
  ENDIF   
ENDFOR

END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
