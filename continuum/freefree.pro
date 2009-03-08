;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME
;
;     FREEFREE
;
;  PURPOSE:
;
;     This routine computes the free-free continuum (bremsstrahlung) 
;     using the fitting formulae of Itoh et al. (ApJS 128, 125, 2000) 
;     and Sutherland (MNRAS 300, 321, 1998).
;
;     The Itoh et al. data are valid for smaller ranges for temperature 
;     and wavelength than Sutherland and so for points outside of their 
;     ranges we use the data of Sutherland.
;
; INPUTS
;
;    TEMP    Temperature (in K). They do not have to be in monotonically
;            increasing order.
;
;    WVL     Wavelengths in angstroms. Can be a scalar or vector. If /keV is
;            set, then WVL is assumed to contain energies in keV units.
;
; OUTPUTS
;
;    INT     Free-free continuum intensity in units of 
;            10^-40 erg cm^3/s/sr/Angstrom  
;            [ integral(N_H N_e dh) in cm^-5 ] if a DEM is not defined. 
;
;            If DEM values are defined, it is assumed that they are given
;            as N_H N_e dh/dT.  The units are 10^-40 erg/cm^2/s/sr/Angstrom. 
;
;            If T is given as a 1-D array, then the output will be a 2-D array,
;            with one element for each temperature and wavelength 
;            (but also see SUMT).
;
;            If the keyword /keV is set, then the units of INT will be 
;            10^-40 erg cm^3/s/sr/keV
;
; OPTIONAL INPUTS
;
;    DEM_INT   This should be the same size as TEMP and contain
;              differential emission measure values. Specifying
;              DEM_INT places a restriction on TEMP (see above). Note
;              that DEM_INT values are multiplied by the factor dT
;              within the routine, whereas EM_INT values are not.
;
;    EM_INT    This should be the same size as TEMP and contain
;              emission measure values. Note
;              that DEM_INT values are multiplied by the factor dT
;              within the routine, whereas EM_INT values are not. If
;              neither EM_INT or DEM_INT are specified, then a EM_INT
;              vector of same size as TEMP with values of unity is
;              assumed.
;
;    MIN_ABUND This keyword allows the specification of a minimum abundance, 
;              such that any elements with an abundance (relative to 
;              hydrogen) less than MIN_ABUND will not be included in the 
;              calculation. E.g., MIN_ABUND=1e-5.
;
; KEYWORDS
;
;    NO_SETUP By default the routine asks the user which ion balance 
;             and abundance files to use via pop-up widgets. If 
;             /no_setup is used then this data is taken from the common 
;             block.
;
;    SUMT     The default is to output the intensity array as an array 
;             of size (nwvl x nT). Setting this keyword performs a sum 
;             over the temperatures to yield a vector of same size as 
;             the input wavelengths, thus producing the complete 
;             free-free spectrum.
;
;    PHOTONS  Gives output emissivity in photon units rather than ergs.
;
;    KEV      If set, then WVL is assumed to contain energies in keV units
;             rather than wavelengths in Angstroms.
;
; CALLS
;
;    SUTHERLAND, ITOH
;
; COMMON BLOCKS
;
;    ELEMENTS
;
; PROGRAMMING NOTES
;
;    The Itoh fitting formula is only valid for (6.0 LE logT LE 8.5). 
;    For temperatures below this, we thus switch to the Sutherland 
;    fitting formula. There is very little (<1%) between the two at 
;    logT=6.
;
;    Itoh also has a constraint on the quantity u=hc/kTl (l=wavelength), 
;    such that (-4 LE log u LE 1.0). The upper limit corresponds to the 
;    continuum being cut-off prematurely at low wavelengths. E.g., for 
;    T=10^6 the cutoff is at 14.39 angstroms. For these low wavelengths 
;    we also use the Sutherland data to complete the continuum. Note that 
;    the continuum at these wavelengths is very weak
;
; MODIFICATION HISTORY
;
;    Ver.1, 5-Dec-2001, Peter Young
;         Completely revised to call the separate itoh.pro and 
;         sutherland.pro routines.
;
;    V. 2, 21-May-2002,  Giulio Del Zanna (GDZ),
;          Corrected the description of the  units.
;          Added verbose keyword and a printout.
;
;    V. 3, 22-May-2002,  Peter Young (PRY)
;          Added MIN_ABUND optional input.
;          Changed ioneq_t to ioneq_logt (GDZ).
;
;    V 4, 25-May-2005, GDZ 
;                  corrected routine header.
;
;    V. 5, 9-Mar-2006, Peter Young
;          Added /kev keyword
;
;    V. 6, 9-May-2008, Peter Young & Kevin Tritz
;          Modifed code so that temperatures can be specified in a unordered
;          list.
;
;    V. 7, 7-Oct-2009, Peter Young
;          Added EM_INT= keyword; simplified how the routine decides
;          which calculation to use (Itoh or Sutherland)
;
; VERSION     :  7, 7-Oct-2009
;
;-

PRO freefree, temp, wvl, int, no_setup=no_setup, sumt=sumt, dem_int=dem_int, $
              photons=photons, min_abund=min_abund, verbose=verbose, $
              kev=kev, em_int=em_int

COMMON elements,abund,abund_ref,ioneq,ioneq_logt,ioneq_ref

IF NOT keyword_set(no_setup) THEN BEGIN
   setup_elements
   no_setup=1
ENDIF

t1=systime(1)

temp=double(temp)
wvl=double(wvl)

;
; * keV units *
; For the calculation of the continuum, it's simply necessary to pass the
; /kev keyword on to itoh.pro and sutherland.pro. However, there's also a
; check in freefree.pro based on wavelength that needs to be corrected. I
; thus define 'chckwvl' for this, with a different definition in the case
; of /kev.
; 
chckwvl=wvl
IF keyword_set(kev) THEN chckwvl=12.398/wvl


IF n_elements(dem_int) NE 0 THEN BEGIN
   IF n_elements(dem_int) NE n_elements(temp) THEN BEGIN
      print,'%FREEFREE: Warning, number of elements of DEM_INT must match'
      print,'  the number of temperatures. Freefree continuum not calculated.'
      return
   ENDIF
ENDIF


;
; I've modified the code below to simply calculate the Itoh
; continuum for the complete range of wavelengths and temperatures,
; and then replace any zero values with the continuum from Sutherland.
; 
; Previously I had split the temperature range into above and below
; 10^6 K, using Sutherland below and Itoh above. I then also had to
; check for wavelengths outside of Itoh's range and replace
; these with Sutherland's data. A potential problem was if
; DEM_INT was specified and there was only one temperature below 10^6
; (e.g., logT=5.9). The routine sutherland.pro would then reject the
; single temperature value as two or more are required for a DEM
; calculation. The new method thus fixes this.
;
; The new method will be slower than the old but since computers are
; now a lot faster this isn't a major problem.
;
itoh,temp,wvl,int,photons=photons, $
     dem_int=dem_int,no_setup=1,min_abund=min_abund,kev=kev,em_int=em_int
;
k=where(int EQ 0.,nk)
IF nk NE 0 THEN BEGIN
  sutherland,temp,wvl,int1,photons=photons, $
             dem_int=dem_int,no_setup=1,min_abund=min_abund,kev=kev,em_int=em_int
  int[k]=int1[k]
ENDIF 

;
; For temperature below 10^6, the Sutherland data is used
;
;; ind=where(alog10(temp) LT 6.0)
;; IF ind[0] NE -1 THEN BEGIN
;;    temp1=temp[ind]
;;    IF n_elements(dem_int) NE 0 THEN dem1=dem_int[ind]
;;    IF n_elements(em_int) NE 0 THEN em1=em_int[ind]
;;    sutherland,temp1,wvl,int1,photons=photons, $
;;      dem_int=dem1,no_setup=1,min_abund=min_abund,kev=kev,em_int=em1
;;    no_setup=1
;; ENDIF


;; ind=where(alog10(temp) GE 6.0)
;; IF ind[0] NE -1 THEN BEGIN
;;    temp2=temp[ind]
;;    IF n_elements(dem_int) NE 0 THEN dem2=dem_int[ind]
;;    IF n_elements(em_int) NE 0 THEN em2=em_int[ind]
;;    itoh,temp2,wvl,int_i,photons=photons, $
;;      dem_int=dem2,no_setup=1,min_abund=min_abund,kev=kev,em_int=em2
;;    no_setup=1
;;   ;
;;   ; For wavelengths which lie outside Itoh's limits, use the Sutherland 
;;   ; data
;;   ;
;;    l_ind=where( alog10(chckwvl) LT alog10(1.439d8/min(temp2))-1.0 )
;;    IF l_ind[0] NE -1 THEN BEGIN
;;       sutherland,temp2,wvl[l_ind],int_s,photons=photons, $
;;         dem_int=dem2,no_setup=1,min_abund=min_abund,kev=kev,em_int=em2
;;       int2=int_i
;;       int2[l_ind,*]=int_s
;;       ind=where(int_i NE 0.)
;;       IF ind[0] NE -1 THEN int2[ind]=int_i[ind]
;;    ENDIF ELSE BEGIN
;;       int2=int_i
;;    ENDELSE
;; ENDIF


;; int = dblarr(n_elements(wvl),n_elements(temp))
;; ind=where(alog10(temp) LT 6.0, n_ind, comp=ind_c, ncomplement=n_indc) 
;; CASE 1 OF 
;;    n_elements(int1) EQ 0: int=int2
;;    n_elements(int2) EQ 0: int=int1
;;    ELSE: begin
;;     IF n_indc NE 0 THEN int[*,ind_c] = int2
;;     IF n_ind NE 0 THEN int[*,ind] = int1
;;   end
;; ENDCASE


;
; PRY, 9-May-2008
;  I've replaced the lines below by the lines above following suggestion
;  by Kevin Tritz
;
; CASE 1 OF 
;    n_elements(int1) EQ 0: int=int2
;    n_elements(int2) EQ 0: int=int1
;    ELSE: int=[[int1],[int2]]
; ENDCASE

siz=size(int)

IF keyword_set(sumt) AND (siz[0] GT 1) THEN int=total(int,2)

IF keyword_set(verbose) THEN BEGIN 
  t2=systime(1)
  print,format='("% FREEFREE: continuum calculated in ",f8.1," seconds")',t2-t1
END  

END
