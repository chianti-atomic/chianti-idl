;+
; NAME:
;	TWO_PHOTON
;
; PURPOSE:
;	Calculate the 2 photon continuum from a hot, low density plasma.
;
;       For the hydrogen isoelectronic sequence, A values
;             Parpia, F. A., and Johnson, W. R., 1982, Phys. Rev. A, 26, 1142.
;       and distribution function as a function of Z is taken from: 
;             Goldman, S.P. and Drake, G.W.F., 1981, Phys Rev A, 24, 183
;
;       For the helium isoelectronic sequence, A values are from:
;             Drake, G.W.F., 1986, Phys. Rev. A, 34, 2871
;       and the distribution function as a function of Z is taken from:
;             Drake, G.W.F., Victor, G.A., Dalgarno, A., 1969, Phys. 
;             Rev. A, 180, 25.
;       in this paper the distribution is only given up to Z=10 but  
;       extrapolation to higher Z appears to be OK.
;
;       Note that, unlike the freefree and freebound routines, two_photon 
;       requies the electron density to be specified. This is because there 
;       is a call to pop_solver.
;
;       To calculate the continuum spectrum the user should specify
;       whether the plasma structure is specified through a differential
;       emission measure (DEM) curve or through emission measure (EM)
;       values. For the former, the DEM values are specified through the
;       DEM_INT keyword and for the latter the EM values are specified
;       through the EM_INT keyword. If neither are specified EM_INT
;       values of unity are assumed.
;
;       If DEM_INT values are specifed they must be specified on a
;       temperature scale that has constant intervals in log10(T). EM_INT
;       values do not have this restriction.
;
;       Note that EM_INT must be used for isothermal plasmas.
;
; CALLING SEQUENCE:
;       TWO_PHOTON, temperature, wavelength, intensity
;
; INPUTS:
;       Temp:   Temperature (in K). Can be an array. If the input
;               DEM_INT= is used then the temperatures must be spaced at
;               equal intervals of 
;               log10(T). E.g., log T=6.0, 6.1, 6.2, 6.3, ...
;       Wavelength:  Wavelengths in Angstroms. If /keV is set, then WVL is
;                    assumed to contain energies in keV units.
;
;
; OPTIONAL INPUTS:
;       Edensity:     Electron density in cm^-3, can also be a 1D array 
;                     of the same size as Temperature. If there are several 
;                     temperatures specified but only one density, then 
;                     density is assumed the same at all temperatures.
;                     If not specified, then default densities of 10^10 
;                     electrons/cm^3 are assumed.
;
;       Dem_Int:  This should be the same size as TEMP and contain
;                 differential emission measure values. Specifying
;                 DEM_INT places a restriction on TEMP (see above). Note
;                 that DEM_INT values are multiplied by the factor dT
;                 within the routine, whereas EM_INT values are not.
;
;       Em_Int:   This should be the same size as TEMP and contain
;                 emission measure values. Note
;                 that DEM_INT values are multiplied by the factor dT
;                 within the routine, whereas EM_INT values are not. If
;                 neither EM_INT or DEM_INT are specified, then a EM_INT
;                 vector of same size as TEMP with values of unity is
;                 assumed.
;
;       Abund_File:  The name of an element abundance file. If not
;                 set, then !abund_file will be used.
;       Ioneq_File:  The name of an ionization balance file. If not
;                 set, then !ioneq_file is used.
;       Sngl_Ion: A string array specifying one or more ions for which
;                 the calculation is required.
;       Element:  Either an integer or string array specifying
;                 elements to include in the calculation.
;
; KEYWORD PARAMETERS:
;
;       MIN_ABUND:  If set, calculates the continuum only from those 
;                   elements which have an abundance greater or equal to min_abund.  
;                   Can speed up the 
;                   calculations.  For example:
;                   abundance (H)  = 1.
;                   abundance (He) = 0.085
;                   abundance (C)  = 3.3e-4
;                   abundance (Si) = 3.3e-5
;                   abundance (Fe) = 3.9e-5
;
;       SUMT        If several temperatures have been specified, then /sumt 
;                   will sum the emissivities over the different 
;                   temperatures, giving an output INTENSITY that has the 
;                   same size as WAVELENGTH.
;
;       PHOTONS     If set the continuum emissivity is output in photon 
;                   units rather than erg units.
;
;       VERBOSE  If set, information is printed to screen as the routine runs.
; 
;       KEV      If set, then WVL is assumed to contain energies in keV units
;                rather than wavelengths in Angstroms.
;
;	NO_SETUP:   This keyword is obsolete now. Retained for
;	            backwards compatibility.
;
; OUTPUTS:
;
;	RAD   2 photon continuum intensity in units of
;             10^-40 erg cm^3/s/sr/Angstrom  per unit emission measure 
;             ( integral(N_H N_e dh) in cm^-5) if a DEM is not defined. 
;
;             If DEM values are defined, it is assumed that they are given
;             as N_H N_e dh/dT.  The units are 
;             10^-40 erg/cm^2/s/sr/Angstrom 
;
;             If T is given as a 1-D array, then the output will be a 2-D
;             array, with one element for each temperature and wavelength 
;             (but also see SUMT).
;
;             If the keyword /keV is set, then the units of INT will be 
;             10^-40 erg cm^3/s/sr/keV
;
; OPTIONAL OUTPUTS:
;       H_Rad:  The same as RAD, but including only the emission from
;               the H-like ions.
;       He_Rad: The same as RAD, but including only the emission from
;               the He-like ions.
;
; CALLS
;       POP_SOLVER, CH_SETUP_ION, READ_MASTERLIST, READ_ABUND, 
;       CONVERTNAME, READ_IONEQ
;
; EXAMPLE:
;        IDL> two_photon,1.e+6,,wvl,int
;        IDL> two_photon,1.e+6,wvl,int,min_abund=3.e-5
;        IDL> two_photon,1.e+6,wvl,int,edens=3e10
;        IDL> two_photon,1.e+6,wvl,int,sngl_ion=['o_7','o_8']
;        IDL> two_photon,1.e+6,wvl,int,element=8
;        IDL> two_photon,1.e+6,wvl,int,element='O'
;        IDL> two_photon,1.e+6,wvl,int,element=['C','N','O']
;
; PROGRAMMING NOTES
;
;    For He 2-photon decays, the distribution function is from Table II 
;    of Drake et al. (1969), except that the values have been divided by 
;    the A-value from Drake (1986).
;
;
; MODIFICATION HISTORY:
; 	Written by:	Ken Dere
;	February 2001:  Version 1.0
;
;       Ver.2, 19-Dec-2001, Peter Young
;           Now allows an array of temperatures.
;           Added /SUMT keyword.
;           Added DEM_INT= optional input.
;           Switched to using spl_init & spl_interp for the spline fits.
;           Corrected a small bug.
;
;       Ver.3, 20-Dec-2001, Peter Young
;           Corrected a bug related to density indices.
;
;       Ver.4, 23-Apr-2002, Peter Young
;           Added /photons keyword.
;
;       Ver.5, 28-May-2002, Peter Young
;           Corrected way in which DEM_INT is handled.
;
;       V. 6, 28-May-2002, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;             Corrected the description of the  units and various
;             inaccuracies in the header.
;    
;       V.7, 14-Aug-02, GDZ
;             Added keyword VERBOSE
;
;       V.8,  4-Oct-2003, GDZ
;                  modified the input to POP_SOLVER, now it is dealt with an
;                  input structure.
;
;       V.9,  8-Jun-2004, EL
;                  modified the input to POP_SOLVER, now it includes ion/rec
;
;       V.10, 5-Jul-2005
;                  corrected problems with the input structure for pop_solver
;
;       v.11 29-Jul-2005, GDZ
;                  fixed bug, only define ioneq_ionrec when files are
;                  there (it was failing with neutrals)
;
;       v.12,  31-Aug-2005, GDZ
;                  missed one correct definition of ioneq_ionrec. Fixed it now.
;
;       v.13, 10-Mar-2006, Peter Young
;                  added /keV keyword 
;
;       v.14, 12-Jun-2009, Enrico Landi
;                  Changed the definition of the temperature array for ion fractions
;                  in the IONREC variable, now taken directly from the output of
;                  READ_IONEQ.PRO
;
;       v.15, 11-Sep-2009, Peter Young
;            Routine now checks that temperatures are evenly spaced
;            when DEM_INT is specified. Also 0.1 dex spacing is no
;            longer assumed.
;
;       v.16, 6-Oct-2009, Peter Young
;            Introduced EM_INT keyword for specifying emission measure
;            values for isothermal plasmas.
;
;       v.17, 13-Apr-2010, Enrico Landi
;            Remove a bug in the definition of temperature for ion/rec effects when
;            dealing with He-like ions
;
;       v.18, 6-Oct-2010, Peter Young
;            Changed the 'list' array to 'mlist' for compatibility
;            with IDL 8.
;
;       v.19, 5-Sep-2012, Peter Young
;            Updated header, and changed check on min_abund to make it
;            consistent with freefree.pro and freebound.pro (check
;            used to be gt min_abund, now it's ge min_abund).
;
;       v.20, 21-Aug-2015, Ken Dere
;               the index for the upper level has changed for several of the
;               helium-like sequence
;
;       v.21, 19-Jan-2018, Peter Young
;            This is a major update. The common blocks have been
;            removed, and I now call ch_setup_ion to load the atomic
;            data for the ions. In the process of doing this I found
;            three bugs in the original code (i) proton rates were not
;            loaded previously, yet they were important for c_6; (ii)
;            the upper levels for the He-like ions are different for
;            most of the dielectronic ions, so these are now defined
;            separately (heseq_lvl2d); (iii) for He-like dielectronic
;            ions, the routine was incorrectly using atomic data from
;            the previous ion in masterlist. I've also added a
;            check on the upper levels of the He-like ions to make
;            sure they really do correspond to the two-photon
;            transition.
;
;       v.22, 12-Aug-2019, Peter Young
;            The abund_file and ioneq_file inputs were not working, so
;            this has been fixed.
;
;       v.23, 10-Mar-2020, Peter Young
;            Added ELEMENT= optional input; changed how SNGL_ION input
;            was interpreted; set /quiet keyword for calls to
;            ch_setup_ion; no longer uses default abundance file if
;            abund_file not specified.
;
;       v.24, 12-Nov-2020, Peter Young
;            The upper levels of the two-photon transitions are no
;            longer hard-coded. Instead they are accessed through the
;            two_photon tag of the ch_setup_ion output.
;-
pro two_photon,temperature,wvl,rad, no_setup=no_setup, $
               min_abund=min_abund, masterlist=masterlist, $
               sngl_ion=sngl_ion, sumt=sumt, dem_int=dem_int, $
               photons=photons, edensity=edensity, verbose=verbose, $
               kev=kev, em_int=em_int, $
               abund_file=abund_file, ioneq_file=ioneq_file, $
               h_rad=h_rad, he_rad=he_rad, $
               element=element
               


if n_params(0) lt 3 then begin
   print,'   IDL> two_photon,temperature,wavelength,intensity '
   print,'              [ min_abund= , dem_int= , /sumt '
   print,'               ,masterlist= , sngl_ion= , edensity=, '
   print,'               /verbose, /kev, /photons, em_int=, $'
   print,'               ioneq_file=, abund_file=, h_rad=, he_rad= $'
   print,'               element= ]'
   return
endif

t1=systime(1)

kb=1.38062d-16                  ;  erg deg-1
h=6.6262d-27
c=2.997925d+10
ryd=2.17992d-11                 ; erg
factor=h*c/(4.*!pi)
rescale=1.d+40
;
; PRY, 18-Jan-2018
; Here the level index of upper level of the two-photon transition of
; the helium sequence ions is hard-coded. Note that -1 indicates the
; ion does not have a model.
; Note that the index can be different for the dielectronic ions.
; PRY, 12-Nov-2020
; These are no longer used.
;
heseq_lvl2 = [-1, 3,-1,-1,-1,5,6,6,-1,6, 6,6,5,5, 3,5, 3,5, 3,5,-1,-1,-1,-1,-1,4,-1,4,-1, 4]
heseq_lvl2d =[-1,-1,-1,-1,-1,3,3,3,-1,3,-1,3,3,3,-1,3,-1,3,-1,3,-1,-1,-1,-1,-1,3,-1,3,-1,-1]
;
temperature=double(temperature)
wvl=double(wvl)
nwvl=n_elements(wvl)
nt=n_elements(temperature)
;
rad=dblarr(nwvl,nt)
h_rad=rad
he_rad=rad
distr=dblarr(nwvl)

;
; * keV units *
;
IF keyword_set(kev) THEN BEGIN
  wvl_save=wvl
  efact=12.398/wvl^2   ; hc/E^2 factor to multiply continuum by
  efact=reverse(efact)
  wvl=12.398/wvl   ; convert energies to angstroms
  wvl=reverse(wvl)
ENDIF ELSE BEGIN
  efact=1d0
ENDELSE


IF keyword_set(no_setup) THEN BEGIN
  print,'% TWO_PHOTON: the keyword /no_setup is obsolete now. Please check the routine header.'
ENDIF 

;
; Read the abundance file
;
IF n_elements(abund_file) NE 0 THEN BEGIN
  chck=file_search(abund_file,count=count)
  IF count EQ 0 THEN BEGIN
    print,'% FREEBOUND: the specified abund_file was not found. Returning...'
    return
  ENDIF 
ENDIF ELSE BEGIN
  abund_file=''
ENDELSE 
read_abund,abund_file,abund,ab_ref, element=element



;
; Read the ioneq file
;
IF n_elements(ioneq_file) EQ 0 THEN BEGIN
  ioneq_file=!ioneq_file
ENDIF
read_ioneq,ioneq_file,ioneq_t,ioneq,ioneq_ref,element=element,sngl_ion=sngl_ion


IF n_elements(edensity) EQ 0 THEN BEGIN
  edensity=1d10
  print,'% TWO_PHOTON: EDENSITY not specified, assuming 10^10.'
ENDIF

IF n_elements(edensity) EQ 1 THEN edens=dblarr(nt)+edensity $
ELSE edens=edensity

;
; DEM_INT and EM_INT implementation
; ---------------------------------
;  If DEM_INT is specified then it's necessary to multiply by the
;  DEM values by the factor dT (=ln(10)*T*d(log10(T)).
;
;  I implement this by introducing a keyword dT. If dT=1 then a DEM is
;  assumed, while dt=0 implies EM values. In the latter case DEM_INT
;  is set to be EM_INT.
;
;  EM_VALS is set to be DEM_INT or EM_INT and is used in the rest of
;  the code.
;
CASE 1 OF 
  n_elements(dem_int) NE 0 AND n_elements(em_int) NE 0: BEGIN
    print,'%TWO_PHOTON: Please specify only one of DEM_INT and EM_INT. Returning...'
    rad=0.
    return
  END
 ;
  n_elements(dem_int) EQ 0 AND n_elements(em_int) NE 0: BEGIN
    em_vals=em_int
    dt=0
  END 
 ;
  n_elements(dem_int) NE 0 AND n_elements(em_int) EQ 0: BEGIN
    em_vals=dem_int
    dt=1
  END 
 ;
  n_elements(dem_int) EQ 0 AND n_elements(em_int) EQ 0: BEGIN
    em_vals=dblarr(nt)+1.
    dt=0
  END 
ENDCASE

;
; Make sure EM_VALS has same size as TEMP
;
em_vals=double(em_vals)
IF n_elements(em_vals) NE nt THEN BEGIN
  print,'%TWO_PHOTON: Warning, number of elements of DEM_INT (EM_INT) must match the'
  print,'  number of temperatures. Returning...'
  rad=0.
  return
ENDIF 


;
; If /DT set then need to multiply by dT for each element of EM_VALS.
;
IF keyword_set(dt) THEN BEGIN
  IF nt LT 2 THEN BEGIN
    print,'%TWO_PHOTON: Two or more temperatures must be specified if DEM_INT is specified.'
    print,'             Returning...'
    rad=0.
    return
  ENDIF 
 ;
  lt1=alog10(temperature[0:nt-2])
  lt2=alog10(temperature[1:*])
  mean_lt=mean(lt2-lt1)
  max_LT=max(lt2-lt1)
  min_LT=min(lt2-lt1)
  IF max_LT-mean_LT GE 0.05*mean_LT OR mean_LT-min_LT GE 0.05*mean_LT THEN BEGIN
    print,'%TWO_PHOTON: Error, LOG10 TEMP needs to be evenly spaced if DEM_INT is specified.'
    print,'             E.g., LOG10 TEMP should be spaced at 0.1 dex intervals.'
    rad=0.
    return
  ENDIF 
  demfactor=temperature*mean_lt/alog10(exp(1))
END ELSE BEGIN
  demfactor=dblarr(nt)+1.
ENDELSE 
;
; Create dem_arr
;
dem_arr= em_vals*demfactor



;  read master list of ions
;
if datatype(masterlist,2) gt 0 then begin
   mname=masterlist
endif else begin
   mname=concat_dir(concat_dir(!xuvtop,'masterlist'), 'masterlist.ions')
endelse
;
read_masterlist,mname,mlist
nlist=n_elements(mlist)
;
n_ioneq_t=n_elements(ioneq_t)


;
;   read in hydrogen sequence two photon data  *************
;
datafile=concat_dir(concat_dir(!xuvtop,'continuum'), 'hseq_2photon.dat')
;
openr,luw,datafile,/get_lun
;
y0=fltarr(17)
readf,luw,y0
;
z0=fltarr(4)
readf,luw,z0
;
nz=30
avalue=fltarr(nz)
asum=fltarr(nz)
psi0=fltarr(17,nz)              ;  psi is the 2 photon distribution function
;                      psi for H- and He-like are normalized differently
f1=fltarr(19)
i1=1
;
for iz=0,nz-1 do begin
   readf,luw,i1,f1
   avalue(iz)=f1(0)
   asum(iz)=f1(1)
   psi0(0,iz)=f1(2:*)
endfor
;
free_lun,luw
;
;   psi is properly normalized to give an integral of 2.0
;
;   finished reading in 2 photon data for hydrogen
;
;
;  ********  first go through for hydrogen   *******************************
;
for ilist=0,nlist-1 do begin
;
   gname=mlist(ilist)
   convertname,gname,iz,ion
;
;
   hydrogentest = (iz-ion) eq 0
;
;
   if hydrogentest then begin
;
;
      this_abund=abund(iz-1)
;
      if keyword_set(min_abund) then begin
         abundtest = this_abund GE min_abund
      endif else abundtest = this_abund gt 0.
;
      if abundtest then begin
                                ;
         this_ioneq=ioneq(*,iz-1,ion-1)
         
         ind_i=where(this_ioneq gt 0.,ngt)
         IF ngt GE 2 THEN BEGIN
            goodt=ioneq_t[ind_i]
            goodi=this_ioneq[ind_i]
                                ;
            ltemp=alog10(temperature)
            t_ind=where( (ltemp GE min(goodt)) AND (ltemp LE max(goodt)) )
            IF t_ind[0] NE -1 THEN BEGIN
               y2=spl_init(goodt,alog10(goodi))
               ioneq1=spl_interp(goodt,alog10(goodi),y2,ltemp[t_ind])
               ioneq1=10.^ioneq1
               temps=temperature[t_ind]
            ENDIF ELSE BEGIN
               ioneq1=0.
            ENDELSE
         ENDIF ELSE BEGIN
            ioneq1=0.
         ENDELSE  
                                ;
         ioneqtest=ioneq1[0] gt 0.
                                ;
                                ;
         if ioneqtest then BEGIN
                                ;
           input=ch_setup_ion(gname,ioneq_file=ioneq_file,abund_file=abund_file,/quiet)
           ecm=input.ecm
           two_photon=input.two_photon
           
           wvl0=1.d+8/(ecm(1)-ecm(0))
            w_ind = where(wvl gt wvl0,wvltest)
                                ;
            if wvltest gt 0 AND n_tags(two_photon) NE 0 then begin
               pop_idx=two_photon.lvl-1

               y=wvl0/wvl[w_ind]

               y2=spl_init(y0,psi0[*,iz-1])
               distr=y*spl_interp(y0,psi0[*,iz-1],y2,y)/asum[iz-1]/wvl[w_ind]

               nt=n_elements(temps)
               FOR i=0,nt-1 DO BEGIN
                 pop_solver,input, temps[i],edens[t_ind[i]],pop
                  IF keyword_set(photons) THEN BEGIN
                     distr1=rescale/4d0/!pi*avalue[iz-1]*this_abund* $
                       distr * $
                       (ioneq1[i]*pop[pop_idx]/edens[t_ind[i]]) * dem_arr[t_ind[i]]
                  ENDIF ELSE BEGIN
                     distr1=rescale*factor*1d8*avalue[iz-1]*this_abund* $
                       (distr/wvl[w_ind]) * $
                       (ioneq1[i]*pop[pop_idx]/edens[t_ind[i]]) * dem_arr[t_ind[i]]
                  ENDELSE
                  h_rad[w_ind,t_ind[i]]=h_rad[w_ind,t_ind[i]]+distr1
               ENDFOR
                                ;
            ENDIF               ;  wvl test
;
         endif                  ;  ioneq test
;
      endif                     ;  abundance test
;
   endif                        ;  hydrogentest
;
endfor                          ; ilist

t2=systime(1)

;
;  now go through for the  helium iso-electronic sequence
;
;
datafile=concat_dir(concat_dir(!xuvtop,'continuum'), 'heseq_2photon.dat')
;
openr,luw,datafile,/get_lun
;
y0=fltarr(41)
readf,luw,y0
;
nz=30
avalue=fltarr(nz)
psi0=fltarr(41,nz)
;
f1=fltarr(42)
i1=1
;
FOR iz=1,nz-1 DO BEGIN
   readf,luw,i1,f1
   avalue(iz)=f1(0)
   psi0(0,iz)=f1(1:*)
ENDFOR 
;
free_lun,luw



for ilist=0,nlist-1 do BEGIN

   gname=mlist(ilist)  ; This line was missing prior to 18-Jan-2018.
   convertname,gname,iz,ion,dielectronic=dielectronic

   heliumtest=(iz-ion EQ 1)

   if heliumtest then BEGIN
;
      this_abund=abund(iz-1)
;
      if keyword_set(min_abund) then begin
         abundtest = this_abund GE min_abund
      endif else abundtest = this_abund gt 0.
;
      if abundtest then begin
;
         this_ioneq=ioneq(*,iz-1,ion-1+dielectronic)
         
         ind_i=where(this_ioneq gt 0.,ngt)
         IF ngt GE 2 THEN BEGIN
            goodt=ioneq_t[ind_i]
            goodi=this_ioneq[ind_i]
                                ;
            ltemp=alog10(temperature)
            t_ind=where( (ltemp GE min(goodt)) AND (ltemp LE max(goodt)) )
            IF t_ind[0] NE -1 THEN BEGIN
               y2=spl_init(goodt,alog10(goodi))
               ioneq1=spl_interp(goodt,alog10(goodi),y2,ltemp[t_ind])
               ioneq1=10.^ioneq1
               temps=temperature[t_ind]
            ENDIF ELSE BEGIN
               ioneq1=0.
            ENDELSE
         ENDIF ELSE BEGIN
            ioneq1=0.
         ENDELSE  
                                ;
         ioneqtest=ioneq1[0] gt 0.
                                ;
                                ;
         if ioneqtest then BEGIN
            input=ch_setup_ion(gname,ioneq_file=ioneq_file,abund_file=abund_file,/quiet)
            ecm=input.ecm
            two_photon=input.two_photon

            wvl0=1.d+8/(ecm(2)-ecm(0))
            w_ind = where(wvl gt wvl0,wvltest)

            if wvltest gt 0 AND n_tags(two_photon) NE 0 then begin

               pop_idx=two_photon.lvl-1
               
               y=wvl0/wvl[w_ind]

               y2=spl_init(y0,psi0[*,iz-1])
               distr=y*spl_interp(y0,psi0[*,iz-1],y2,y)/wvl[w_ind]

               nt=n_elements(temps)
               FOR i=0,nt-1 DO BEGIN
                  pop_solver, input, temps[i],edens[t_ind[i]],pop
                  IF keyword_set(photons) THEN BEGIN
                     distr1=rescale/4d0/!pi*avalue[iz-1]*this_abund* $
                       distr * $
                       (ioneq1[i]*pop[pop_idx]/edens[t_ind[i]]) * dem_arr[t_ind[i]]
                  ENDIF ELSE BEGIN
                     distr1=rescale*factor*1d8*avalue[iz-1]*this_abund* $
                       (distr/wvl[w_ind]) * $
                       (ioneq1[i]*pop[pop_idx]/edens[t_ind[i]]) * dem_arr[t_ind[i]]
                  ENDELSE
                  he_rad[w_ind,t_ind[i]]=he_rad[w_ind,t_ind[i]]+distr1
               ENDFOR
                                ;
            ENDIF               ;  wvl test
         ENDIF
;
      endif                     ;  abundance test
;
   endif                        ;  helium sequence test
;
endfor                          ;  ilist


rad=h_rad+he_rad

;
; multiply int by the hc/E^2 factor if /kev set
;
IF keyword_set(kev) THEN BEGIN
  nt=n_elements(temperature)
  e_arr=efact#(dblarr(nt)+1d0)
  rad=rad*e_arr
 ;
 ; reverse the wavelength dimension to match energy ordering rather than
 ; wavelength ordering
  rad=rad[reverse(indgen(nwvl)),*]
 ;
 ; and set wavelengths back to energy units
  wvl=wvl_save
ENDIF


IF keyword_set(sumt) AND (n_elements(temperature) GT 1) THEN rad=total(rad,2)

t3=systime(1)

IF keyword_set(verbose) THEN $
  print,format='("% TWO_PHOTON: continuum calculated in ",f8.1," seconds")',t3-t1

end

