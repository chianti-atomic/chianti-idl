
;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving
;       George Mason University USA), the University of Michigan (USA),
;       and Cambridge University (UK).
;
;
; NAME
;
;     IONIZ_CROSS
;
;  PURPOSE:
;
;     This routine computes the ionization cross section of an ion in cm^2.
;     The cross sections are from Dere, K. P., 2007, A&A, submitted
;
; INPUTS
;
;    GNAME     A string, specifying the ion name in Chianti style.  For example
;              'o_6' specifies O VII or O+5.
;
;              if the energy is not specified, a set of energies in eV above the
;              ionization potential (IP) is created
;
;
;
; OUTPUTS
;
;    IONIZATION CROSS SECTION
;
; OPTIONAL INPUTS
;
;    ENERGY    Incident electron energy in eV.
;
; KEYWORDS
;
;    Z, ION  Z specified the nuclear charge and ion specifies the ionization state in
;            Chianti style.  Both must be set.  For example, Z=8 and ION=6 specifies the
;            ion 'o_6'
;
; CALLS
;
;    SCALE_BT, DESCALE_BT, READ_SPLUPS, DESCALE_SPLOM, QRP
;
; COMMON BLOCKS
;
;    NONE
;
; PROGRAMMING NOTES
;
;    NONE
;
; MODIFICATION HISTORY
;
;    Ver.1, 17-Nov-2006, Ken Dere
;    Ver.2, 24-Sep-2015, Peter Young
;       Completed pairs of quotes in string expressions to help with
;       viewing code in certain types of text editor. Otherwise no
;       change to code.
;
; VERSION     :  2, 24-Sep-2015
;
;-
function qrp,z,u
   ;
   ;  calculate Qr-prime (equ. 2.12) of Fontes, Sampson and Zhange 1999
   ;
   ;
   aa=1.13d  ; aa stands for A in equ 2.12
   ;

   if z ge 16 then begin
      ; use Fontes Z=20, N=1 parameters
      dd=3.70590d
      c=-0.28394d
      d=1.95270d
      cc=0.20594d
   endif else begin
      ; use Fontes Z=10, N=2 parameters
      dd=3.82652d
      c=-0.80414d
      d=2.32431d
      cc=0.14424d
   endelse
   ;
   if z gt 20 then cc=cc+((z-20.)/50.5)^1.11
   ;
   q=fltarr(n_elements(u))
   ;
   gu=where(u ge 1.,ngu)
   ;
   for igu=0,ngu-1 do begin
      iu=gu(igu)
      q(iu)=( $
         aa*alog(u(iu)) $
         + dd*(1.-1./u(iu))^2 $
         + cc*u(iu)*(1.-1./u(iu))^4 $
         + (c/u(iu) + d/u(iu)^2)*(1.-1./u(iu)) $
         )/u(iu)
   endfor
   ;
   return,q
end
;
FUNCTION ioniz_cross,gname,energy_in,z=z,ion=ion,verbose=verbose
   ;
   ;  calculate ionization cross sections
   ;
   ;
   if n_params() eq 0 then begin
      print,' '
      print,' > cross=ioniz_cross(gname,energy,[z=z,ion=ion]) '
      print,'    calculate the ionization cross section in cm^2 '
      print,'    as a function of energy (eV) '
      print,' '
      return,-1
   endif else if n_params() eq 1 then begin
      read_ip, !xuvtop+'/ip/chianti.ip' ,ip, ipref
      convertname, gname, z, ion
      thisIp = ip[z-1, ion-1]
      thisIpEv = thisIp/8.06554465e+3
      if keyword_set(verbose) then begin
         formatip = '("this IP cm^-1, eV= ", E12.3, E12.3)'
         print, FORMAT=formatip, thisIp, thisIpEv
      endif
      energy_in = thisIpEv*1.01*10.^(0.02*findgen(101))
   endif
   ;
   ;
   ;
   if not keyword_set(iz) and not keyword_set(ion) then convertname,gname,z,ion
   ;
   energy=energy_in
   nenergies=n_elements(energy)
   ;
   iso=z-ion+1
   ;
   if iso eq 1 and z ge 6 then begin
      ; this is the hydrogen sequence and will use Fontes cross sections
      ;
      read_ip,!xuvtop+'/ip/chianti.ip',ip,ref
      ipev=ip/8.06554445d+3

      ;
      ryd=27.2113845d/2.d
      u=energy/ipev(z-1,ion-1)
      ev1ryd=ipev(z-1,ion-1)/ryd
      ;
      a0=0.5291772108d-8
      ;
      a_bohr=!pi*a0^2   ; area of bohr orbit

      if z gt 20 then begin
         ff=(140.+(z/20.)^3.2)/141.
      endif else ff=1.
      ;
      qr=qrp(z,u)*ff
      bb=1.  ; hydrogenic
      qh=bb*a_bohr*qr/ev1ryd^2
      ;
      return, qh
      ;
   endif else if iso eq 2 and z ge 10 then begin
      ;
      ; this is the helium sequence and will use Fontes cross sections
      ;
      read_ip,!xuvtop+'/ip/chianti.ip',ip,ref
      ipev=ip/8.06554445d+3
      ;
      u=energy/ipev(z-1,ion-1)
      ryd=27.2113845d/2.d
      ev1ryd=ipev(z-1,ion-1)/ryd
      ;
      a0=0.5291772108d-8
      a_bohr=!pi*a0^2   ; area of bohr orbit

      if z gt 20 then begin
         ff=(140.+(z/20.)^3.2)/141.
      endif else ff=1.
      ;
      qr=qrp(z,u)*ff
      bb=2.  ; helium sequence
      qh=bb*a_bohr*qr/ev1ryd^2
      ;
      return, qh
      ;
   endif else begin

      zion2filename,z,ion,fname
      ;
      difile=fname+'.diparams'
      ;
      tst=findfile(difile)
      if tst(0) eq '' then begin
         print,' file does not exist',difile
         return,-1
      endif
      ;
      openr,lur,difile,/get_lun
      ;
      idum=1.  & odum=1.
      i5=intarr(5)
      f3=fltarr(3)
      eastr=''
      f1=1.
      f11=1.
      str1=''
      ;
      readf,lur,i5,format='(5i5)'
      kz=i5(0)
      kon=i5(1)
      nspl=i5(2)
      nfac=i5(3)

      if kz ne z and kon ne ion then begin
         print,' not the correct file ',difile
         return,-1
      endif

      ff1=fltarr(nspl+1)
      ff2=fltarr(nspl+1)
      x_spline=fltarr(nspl,nfac)
      y_spline=fltarr(nspl,nfac)
      ev1=fltarr(nfac)
      bt_bethe=fltarr(nfac)
      btf=fltarr(nfac)
      for ifac=0,nfac-1 do begin
         readf,lur,ff1
         btf(ifac)=ff1(0)
         x_spline(0,ifac)=ff1(1:*)
         readf,lur,ff2
         ev1(ifac)=ff2(0)
         ;		print,' ev1 = ',ev1
         y_spline(0,ifac)=ff2(1:*)*1.d-14
         bt_bethe(ifac)=y_spline(nspl-1)
      endfor
      neaev=i5(4)
      if neaev gt 0 then begin
         readf,lur,eastr  ; lur,str1   ;   ,format='(e12.3)'
         eastra=strsplit(eastr,' ',/extract)
         f1=float(eastra)
         ;
         ;        read in EA collision strengths -----------
         ;
         eaname=fname+'.easplom'
         read_splups,eaname,ea_splom,ea_splomref
      endif
      ;
      free_lun,lur
      ;
      energy_in=energy
      ;
      goode=where(energy ge ev1(0),nev)
      if nev eq 0 then begin
         print,' no energies above ionization potential = ',ev1(0)
         return,-1
      endif
      energy=energy(goode)
      ;
      dum=1.
      scale_bt,energy,dum,dum,btf(0),ev1(0),btev,dum,dum
      ;
      y2=spl_init(x_spline(*,0),y_spline(*,0))
      btcross=spl_interp(x_spline(*,0),y_spline(*,0),y2,btev)
      descale_bt,btev,btcross,dum,btf(0),ev1(0),evd,direct_cross,odum
      ;
      direct_cross=direct_cross
      ;
      nenergy=n_elements(energy)
      ;
      for ifac=1,nfac-1 do begin
         good=where(energy ge ev1(ifac),ngood)
         if ngood gt 0 then begin
            scale_bt,energy(good),dum,dum,btf(ifac),ev1(ifac),btev,dum,dum
            y2=spl_init(x_spline(*,ifac),y_spline(*,ifac))
            btcross=spl_interp(x_spline(*,ifac),y_spline(*,ifac),y2,btev)
            descale_bt,btev,btcross,dum,btf(ifac),ev1(ifac),evd2,dcross,odum
            newcross=[fltarr(nenergy-ngood),dcross]
            direct_cross=direct_cross+newcross
         endif
      endfor
      ;
      ea_cross=fltarr(nev)
      ;
      if neaev gt 0 then begin
         descale_splom,ea_splom,energy,ea_om
         nsplom=n_elements(ea_splom.lvl1)
         if n_elements(f1) lt nsplom then f1=replicate(f1,nsplom)
         for isplom=0,nsplom-1 do begin
            ea_cross(0)=ea_cross(*)+f1(isplom)*8.797e-17*ea_om(*,isplom)/(energy/13.6056981d)
         endfor
         cross=direct_cross + ea_cross
      endif else cross=direct_cross
      ;
      if nenergies - nev eq 0 then begin
         return, cross
      endif else begin
         return,[fltarr(nenergies-nev),cross]
      endelse
      ;
   endelse
   ;
end
