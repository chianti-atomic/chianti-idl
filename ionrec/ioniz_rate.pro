
;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;               Florence (Italy), the University of Cambridge and the Rutherford Appleton
;               Laboratory (UK).
;
;
; NAME
;
;     IONIZ_RATE
;
;  PURPOSE:
;
;     This routine computes the ionization rate coefficient of an ion in cm^-3 s^1.
;     The rates are from Dere, K. P., 2007, A&A, 466, 771
;
; INPUTS
;
;    GNAME     A string, specifying the ion name in Chianti style.  For example
;              'o_6' specifies O VII or O+5.
;
;    TEMPERATURE    Temperature (K).
;
;
; OUTPUTS
;
;    IONIZATION RATE COEFFICIENT (cm^3 s^-1)
;
; OPTIONAL INPUTS
;
;    NONE
;
; KEYWORDS
;
;    Z, ION  Z specified the nuclear charge and ion specifies the ionization state in
;            Chianti style.  Both must be set.  For example, Z=8 and ION=6 specifies the
;            ion 'o_6'
;
; CALLS
;
;    SCALE_BT, DESCALE_BT, READ_SPLUPS, DESCALE_ALL
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
FUNCTION ioniz_rate,gname,temperature_in,z=z,ion=ion,verbose=verbose
;
;  calculate ionization cross sections 
;   
;
if n_params() lt 2 then begin
   print,' '
   print,' > cross=ioniz_rate(gname,temperature,[z=z,ion=ion]) '
   print,'    calculate the ionization rate coefficient in cm^3 s^-1 '
   print,'    as a function of temperature (K) '
   print,' '
   return,-1
endif
;
;
alpha=5.287d+13   ; not the fine structure  constant
kb=1.380626d-16
kbev=8.617343d-5
;
tev=kbev*temperature_in
temperature=temperature_in
ntemps=n_elements(temperature)
;
dirate=dblarr(ntemps)
earate=dblarr(ntemps)
totrate=dblarr(ntemps)
;
;   gauss laguerre coefficients for n=12
;
ngl=12
xgl=[0.115722117358021d,0.611757484515131d,1.512610269776419d,2.833751337743509d     $
    ,4.599227639418353d,6.844525453115181d,9.621316842456871d,13.006054993306350d    $
    ,17.116855187462260d,22.151090379396983d,28.487967250983992d,37.099121044466926]


wgl=[2.647313710554435d-01,3.777592758731382d-01,2.440820113198774d-01,9.044922221168074d-02  $
    ,2.010238115463406d-02,2.663973541865321d-03,2.032315926629993d-04,8.365055856819753d-06  $          
    ,1.668493876540914d-07,1.342391030515027d-09,3.061601635035012d-12,8.148077467426124d-16]
;
;
if not keyword_set(iz) and not keyword_set(ion) then convertname,gname,z,ion
;
;
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
;
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
		eaname=fname+'.easplups'
		read_splups,eaname,ea_splups,ea_splupsref
endif
;
free_lun,lur
;
;
;
;
dum=1.
;
for it=0,ntemps-1 do begin
   ;
	for ifac=0,nfac-1 do begin
		x0=double(ev1(ifac)/tev(it))
		beta=sqrt(kb*temperature(it))
		egl=ev1(ifac)+xgl*tev(it)
		scale_bt,egl,dum,dum,btf(ifac),ev1(ifac),btgl,dum,dum
		;  di splines were made with a cubic spline fit
		y2=spl_init(x_spline(*,ifac),y_spline(*,ifac))
		btcross=spl_interp(x_spline(*,ifac),y_spline(*,ifac),y2,btgl)
		btcross=double(btcross)
		descale_bt,btgl,btcross,dum,btf(ifac),ev1(ifac),evd2,crossgl,odum
		newcross=alpha*beta*exp(-x0)*(total(wgl*xgl*crossgl)+x0*total(wgl*crossgl))
		dirate(it)=dirate(it)+newcross
	endfor
endfor  ;  it
;
;
;
;
if neaev gt 0 then begin
	nups=n_elements(ea_splups.lvl1)
	if n_elements(f1) lt nups then f1=replicate(f1,nups)
	for iups=0,nups-1 do begin
		descale_all,temperature,ea_splups,iups,ups
		xt=tev/(ea_splups(iups).de*13.6056981)
		earate=earate+f1(iups)*8.63d-6*ups*exp(-1./xt)/(sqrt(temperature))
	endfor
	totrate=dirate + earate
endif else totrate=dirate
;
;
return,totrate
;
end
