
;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;		Florence (Italy), the University of Cambridge and the Rutherford Appleton
;		Laboratory (UK).
;       
;
; NAME
;
;     RECOMB_RATE
;
;  PURPOSE:
;
;     This routine computes the ionization rate coefficient of an ion in cm^-3 s^1.
;     The rates are from Dere, K. P., 2007, A&A, submitted
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
;    RECOMBINATION RATE COEFFICIENT (cm^3 s^-1)     
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
;    NONE
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
;
; VERSION     :  1, 17-Nov-2006
;
;-
FUNCTION recomb_rate,gname,temperature,z=z,ion=ion,verbose=verbose
; 
;   
;
if n_params() lt 2 then begin
   print,' '
   print,' > rate = recomb_rate(gname,temperature,[z=z,ion=ion])
   print,'    calculate the recombination rate coefficient in cm^3 s^-1
   print,'    as a function of temperature (K)
   print,' '
   return,-1
endif
;
if not keyword_set(iz) and not keyword_set(ion) then convertname,gname,z,ion
;
str=''
;
t=temperature
;
zion2filename,z,ion,fname

trfile=fname+'.trparams'
if findfile(trfile) ne '' then begin
	;  using the total ionization rates of Nahar et al
	if keyword_set(verbose) then print, 'found trparams file = ',trfile
	openr,lur,trfile,/get_lun
	readf,lur,nt,format='(i5)'
	totalt=fltarr(nt)
	totalr=fltarr(nt)
	for it=0,nt-1 do begin
		readf,lur,tt,tr,format='(2e12.4)'
		totalt(it)=tt
		totalr(it)=tr
	endfor
	free_lun,lur
	recomb=exp(spline(alog(totalt),alog(totalr),alog(t)))
	return,recomb
endif else begin
	;  use both rr and dr files
	rrfile=fname+'.rrparams'
	openr,lur,rrfile,/get_lun
	if keyword_set(verbose) then print, ' using rrparams file = ',rrfile
	rrtype=1
	;  print,' rrtype = ',rrtype
	readf,lur,rrtype
	;
	if rrtype eq 1 then begin
		; a Badnell type
		;
		z=1 
		n=1
		m=1
		w=1
		a=1.
		b=1.
		t0=1.
		t1=1.
		fmt1a='(3i5,e12.4,f10.5,2e12.4)'
		readf,lur,z,ion,m,a,b,t0,t1,format=fmt1a
		free_lun,lur
		recomb=a/(sqrt(t/t0)*(1.+sqrt(t/t0))^(1.-b)*(1.+sqrt(t/t1))^(1.+b))
		;
	endif else if rrtype eq 2 then begin
		; a badnell type
		;
		z=1 
		n=1
		m=1
		w=1
		a=1.
		b=1.
		t0=1.
		t1=1.
		fmt2a='(3i5,e12.4,f10.5,2e11.4,f10.5,e12.4)'
		readf,lur,z,ion,m,a,b,t0,t1,c,t2,format=fmt2a
		b=b+c*exp(-t/t2)
		recomb=a/(sqrt(t/t0)*(1.+sqrt(t/t0))^(1.-b)*(1.+sqrt(t/t1))^(1.+b))
		free_lun,lur
		;
	endif else if rrtype eq 3 then begin
		;  shull type
		z=1
		ion=1
		arad=1.
		xrad=1.
		readf,lur,z,ion,arad,xrad,format='(2i5,2e12.4)'
		recomb=arad/(temperature/1.e+4)^xrad
		free_lun,lur
	endif else begin
		print,' rrtype not defined'
		recomb=0.
	endelse
	;
	
	;  now get dielectronic recombination rate
	;
	drfile=fname+'.drparams'
	;
	z=1
	ion=1
	if findfile(drfile) ne '' then begin
	;
		dr1=fltarr(n_elements(temperature))
		openr,lur,drfile,/get_lun
		if keyword_set(verbose) then print, ' using drparams file = ',drfile
	;
		drtype=1
	;
		readf,lur,drtype
	;
		if drtype eq 1 then begin
			; a Badnell type
		if keyword_set(verbose) then print, ' drtype = ',drtype, ' a Badnell type 1 file'
			ncoef=8
			ener=fltarr(8)
			coef=fltarr(8)
			readf,lur,z,ion,ener,format='(2i5,8e12.4)'
			readf,lur,z,ion,coef,format='(2i5,8e12.4)'
			for icoef=0,ncoef-1 do begin
				dr1(0)=dr1+coef[icoef]*exp(-ener[icoef]/temperature)
			endfor
			dr1=dr1*temperature^(-1.5)
			if keyword_set(verbose)  and keyword_set(dr) then begin
				for icoef=0,ncoef-1 do begin
					print,z,ion,coef[icoef],ener[icoef],format='(2i5,8e12.4)'
				endfor
				print ,' max of dr = ',max(dr1)
			endif
			readf,lur,str
			dr_ref=''
			while not eof(lur) do begin
				readf,lur,str
				dr_ref=dr_ref+str
			endwhile
			free_lun,lur
		endif else if drtype eq 2 then begin
		if keyword_set(verbose) then print, ' drtype = ',drtype, ' a Shull type 1 file'
		;   a Shull type file
			readf,lur,z,ion,adi,bdi,t0,t1,format='(2i5,4e12.4)
			dr1=adi*exp(-t0/temperature)*(1.+bdi*exp(-t1/temperature))/temperature^1.5
			readf,lur,str
			dr_ref=''
			while not eof(lur) do begin
				readf,lur,str
				dr_ref=dr_ref+str
			endwhile
			free_lun,lur
		endif else begin
			if keyword_set(verbose) then print,' drtype not understood'
			dr1=0.
		endelse
	endif else begin
		if (keyword_set(verbose)) then print,' dr file does not exist'
		dr_ref=''
		dr1=fltarr(n_elements(temperature))
	endelse
;
	return,recomb+dr1
;
endelse
;
end
