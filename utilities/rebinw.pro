function rebinw,f,x,y,perbin=perbin,nbin=nbin,wrange=wrange, _extra=e
;+
;function	rebinw
;	rebins the input array on a new grid and returns the output array.
;	this routine is preferred over simple INTERPOL or SPLINE because
;	this is guaranteed to e.g., conserve energy while rebinning spectra.
;
;syntax
;	ff=rebinw(f,x,y,/perbin,nbin=nbin,wrange=wrange)
;
;parameters
;	f	[INPUT; required] array values
;	x	[INPUT; required] absissae for F
;		1: if size matches that of F, assumed to be mid-bin values
;		2: if size is size(F)+1, assumed to be bin boundaries
;		   (i.e., all the bin-beginners and a final bin-ender)
;		3: if neither of the above, this is assumed to be the
;		   desired *output* grid, and the absissae of the input
;		   are assumed to span the linear range 1..N(F)
;		   * Y is ignored on input, but will be overwritten on output
;	y	[I/O] the new absissa values
;		1: if scalar or 1-element vector on input, assumed to be
;		   number of bins in output
;		   * linear binning if +ve, log binning if -ve
;		2: if >1-element vector on input, assumed to be the desired
;		   output grid of bin boundaries, with the last element
;		   defining the final bin-ending value.
;		3: if not defined, takes value from NBIN
;		NOTE: F extrapolates as zeros.
;
;keywords
;	perbin	[INPUT] if set, assumes that units of F(X) are [.../bin]
;	nbin	[I/O] number of bins in output
;		* overwritten if defined via X or Y
;		* linear binning if +ve, log binning if -ve
;		* default: 1
;	wrange	[INPUT] output grid range
;		* ignored if defined via array Y (or X -- case 3)
;		* overrides if determined via X and CASE 1 of Y
;		* default: [1.,500.]
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	beware that if Y is not sorted in ascending order, on output
;	it will be!
;	requires external procedure FINDEX.PRO
;
;history
;	vinay kashyap (Jan98)
;	forced x to be required inputs (Jun98)
;	changed call to INTERPOL to combo call to INTERPOLATE/FINDEX (VK;99Jul)
;	corrected long-standing "feature" of NBIN being useless, changed
;	  default behavior of Y=UNDEFINED, NBIN, and X=/=F (VK; MarMM)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(f) & nx=n_elements(x)
ny=n_elements(y)
if np lt 2 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='Input function undefined' else $
  if nx eq 0 then ok='bin boundaries undefined' else $
   if nf eq 1 then ok='Cannot make a spectrum out of one element'
if ok ne 'ok' then begin
  print,'Usage: g=rebinw(f,x,y,/perbin,nbin=nbin,wrange=wrange)'
  print,'  rebins F(X) into G(Y) on a different grid while conserving norm'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	initialize via keywords
if n_elements(nbin) eq 0 then nbin=1L else nbin=long(nbin(0))
if nbin eq 0 then nbin=1L	;NBIN cannot be 0
if n_elements(wrange) ge 2 then wr=[wrange(0),wrange(1)] else wr=[1.,500.]

;	compatibility checks
xx=x & wrr=[min(xx),max(xx)]		;case 2: default
if nx eq nf then begin			;(case 1: mid-bin values
  dx=x(1:*)-x & dx=0.5*dx & xx=[x-dx,x(nx-1)+dx(nx-2)*[-1,1]]
  wrr=[xx(0),xx(nx)]
endif else if nx ne (nf+1) then begin 	;)(case 3: X is the output grid
             xx=lindgen(nf+1)+1
	     if nx gt 0 then yy=x else yy=1L
	     wrr=[xx(0),xx(nf)]
	     y=yy & ny=n_elements(yy)		;ignore Y on input
           endif			;handled X)
nx=n_elements(xx)
if n_elements(wrange) eq 2 then wrr=wr

;	figure out default output grid
if nbin gt 0 then begin
  mbin=nbin
  if n_elements(wrr) eq 0 then wrr=[min(xx),max(xx)]
  wmin=min(wrr,max=wmax)
  dw=float(wmax-wmin)/float(mbin)
  yydef=[wmin+findgen(mbin)*dw,wmax]
endif else begin
  mbin=-nbin & wmin=min(alog10(abs(wrr)>(1e-10)),max=wmax)
  dw=(wmax-wmin)/mbin
  yydef=[wmin+findgen(mbin)*dw,wmax] & yydef=10.^(yydef)
endelse

;	now figure out which output grid to use
case ny of
  0: begin				;(case 3 -- not case 3 of X
     yy=yydef & ny=n_elements(yy) & nbin=ny-1L
  end					;NY=0)
  1: begin				;(case 1
    nbin=long(y(0))
    if nbin eq 0 then nbin=1L	;NBIN cannot be 0
    if nbin gt 0 then begin
      mbin=nbin & wmin=min(wrr,max=wmax)
      dw=float(wmax-wmin)/float(mbin)
      yy=[wmin+findgen(mbin)*dw,wmax]
    endif else begin
      mbin=-nbin & wmin=min(alog10(abs(wrr)),max=wmax)
      dw=(wmax-wmin)/mbin
      yy=[wmin+findgen(mbin)*dw,wmax] & yy=10.^(yy)
    endelse
    ny=n_elements(yy) & nbin=ny-1L
    ;if nbin eq 0 then begin	;(case 3 above
    ;  ny=n_elements(yy)
    ;  if ny eq 1 then nbin=yy(0) else nbin=ny-1
    ;endif else yy=y		;NBIN=0)
  end					;NY=1)
  else: begin				;(case 2 -- also case 3 of X
    yy=y & nbin=ny-1 & wrr=[yy(0),yy(ny-1)]
    ;if keyword_set(wrange) then message,'ignoring specified WRANGE',/info
  end					;NY>1)
endcase
wr=wrr

;if ny eq 1 then begin			;(make regular grid
;  if nbin lt 0 then begin		;(log grid
;    mbin=-nbin & wmin=min(alog10(abs(wr)),max=wmax)
;    dw=(wmax-wmin)/mbin
;    yy=[wmin+findgen(mbin)*dw,wmax] & yy=10.^(yy)
;  endif else begin			;)(linear grid
;    mbin=nbin & wmin=min(wr,max=wmax)
;    dw=float(wmax-wmin)/float(mbin)
;    yy=[wmin+findgen(mbin)*dw,wmax]
;  endelse				;NBIN>0)
;endif					;NY=1)

;	get bin widths
ff=f & dx=xx(1:*)-xx & dy=yy(1:*)-yy

;	make sure everything's in the right order
o1=where(dx lt 0,mo1)
if mo1 gt 0 then begin
  oz=sort(xx) & xx=xx(oz) & ff=ff(oz) & dx=xx(1:*)-xx
endif
o1=where(dy lt 0,mo1)
if mo1 gt 0 then begin
  oz=sort(yy) & yy=yy(oz) & dy=yy(1:*)-yy
endif

;	and correct for stupid user tricks
o1=where(dx eq 0,mo1)
if mo1 gt 0 then begin
  o2=where(dx gt 0,mo2)
  if mo2 eq 0 then begin
    message,"Give us this day sliced bread ..er, array for gauss' sake",/info
    return,total(ff)+0*yy
  endif
  message,'forcing non-uniqueness on the absissae',/info
  mdx=min(dx(o2))
  dx(o1)=((spline(o2,dx(o2),o1,0.1)<(mdx/10.))>(mdx/1e6))
  ;dx(o1)=min(dx(o2))/10.
  for i=1L,mo1-1L do xx(o1(i))=xx(o1(i)-1)+dx(o1(i)-1)
endif

if not keyword_set(perbin) then ff=ff*dx		;(units correction

;	make cumulative function
top=total(ff) & cf=0*xx+0*ff(0)
for i=0L,nf-1 do cf(i+1)=cf(i)+ff(i)
base=min(cf)

cg=interpolate(cf,findex(xx,yy))			;interpolate
	;(slower) cg=interpol(cf,xx,yy)			;interpolate
cg = ((cg > base) < top)		;extrapolate

;	undo accumulation
gg=cg(1:*)-cg

y=yy					;output
if not keyword_set(perbin) then gg=gg/dy		;units correction)

return,gg
end
