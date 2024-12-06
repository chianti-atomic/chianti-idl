
PRO ION_INTERP, T, G_T, TI, G_TI, N, PLOT=PLOT

;+
; NAME
;
;       ION_INTERP
;
; PROJECT
;
;       CHIANTI
;
; EXPLANATION:
;
;	When accessing ionisation balance data from the CHIANTI database, it 
;	is output at log temperature values from 4 to 8 at 0.1 intervals. 
;	For display purposes it is often nice to perform spline interpolation 
;	to give a smoother function. This interpolation should be performed 
;	on the logarithm of the data as this is more slowly varying.
;
;	This routine performs such an interpolation. It is only essential to 
;	input G_T and N. Specifying T is optional (see below), while TI and 
;	G_TI are output. N gives the number of points within an interval over 
;	which interpolation is performed.
;
; INPUTS:
;
;       T:      Log10 temperatures at which G_T is specified. Must be
;               same size as G_T.
;	G_T:	Data for which interpolation required
;
; OPTIONAL INPUTS:
;
;	N:	Specifies 'smoothness' of fit. E.g., N=5 means that single 
;		intervals are divided into 5 intervals.
;	TI:	If some temperatures are specified, then the SPLINE call is 
;		used to work out the values of G_TI that correspond to these 
;		temperatures. Only works if N is not specified.
;
; OUTPUTS
;
;	TI:	Interpolated T values. If the TI are specified on the command 
;		line, then these values are returned. 
;	G_TI:	Interpolated G_T values
;
; KEYWORDS
;
;	PLOT -	The interpolated data is displayed together with the original 
;		data.
;
; EXAMPLE:
;
;	g_t=g_of_t(26,13,/ray)
;	ion_interp,t,g_t,ti,g_ti,5,/plot
;	plot,t,g_t,psym=2,xra=[5.8,6.7]
;	oplot,ti,g_ti
;
;	ion_interp,t,g_t,[6.112,6.254],g_ti,/plot
;
; CALLS:
;
;	SPLINE, ION_FRAC_INTERP
;
; HISTORY:
;
;	Ver 1, PRY 17-JUN-98
;	Ver 2, PRY 10-AUG-98, corrected minor bug
;	Ver 3, PRY 22-SEP-98, corrected interpolation error
;	Ver 4, PRY 3-MAR-99,  now allows TI to be specified
;       Ver 5, PRY 9-Mar-12, now calls ION_FRAC_INTERP when N is not
;               specified; T can no longer be used as an output; added
;               some checks on inputs; switched to using spl_init and
;               spl_interp instead of spline.
;
; CONTACT:
;
;	Peter Young, George Mason University (pyoung9@gmu.edu)
;-

IF N_PARAMS() LT 4 THEN BEGIN
  PRINT,'Use:  IDL> ion_interp, t, g_t, ti, g_ti, [ n, /plot]'
  print,'   t - input temperatures'
  print,'   g_t - input ion fractions'
  print,'   g_ti - output ion fractions'
  print,'  either set number of spline points (n)'
  print,'    or   set temperature array (ti)'
  RETURN
ENDIF

IF n_elements(n) NE 0 THEN BEGIN
  IF n NE fix(n) THEN BEGIN
    print,'%ION_INTERP: the input N must be an integer. Returning...'
    return
  ENDIF 
ENDIF

nt=n_elements(t)

;
; The following checks to make sure temperatures are evenly spaced.
;
ints=t[1:nt-1]-t[0:nt-2]
ints_mean=mean(ints)
getmax=max(abs(ints-ints_mean))
IF getmax/ints_mean GE 0.001 THEN BEGIN
  print,'%ION_INTERP: the input temperatures are spaced unequally. This routine will not work with this data.'
  print,'Returning...'
  return
ENDIF 

;---------------------<
IF n_elements(t) NE n_elements(g_t) THEN BEGIN
  print,'%ION_INTERP: the number of elements of T must be the same as those of G_T'
  print,'       n_elements(t)=  '+trim(n_elements(t))
  print,'       n_elements(g_t)='+trim(n_elements(g_t))
  print,'Returning...'
  return
ENDIF 
;---------------------<

;-------------------------------------[0]
; Extract some basics; y contains only non-zero elements of g_t
;                      x contains corresponding t values
;
ngt=N_ELEMENTS(g_t)
y=g_t(WHERE(g_t NE 0.))      ; define y
x=t(WHERE(g_t NE 0.))        ;    and x
;
n_1=WHERE(g_t EQ y(0)) & n_1=n_1(0)
;
npts=N_ELEMENTS(y)
;-------------------------------------[0]


IF keyword_set(plot) THEN BEGIN
  PLOT,x,ALOG10(y),psym=7,xsty=1, $
       xtit='Log!d10!n ( Temperature / K )', $
       ytit='Log!d10!n ( Ionization fraction )'
ENDIF 

;-----------------------------------------------<I>
; If n has been specified, then interpolate g_t over small intervals.
;
IF N_ELEMENTS(n) NE 0 THEN BEGIN
  xmin=min(x)
  xmax=max(x)
  dx=(x[1]-x[0])/float(n)
  n_xi=round((xmax-xmin)/dx)
  xi=findgen(n_xi+1)*dx+xmin

;  xi=FINDGEN( ROUND( (x(npts-1)-x(0))*10.*n +1) ) /10./n + x(0)
  y2=spl_init(x,alog10(y))
  yi=spl_interp(x,alog10(y),y2,xi)
;  yi=SPLINE(x,ALOG10(y),xi)                    ; evaluate yi
 ;
  IF KEYWORD_SET(plot) THEN BEGIN
    OPLOT,xi,yi
  ENDIF
 ;
  g_ti=fltarr((ngt-1)*n +1)
  ti=g_ti
  ni1=n_1*n  &  ni2=n_1*n + (npts-1)*n
  g_ti(ni1:ni2)=10.^yi
  ti=findgen((ngt-1)*n+1)/(ngt-1)/n*(t(ngt-1)-t(0)) + t(0)
 ;
ENDIF
;-----------------------------------------------<I>


;-----------------------------------------------(H)
; n not specified
;
IF N_ELEMENTS(n) EQ 0 THEN BEGIN
  nti=N_ELEMENTS(ti)
 ;
  IF nti EQ 0 THEN BEGIN
    PRINT,' ** Please specify N, or give a non-zero TI **'
    RETURN
  ENDIF 

  g_ti=ion_frac_interp(10.^ti,t,g_t)

  IF KEYWORD_SET(plot) THEN BEGIN
    OPLOT,ti,alog10(g_ti)
  ENDIF
 ;
ENDIF
;-----------------------------------------------(H)

END
