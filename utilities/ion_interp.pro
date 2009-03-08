
PRO ION_INTERP, T, G_T, TI, G_TI, N, PLOT=PLOT

;+
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
;	G_T:	Data for which interpolation required
;
; OPTIONAL INPUTS:
;
;	T:	Values at which is G_T is tabulated. If not specified, then 
;		taken as 4 to 8 in 0.1 intervals.
;	N:	Specifies `smoothness' of fit. E.g., N=5 means that single 
;		intervals are divided into 5 intervals.
;	TI:	If some temperatures are specified, then the SPLINE call is 
;		used to work out the values of G_TI that correspond to these 
;		temperatures. Only works if N is not specified.
;
; OUTPUTS
;
;	T:	If not specified, then output as an array from 4 to 8 in 0.1 
;		intervals
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
;	SPLINE
;
; HISTORY:
;
;	Ver 1, PRY 17-JUN-98
;	Ver 2, PRY 10-AUG-98, corrected minor bug
;	Ver 3, PRY 22-SEP-98, corrected interpolation error
;	Ver 4, PRY 3-MAR-99,  now allows TI to be specified
;
; CONTACT:
;
;	Peter Young, Cambridge University (P.R.Young@damtp.cam.ac.uk)
;-

IF N_PARAMS() LT 4 THEN BEGIN
  PRINT,'Use:  IDL> ion_interp, t, g_t, ti, g_ti, [ n, /plot]
  RETURN
ENDIF

;---------------------<
; If T not specified, then assume that we have the standard CHIANTI 
; temperature scale
;
IF N_ELEMENTS(t) EQ 0 THEN t=findgen(41)/10. +4.
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


;-----------------------------------------------<I>
; If n has been specified, then interpolate g_t over small intervals.
;
IF N_ELEMENTS(n) NE 0 THEN BEGIN
  xi=FINDGEN( ROUND( (x(npts-1)-x(0))*10.*n +1) ) /10./n + x(0)
  yi=SPLINE(x,ALOG10(y),xi)                    ; evaluate yi
 ;
  IF KEYWORD_SET(plot) THEN BEGIN
    PLOT,x,ALOG10(y),psym=2,xsty=1
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
 ;
 ; note: I've had to add the 0.001's below because IDL dosen't always think 
 ; that, e.g., 5.80=5.80 !
  ind=WHERE( (ti GE x(0)-0.001) AND (ti LE x(npts-1)+0.001) )
  IF ind(0) EQ -1 THEN BEGIN
    g_ti=ti*0.    ; set all gti's to zero
    RETURN
  ENDIF ELSE BEGIN
    g_ti=MAKE_ARRAY(nti,value=0.)
    tix=ti(ind)
  ENDELSE
 ;
  yi=SPLINE(x,ALOG10(y),tix)        ; evaluate yi
  g_tix=10.^yi
  g_ti(ind)=g_tix
 ;
  IF KEYWORD_SET(plot) THEN BEGIN
    PLOT,x,ALOG10(y),xsty=1
    OPLOT,ti,yi,psym=4,symsiz=2
  ENDIF
 ;
ENDIF
;-----------------------------------------------(H)

END
