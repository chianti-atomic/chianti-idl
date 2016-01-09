
PRO burly_fit_fn, x, a, f, pder

;+
; EXPLANATION
;
;     This routine is called from burly_ups.pro. Based on the data in X, 
;     A and the common block it constructs a function F and an array of 
;     partial derivatives, PDER. It used to optimise the spline values 
;     that reproduce the original collision data, through the IDL routine 
;     curvefit.pro.
;
;     This routine should be called by curvefit with the /NODERIV keyword 
;     set. This forces curvefit to compute the partial derivatives 
;     itself. If /NODERIV is not set, then I compute the partial 
;     derivatives in a basic way below (see PROGRAMMING NOTES), but this 
;     should be avoided if at all possible.
;
; INPUTS
;
;     X    Values at which F is defined, the X are expected to be scaled 
;          temperatures.
;
;     A    The spline values.
;
; OUTPUTS
;
;     F    A function defined at the values X and dependent on the A 
;          parameters.
;
;     PDER The partial derivatives dF/dA, see programming note below.
;
; PROGRAMMING NOTES
;
;     The partial derivatives are derived simply by increasing A(i) to 
;     1.001*A(i) and seeing what the change in F(A,X) is.
;
; HISTORY
;
;     Ver.1, 27-Feb-01, Peter Young
;     Ver.2, 8-May-05, Peter Young
;          9 point splines no longer hard-coded
;     Ver.3, 11-Jul-06, Peter Young
;          Added special case for 11 point splines.
;
; CONTACT
;
;     Peter Young, Rutherford Appleton Laboratory, UK
;-


nx=n_elements(x)
na=n_elements(a)

IF na EQ 15 THEN BEGIN
  xs1=dindgen(5)*0.05
  xs2=dindgen(5)*0.125 + 0.25
  nodes=[xs1,xs2,xs1+0.8]
;  nodes=[0,0.125/2.,0.125,0.125*3./2.,0.25,0.375,0.50,0.625,0.75,0.875,1.0]
ENDIF ELSE BEGIN
  nodes=dindgen(na)/(na-1)
ENDELSE


; f is simply obtained by running a spline through the spline values, 
; and obtaining the values at X
;
y2=spl_init(nodes,a)
f=spl_interp(nodes,a,y2,x)

b=a
IF n_params() GT 3 THEN BEGIN
  pder=dblarr(nx,na)
  FOR i=0,na-1 DO BEGIN
    b[i]=a[i]*1.001
    y2=spl_init(nodes,b)
    f2=spl_interp(nodes,b,y2,x)
    pder[*,i]=(f2[i]-f[i])/(0.001*a[i])
  ENDFOR
ENDIF

END
