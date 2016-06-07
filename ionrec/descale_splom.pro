
 ;+
; PROJECT     :  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory, George Mason University (USA), the University of
;	Florence (Italy), and the University of Cambridge (UK).
;       
;
; NAME
;
;     DESCALE_SPLOM	
;
;  PURPOSE:
;
;     This routine returns ionization energies and cross sections from scaled ionization energies 
;     and cross sections using a Burgess-Tully (1992) type scaling developed for ionization
;     
;     
;
; INPUTS
;
;
;    SPLOM     A structure returned when a .splom file is read by read_splups
;    E_IN      The incident energy (Rydbergs) for which the collision strength is returned
;    
;
; OUTPUTS
;
;    OMEGA    The collision strength.
;
; OPTIONAL INPUTS
;
;    NONE 
;
; KEYWORDS
;
;    NONE 
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
;
;    *************************************************
;
pro descale_splom,splom,e_in,omega
;
;  to return collision strengths from scaled collision strength,
;   given the energy e_int
;
;  splom is a structure created by read_splups
;  splom.lvl1, splom.lvl2, splom.t_type, splom.gf, splom.de, 
;      splom.c_ups, splom.nspl, splom.spl
;
nsplom=n_elements(splom.lvl1)
;
omega=fltarr(n_elements(e_in),nsplom)
;
;  to do spline interpolation, e_in must be monotonically increasing
;
esrt=sort(e_in)
e_in=e_in(esrt)
index=indgen(n_elements(e_in))
unsrt=index(esrt)
;

;
for isplom=0,nsplom-1 do begin
	;
	de_in=splom.de
	de_in=de_in(isplom)*13.6056923
	c_ups=splom.c_ups
	c_ups=c_ups(isplom)
	c_curr=c_ups
	t_type=splom.t_type
	t_type=t_type(isplom)
	;
	nspl=splom.nspl
	nspl=nspl(isplom)
	dxs=1./float(nspl-1.)
	sx=dxs*findgen(nspl)
	som=splom.spl
	som=som(0:nspl-1,isplom)
	gen=where(e_in ge de_in,nge)
	if nge le 0 then return	
    x_int=e_in(gen)/de_in
	;
;
CASE t_type OF

1:  begin
      sx_int=1.-alog(c_curr)/alog(x_int-1.+c_curr)
      som2=nr_spline(sx,som)
      som_int=nr_splint(sx,som,som2,sx_int)
	  for ige=0,nge-1 do begin
         omega(gen(ige),isplom)=som_int(ige)*alog(x_int(ige)-1.+exp(1.))
	  endfor
    end
;
2:  begin
      sx_int=(x_int-1.)/(x_int-1.+c_curr)
      som2=nr_spline(sx,som)
	  for ige=0,nge-1 do begin
         omega(gen(ige),isplom)=nr_splint(sx,som,som2,sx_int(ige))
	  endfor
    end
;
3:  begin
      sx_int=(x_int-1.)/(x_int-1.+c_curr)
      som2=nr_spline(sx,som)
      som_int=nr_splint(sx,som,som2,sx_int)
	  for ige=0,nge-1 do begin
         omega(gen(ige),isplom)=som_int/x_int(ige)^2
	  endfor
    end
;
4:  begin
      sx_int=1.-alog(c_curr)/alog(x_int-1.+c_curr)
      som2=nr_spline(sx,som)
      som_int=nr_splint(sx,som,som2,sx_int)
	  for ige=0,nge-1 do begin
         omega(gen(ige),isplom)=som_int(ige)*alog(x_int(ige)-1.+c_curr)
	  endfor
   end
;
else:  print,' for descale_splom, t_type ne 1,2,3,4=',t_type
;
ENDCASE
;
endfor  ;  loop over nsplom
;
e_in=e_in(unsrt)
for isplom=0,nsplom-1 do begin
	omega(0,isplom)=omega(unsrt,isplom)
endfor
;
return
end
;
;   ***************************************************
