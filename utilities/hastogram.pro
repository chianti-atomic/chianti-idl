function hastogram,list,x,wts=wts, _extra=e
;+
;function	hastogram
;	returns a frequency histogram over a specified grid,
;	calculated in a fast and clever way.  really.
;	works in haste, without waste.
;
;syntax
;	f=hastogram(list,x,wts=wts)
;
;parameters
;	list	[INPUT; required] list of numbers to bin
;	x	[I/O] the required binning scheme
;		1: if not set, assumes a linear grid of size 100 from
;		   min(LIST) to max(LIST)
;		2: if scalar or 1-element vector, assumes this to be the
;		   number of bins
;		   * if -ve, log grid, else linear
;		3: if vector with more than 1 element, assumes it to be
;		   all the bin-beginning values and the final bin-ending
;		   value.
;		*  overwrites variable by adopted grid
;
;keywords
;	wts	[INPUT] if given, appropriately weights each element of LIST
;		* default is unity
;		* if size does not match size(LIST), will be appropriately
;		  interpolated
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	the problem is the following: if the number of bins is large,
;	accumulating a frequency histogram by linear search takes too
;	long, esp. in IDL if one uses for-loops.
;	so, first create a monster array containing both the list elements
;	and the grid values.  then sort this array.  this results in an
;	array where the list elements are all nicely placed within the
;	correct bins.  now, if we've been keeping track of the positions
;	of the list elements, it's an easy job to count up the number in
;	each bin of the grid.  to do the latter, we first create a new
;	array made up of -1s in positions of list elements and position
;	indices for grid values and reorder according to the above sort.
;	then, replace each -1 by the nearest non-(-1) from the left.  now
;	each list element is assigned the correct bin number.  voila!
;
;restrictions
;	X must be sorted in increasing order.
;
;subroutines
;
;
;history
;	vinay kashyap (Jan98)  PINT_of_ALE
;	added kludge to speed up in case max(X) < max(LIST) (VK; Sep98)
;	added quit in case X is not sorted in ascending order (VK; JanMMI)
;-

;	usage
nn=n_elements(list)
if nn eq 0 then begin
  print,'Usage: f=hastogram(list,x,wts=wts)'
  print,'  returns frequency histogram accumulated over irregular grid'
  return,-1L
endif

;	decode X
nx=n_elements(x)
case nx of					;{decode X
  0: begin				;(use hardcoded default
     xmin=min(list,max=xmax) & nbin=100L & dx=float(xmax-xmin)/float(nbin)
     xx=[xmin+findgen(nbin)*dx,xmax]
  end					;NX=0)
  1: begin				;(number of bins given
     nbin=x(0)
     if nbin eq 0 then nbin=100L
     if nbin lt 0 then begin		;(log grid
       xmin=min(alog10(abs(list)),max=xmax)
       dx=float(xmax-xmin)/float(nbin)
       xx=[xmin+findgen(nbin)*dx,xmax] & xx=10.^(xx)
     endif else begin			;)(linear grid
       xmin=min(list,max=xmax)
       dx=float(xmax-xmin)/float(nbin)
       xx=[xmin+findgen(nbin)*dx,xmax]
     endelse				;)
  end					;NX=1)
  else: xx=x				;all is well
endcase						;X}
nx=n_elements(xx)

;	check to see that XX are sorted in increasing order
if xx(1) lt xx(0) then begin
  message,'binning grid must be sorted in increasing order',/info
  return,-1L
endif

;	poundage
nw=n_elements(wts) & ww=intarr(nn)+1		;default=1
if nw eq 1 then ww(*)=wts(0)			;all set to same
if nw gt 1 and nw ne nn then $			;interpolate over range
	ww=interpol(wts,findgen(nw),findgen(nn)*(nw-1.)/(nn-1.))
if nw eq nn then ww=wts				;specified

;	declare the output
h=intarr(nx-1) +$			;regular frequency histogram
	0*ww(0)				;changed to appropriate type

;	find the indices
allX=[xx(0),list,xx(1:*)]	;the large array containing LIST and grid
	;NOTE: XX is split so that min(LIST)=xx(0) is not lost in the sort
oX=sort(allX)		;sort above to move list elements into correct bins
allI=[0L,lonarr(nn)-1L,lindgen(nx-1)+1]	;the position remembering array
allW=[0,ww,intarr(nx-1)]		;and the weight remembering array
allI=allI(oX)			;positions have been reordered per above sort
allW=allW(oX)				;same with the weights
    ;kludge in case max(X) < max(LIST)
    ugh=where(allI ge 0,mugh)
    if ugh(mugh-1L) lt nn+nx-1 then allI(ugh(mugh-1L)+1:*)=allI(ugh(mugh-1L))
oI=where(allI lt 0,moI)			;remember where the LIST elements are
oJ=oI					;copy of same, will be modified
while moI gt 0 do begin		;{while there are list elements w/o bin indices
;  kilroy; was here - commented out -- GDZ
  allI(oJ)=allI(oJ-1L)			;assign bin value
  moJ=moI & oJ=where(allI lt 0,moI)	;any left?
  if moI eq moJ then moI=0		;in case XMIN > min(LIST)
endwhile			;MOI>0}

;	extract the indices
oo=allI(oI)			;these are the bin indices of LIST elements
f=allW(oI)			;extracting the weights in the right order

;	now make the frequency histogram
;	we expect to use this routine when the number of bins are much
;	larger than the number of LIST elements, i.e., when very few of
;	the bins are actually populated.  so minimal qualms about stepping
;	through only the populated bins using a (ugh) for-loop
oq=oo(uniq(oo,sort(oo))) & noq=n_elements(oq)	;all the populated bins
for i=0L,noq-1L do begin		;{for each populated bin
  k=oq(i) & ok=where(oo eq k)		;find the LIST elements
  if k ge 0 and k lt nx-1 then h(k)=total(f(ok))	;and approp. weights
endfor					;I=0,NOQ-1}

;	outputs
x=xx
return,h

end
