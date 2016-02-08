;+
; HISTORY
;
;     Ver.1, 18-May-2013, Ken Dere
;-
function chianti_ion::CONSTANTS
  ;
  ; cgs
  constants_str = {light:29979245800.d, $ 
  boltzmann:1.3806504d-16, $
  planck:6.6260693d-27, $
  emass:9.10938215d-28, $
  pi:3.1415926535897931d}
  return, constants_str
end