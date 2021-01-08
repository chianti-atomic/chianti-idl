
;+
; NAME:
;     CHIANTI_HTML
;
; PURPOSE:
;     Creates a series of html files and tables, organised by element,
;     under the 'your_top_WWW_page' directory.
;
;     The files are created locally and should be ftp'ed to the
;     CHIANTI webpage site.
;
;     When asked to create the index.html file, you should type 'y'
;     for yes.
;
; CATEGORY:
;     CHIANTI; webpage.
;
; CALLING SEQUENCE:
;     CHIANTI_HTML
;
; INPUTS:
;     None.
;
; OPTIONAL INPUTS:
;     Directory:  The top directory where the files will go. If not
;                 set, then the files will go in the sub-directory
;                 "tree" of the current working directory.
;     Www_Address: The http address of the 'dbase' directory
;                  containing the Solarsoft database directories. If
;                  not specified, then the NASA sohoftp address is
;                  used.
;     XuvTop:  The path to the CHIANTI dbase directory on the
;              user's computer. If not specified, then the
;              system variable !xuvtop is used.
;	
; KEYWORD PARAMETERS:
;     None.
;
; OUTPUTS:
;     Creates a sub-directory 'tree' that contains further
;     sub-directories for each element. A html file for each ion is
;     created containing the file comments for the ion's data
;     file and stored in the element directory. A file tree/index.html
;     is also created containing a table with links to all the ion
;     files. 
;
; EXAMPLE:
;     IDL> chianti_html
;
; MODIFICATION HISTORY:
;     Ver.1, Aug-2002, GDZ
;     Ver.2, Sep-2015, GDZ
;       Modified for v.8.
;     Ver.3, 17-Nov-2020, GDZ
;       Removed the 'dielectronic' table.
;     Ver.4, 18-Nov-2020, GDZ
;       Added a warning and removed a double if.
;     Ver.5, 06-Jan-2021, Peter Young
;       Changed default value of WWW_ADDRESS to point to Solarsoft;
;       updated header.
;-



PRO chianti_html, directory =directory, www_address=www_address,  xuvtop =xuvtop 


  IF n_elements(directory) EQ 0 THEN directory= 'tree' 


 ;
 ; The following means that links to the data files go to the
 ; SolarSoft master tree.
 ;
  IF n_elements(www_address) eq 0 then www_address='https://sohoftp.nascom.nasa.gov/solarsoft/packages/chianti/dbase/'
  
  IF n_elements(xuvtop) EQ 0 THEN   xuvtop = !xuvtop

  CASE os_family() OF
     'vms': path_sep = ','
     'Windows': path_sep = ';'
     ELSE: path_sep = ':'
  ENDCASE


  dirs_elements = ['al', 'ar', 'b','be','c', 'ca', 'cl', 'co', 'cr', $
                   'cu','f', 'fe', 'h', 'he', 'k', 'li','mg', 'mn', 'n', 'na', $
                   'ne', 'ni', 'o', 'p', 's','sc', 'si','ti', 'v','zn']

  top_dirs = '+'+xuvtop+'/'+dirs_elements(0)
  FOR el = 1, N_ELEMENTS(dirs_elements)-1 DO $ 
     top_dirs = top_dirs+path_sep+'+'+xuvtop+'/'+dirs_elements(el)

  all_sub_dirs = expand_path(top_dirs, /all_dirs, /ARRAY)
  break_file,  all_sub_dirs, disk_log, dir, ionst, ext


  element=['H','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',$
           'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn']

  z_lbl=['h','he','li','be','b','c','n','o','f','ne','na','mg','al','si','p','s','cl',$
         'ar','k','ca','sc','ti','v','cr','mn','fe','co','ni','cu','zn']

  ionstage=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV',$
            'XVI','XVII','XVIII','XIX','XX','XXI',' XXII','XXIII','XXIV','XXV','XXVI','XXVII',$
            'XXVIII']           ;,'XXIX','XXX','XXXI']


;      snote=element(iz-1)+' '+ionstage(ion-1)

  yes_no, ' create the index.html file with the front page containing all ions ? ', yesno, 'N'

  IF yesno THEN BEGIN 

     openw, luo, directory+'/'+'index.html', /get_lun


     openr, 1, xuvtop+'/VERSION'
     version = ''
     readf, 1, version
     close, 1

     printf, luo, '<center><h2> VERSION '+version+' </h2></center>'
     printf, luo, '<p><hr><p>'

     printf, luo, 'Below are links to the data, references, and tables of most prominent lines. <BR>'

     printf, luo, '<h3> Different colours indicate different isoelectronic sequences.</h3>'
     

     printf, luo, '<font size=-2>'

     printf, luo, '<table  WIDTH="100%" NOSAVE >'
     printf, luo, '<tr><td></td>'
     text = ''

     FOR  i=0, n_elements(ionstage)-1 DO text = text+ '<td><B>  </B></td> '
     printf, luo, text
     printf, luo, '</tr> '


     FOR  el=0, n_elements(element)-1 DO BEGIN 

;first check that we have the element:
;-------------------------------------
        if dir_exist(xuvtop+'/'+trim(strlowcase(element(el)))) then begin 

           printf, luo, '<tr>'
           printf, luo, '<td><B> '+element(el)+' </B></td>'

;loop throught the ions. 
;------------------------

           FOR j=0, n_elements(ionstage)-1 DO BEGIN 

              iz = el+1 &  ion=j+1

              ion_string=strtrim(z_lbl(iz-1),2)+'_'+strtrim(string(ion,'(i2)'),2)

              dir_name=xuvtop+'/'+strtrim(z_lbl(iz-1),2)+'/'+$
                       ion_string

              wname=dir_name+'/'+ion_string+'.wgfa'

;check if we have the ion. 
;-------------------------        
              
              index = where(ionst EQ ion_string, nn)

;we have the directory 
              IF nn EQ 1 THEN BEGIN 

;check that we have at least the main file:

                 IF  file_exist(wname) then begin 

;         filename = directory+'/'+trim(z_lbl(iz-1))+'/'+$
;           trim(z_lbl(iz-1))+'_'+trim(ion)+'.html'

                    filename = trim(z_lbl(iz-1))+'/'+$
                               trim(z_lbl(iz-1))+'_'+trim(ion)+'.html'

                    text_color = '' 
                    isoel = iz-ion+1

                    IF isoel EQ 1 THEN text_color = 'BGCOLOR=yellow'
                    IF isoel EQ 3 THEN text_color = 'BGCOLOR=lightgreen'
                    IF isoel EQ 5 THEN text_color = 'BGCOLOR=orange'
                    IF isoel EQ 7 THEN text_color = 'BGCOLOR=red'
                    IF isoel EQ 9 THEN text_color = 'BGCOLOR=lightblue'
                    IF isoel EQ 11 THEN text_color = 'BGCOLOR=magenta'
                    IF isoel EQ 13 THEN text_color = 'BGCOLOR=pink'
                    IF isoel EQ 15 THEN text_color = 'BGCOLOR=lightgreen'
                    IF isoel EQ 17 THEN text_color = 'BGCOLOR=yellow'
                    IF isoel EQ 19 THEN text_color = 'BGCOLOR=red'


                    printf, luo, '<td '+text_color+' > <a href="'+filename+'"> '+ionstage(j)+' </a></td>' 

                 endif  else begin 
                    
                    print,'no wgfa file?? ', wname
                                ;   stop
                    
                    printf, luo, '<td>  </a></td>'
                 end 

              ENDIF  ELSE     printf, luo, '<td>   </a></td>'

           ENDFOR               ;loop  throught the ionization stages
           printf, luo, '</tr> '

        endif else begin 

           print,'Element: '+element(el)+' not present'

        end 

     ENDFOR


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     printf, luo, '<em> Page designed and last modified by '+$
             ' Giulio Del Zanna on '+ systime() +'</em>'


     free_lun, luo

;    stop
  ENDIF  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ioneq_name=!ioneq_file
  read_ioneq,ioneq_name,temp_all,ioneq,ioneq_ref

  
  FOR el = 0, N_ELEMENTS(dirs_elements)-1 DO BEGIN

     dirs_ions = expand_path('+'+xuvtop+'/'+dirs_elements(el), /all_dirs, /ARRAY)
;remove the top directory
     nn =  N_ELEMENTS(dirs_ions) 

     IF nn EQ 1 THEN GOTO, finish

     dirs_ions = dirs_ions (0:nn-2)

     FOR io=0, N_ELEMENTS(dirs_ions)-1 DO BEGIN 

        break_file, dirs_ions(io) , disk_log, dir, ion_string, ext

        convertname,ion_string , iz ,ion

        ion2spectroscopic, ion_string ,snote, dielectronic=dielectronic
        
        wname=dirs_ions(io)+'/'+ion_string+'.wgfa'
        elvlcname=dirs_ions(io)+'/'+ion_string+'.elvlc'
        upsname=dirs_ions(io)+'/'+ion_string+'.scups'
        pupsname = dirs_ions(io)+'/'+ion_string+'.psplups'


        IF  file_exist(wname) AND  file_exist(elvlcname) AND  $
           file_exist(upsname)  THEN BEGIN 

           print,'doing ', ion_string
           
;
;  read energy levels
;
           ;; read_elvlc_check, elvlcname, level, config, s_values, spd_values, j_values, $
           ;;   ecm,eryd,ecmth,erydth, eref,$
           ;;   err=err, outfile=outfile, luo=luo
           
           read_elvlc, elvlcname, level, term, config, s_values, spd_values, j_values, $
                       ecm,eryd,ecmth,erydth, eref, elvlcstr=elvlcstr
           
; read things properly:
           config=elvlcstr.data.conf
           spd_values=elvlcstr.data.l_sym
           j_str=elvlcstr.data.j_str
           
;
;   read electron collision data
;
           read_scups,upsname,splstr
           upsref =splstr.info.comments

;
;  read in  wavelengths, gf and A values from .wgfa files
;

           read_wgfa2,wname,lvl1,lvl2,wvl1,gf1,a_value1,wgfaref

           if not dir_exist(directory+'/'+dirs_elements(el)) then $
              mk_dir, directory+'/'+dirs_elements(el)

                                ;           openw, luo, directory+'/'+dirs_elements(el)+'/'+ion_string+'.ht', /get_lun
                                ;           print, 'printing '+directory+'/'+dirs_elements(el)+'/'+ion_string+'.ht'

           openw, luo, directory+'/'+dirs_elements(el)+'/'+ion_string+'.html', /get_lun
           print, 'printing '+directory+'/'+dirs_elements(el)+'/'+ion_string+'.html'
           
           
           openr, 1, xuvtop+'/VERSION'  &  version = ''
           readf, 1, version &  close, 1
           printf, luo, '<center><h2>'+ snote+' -  VERSION '+version+' </h2></center>'
           printf, luo, '<p><hr><p>'
           printf, luo, 'The CHIANTI  database has the following primary ASCII files for this ion: '
           printf, luo, ' <p>'
           break_file, elvlcname , disk_log, dir, name , ext
           elvlcname=www_address+$
                     dirs_elements(el)+'/'+ name+'/'+name+ext


           printf, luo, '<ol>'
           printf, luo, '<h2><li> <a href="'+elvlcname+'"> '+name+ext +' (energy levels)</h2> </a> <p>'
           printf, luo, ' contains the energy levels (in cm-1). '
           printf, luo, '  It includes both experimental  and theoretical values of the levels energies. <p>'
           printf, luo, '<pre>'
           printf, luo, ' '
           FOR n=0, n_elements(eref)-1 DO  printf, luo, eref(n)
           printf, luo, ' '
           printf, luo, '</pre><p>' 


           break_file, wname , disk_log, dir, name , ext

           wn =www_address+$
               dirs_elements(el)+'/'+ name+'/'+name+ext

           printf, luo, '<h2><li> <a href="'+wn+'">'+ name+ext+' (radiative data) </h2> </a>  <p> '
           printf, luo, ' contains  wavelengths, gf and A values of the transitions. ' 
           printf, luo, ' The wavelengths  are based on the experimental energy levels '
           printf, luo, ' and should be the best available.' 
           printf, luo, ' Wavelengths calculated from the theoretical energies are of an' 
           printf, luo, '     indeterminate accuracy and their values are presented as negative'
           printf, luo, '      values of the calculated wavelength. <p>'
           printf, luo, '<pre>'
           printf, luo, ' '
           FOR n=0, n_elements(wgfaref)-1 DO  printf, luo,wgfaref(n)
           printf, luo, '</pre><p>'

           break_file, upsname , disk_log, dir, name , ext

           upsn=www_address+$
                dirs_elements(el)+'/'+ name+'/'+name+ext

           printf, luo, '<h2> <li>  <a href="'+upsn+'"> '+$
                   name+ext+ ' (electron collision data) </h2> </a> <p>'
           printf, luo, ' contains the effective electron collision'
           printf, luo, '       strengths scaled according to the rules formulated by Burgess and'
           printf, luo, '       Tully (1992). ' 
           printf, luo, '<pre>'
           printf, luo, ' '
           
           
           FOR n=0, n_elements(upsref)-1 DO  printf, luo,upsref[n_elements(upsref)-1-n]
           printf, luo, '</pre><p>'


           if  not dielectronic then begin 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              IF file_exist(pupsname) THEN BEGIN 

;
;   read proton collision fits
;
                 read_splups,pupsname,pstr, pupsref, /prot


                 break_file, pupsname , disk_log, dir, name , ext

                 pupsname=www_address+$
                          dirs_elements(el)+'/'+ name+'/'+name+ext

                 printf, luo, ' <h2><li> <a href="'+pupsname+'">  '+$
                         name+ext+ ' (proton collision data) </h2> </a> <p>'
                 printf, luo, ' contains the   spline fits to the scaled proton collision'
                 printf, luo, '       strengths.'

                 printf, luo, '<pre>'
                 printf, luo, ' '
                 FOR n=0, n_elements(pupsref)-1 DO  printf, luo,pupsref(n)
                 printf, luo, '</pre><p>'

              ENDIF ELSE BEGIN 
                 printf, luo, '<h2> <li>  '+ $
                         ' (proton collisional data)  Not available in this VERSION. </h2>  <p>'
              END 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              if  file_exist(dirs_ions(io)+'/'+ion_string+'.reclvl') then begin 

                 read_ionrec,dirs_ions(io)+'/'+ion_string,$ 
                             rec_rate, ci_rate, temp_ionrec, luprec, lupci, status,rec_ref,ci_ref

                 name=www_address+dirs_elements(el)+'/'+ion_string+'/'+ion_string+'.reclvl'
                 printf, luo, ' <h2><li> <a href="'+name+'"> '+$
                         ion_string+'.reclvl'+ ' (total recombination population rates) </h2> </a> <p>'

                 printf, luo, '<pre>' &  printf, luo, ' '
                 FOR n=0, n_elements(rec_ref)-1 DO  printf, luo,rec_ref(n)
                 printf, luo, '</pre><p>'
                 printf, luo, ' '

              endif 

              if   file_exist(dirs_ions(io)+'/'+ion_string+'.cilvl') then begin 

                 read_ionrec,dirs_ions(io)+'/'+ion_string,$ 
                             rec_rate, ci_rate, temp_ionrec, luprec, lupci, status,rec_ref,ci_ref

                 name=www_address+dirs_elements(el)+'/'+ion_string+'/'+ion_string+'.cilvl'
                 printf, luo, ' <h2><li> <a href="'+name+'"> '+$
                         ion_string+'.cilvl'+ ' (collisional ionization population rates) </h2> </a> <p>'

                 printf, luo, '<pre>' &  printf, luo, ' '
                 FOR n=0, n_elements(ci_ref)-1 DO  printf, luo,ci_ref(n)
                 printf, luo, '</pre><p>'
                 printf, luo, ' '

              endif 

           endif  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; calculate emissivities at the peak ion fraction, for a 
; range of densities:

                                ;

           f_all=ioneq(*,iz-1,ion-1+dielectronic)


           ind=WHERE(f_all EQ max(f_all))
           ind=ind(0)           ; in case ind has more than one element (e.g., Ar IX)
           logt=temp_all(ind)

; default density range

           lo_dens=6 & hi_dens=20 & dint=1.0
           n=fix((hi_dens-lo_dens)/dint)+1
           dens=findgen(n)*dint + lo_dens 

; setting logt=log Tmax, we get emissivities for 3 temperatures: log T_max +- 0.15. 

           emiss=emiss_calc(iz,ion,temp=logt,dens=dens, diel=dielectronic,$
                            path=dirs_ions(io), /quiet, ioneq_file=!ioneq_file, $
                            abund_file=!abund_file) ;, /no_de

           cut=1e-3
;;;;;;;;;;;;;;;;;;;;;;;
           
           nd=n_elements(dens)

           dims=n_elements(emiss[0].em[*,0])
           if dims eq 1 then ind_t=0 else if dims eq 3 then ind_t=1 else begin 
              print,'problems with number of temperatures'
              stop
           end 

           index=(WHERE(emiss[*].em[ind_t,0] EQ MAX(emiss[*].em[ind_t,0])))(0)

           for i=0,nd-1 do begin 
;get the index of the strongest line:

              ind_m=WHERE(emiss[*].em[ind_t,i] EQ MAX(emiss[*].em[ind_t,i]))
              ind_m=ind_m[0]

              if emiss[ind_m].em[ind_t,i] ne 0. then begin 

                 ind=WHERE(emiss[*].em[ind_t,i] ge cut*emiss[ind_m].em[ind_t,i],nn)

                 if nn gt 0 then begin 
                    index=[index,ind]
;remove duplicates:
                    index=index(rem_dup(index))
                 end 
              end 
           endfor 

;print,max(index)
;help,emiss,index
;stop

; sort the wavelengths:

           i_sort = sort(emiss[index].lambda)
           index=index[i_sort]


           nl=n_elements(index)

;selected 
           dens_sel=[6,10,15,20]
;*********************

           nsel=n_elements(dens_sel)

           ratios=fltarr(nsel,nl) 

           for l=0,nsel-1 do begin 

              ind_ne=(where(dens_sel[l] eq dens))(0)

;indexes of maximum emissivity:
              ind=WHERE(emiss[*].em[ind_t,ind_ne] EQ MAX(emiss[*].em[ind_t,ind_ne]))

              ratios[l,*]= emiss[index].em[ind_t,ind_ne]/emiss[ind[0]].em[ind_t,ind_ne] 
           endfor 
;stop

           index_w=intarr(nl)

           for i=0L,nl-1 do begin 

              ind=where(lvl1 eq  emiss[index[i]].LEVEL1 and $
                        lvl2 eq emiss[index[i]].LEVEL2,n)

              index_w[i]=ind

           end 

           wvl1 = wvl1[index_w]
           lvl1 = lvl1[index_w]
           lvl2= lvl2[index_w]
           gf1 = gf1[index_w]
           a_value1=a_value1[index_w]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;check how many lines we have:


           IF n_elements(lvl1) LE 200 THEN n_tables = 1 ELSE BEGIN 
              IF n_elements(lvl1)/200. NE n_elements(lvl1)/200 THEN $
                 n_tables =n_elements(lvl1)/200 +1 ELSE n_tables =n_elements(lvl1)/200
           END

           printf, luo, '<p><hr><p>'

           if n_tables eq 1 then $
              printf, luo, '<h3> The  html table below contains the radiative data '+$
                      ' of the brightest lines sorted in wavelength</h3>' else begin 

              printf, luo, '<h3> The  html tables below contain the radiative data '+$
                      'of the brightest lines sorted in wavelength</h3>'
              printf, luo, ' <p> <ul>'
           end 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

           i_start = 0
           i_END = n_elements(lvl1)


           FOR i_tab=0, n_tables-1 DO BEGIN 

              IF n_tables EQ 1 THEN BEGIN 

                 table_name = directory+'/'+dirs_elements(el)+'/'+ion_string+'_table.html' 
                 i_start = 0
                 i_END = n_elements(lvl1)-1

                 printf, luo, ' <a href="'+ion_string+'_table.html'+'">  <h2>'+ ion_string+'_table.html'+$
                         ' (html table) </h2> </a> <p>'

              ENDIF ELSE  BEGIN 

                 table_name = directory+'/'+dirs_elements(el)+'/'+ion_string+'_table'+trim(i_tab+1)+'.html'
                 i_start = 200*i_tab
                 i_END = 200*i_tab+199 
                 IF  i_END GT  n_elements(lvl1)-1 THEN i_END =  n_elements(lvl1)-1


                 tab = ion_string+'_table'+trim(i_tab+1)+'.html'

                 printf, luo, ' <h2><li>  <a href="'+tab+'"> '+ tab+$ 
                         ' (html table -- '+$
                         string(min(abs(wvl1(i_start:i_end))), '(f15.4)')+' -- '+$
                         string(max(abs(wvl1(i_start:i_end))), '(f15.4)')+' A'+$
                         ') </h2> </a> <p>'


              END 

              openw, lup, table_name, /get_lun

              print, 'printing '+table_name

              printf, lup, ' <html> <head>  <title> CHIANTI atomic database</title> '
              printf, lup, '</head>'
              printf, lup, '<BODY  bgcolor=white vlink="#CC33CC">'
              printf, lup, ''
              printf, lup, '<center> <h2><font color="#FF0000" >' 
              printf, lup, ' CHIANTI </font></h2></center>'
              printf, lup, '<center> <h3> An Atomic Database for Spectroscopic Diagnostics' 
              printf, lup, ' of Astrophysical Plasmas.</h3></center> '

              printf, lup, ' '

;      printf, lup, '<pre> http://www.damtp.cam.ac.uk/user/astro/chianti/chianti.html </pre>'

              printf, lup, '<hr><p>'

              openr, 1, xuvtop+'/VERSION'
              version = ''
              readf, 1, version
              close, 1

              printf, lup, '<center><h2>'+ snote+' -  VERSION '+version+' </h2></center>'

              printf, lup, '<p><hr><p>'

              printf, lup, 'Here you find various information on the brightest lines, ordered in wavelength:<BR>'
              printf, lup, 'Wavelengths (Angstroms); weighted oscillator strengths (gf); '
              printf, lup, ' transition probabilities (A); '
              printf, lup, ' intensities (normalised) at '

              text='10<sup>'+trim(dens_sel[0])+'</sup> '
              for l=1,nsel-1 do  text=text+ ', 10<sup>'+trim(dens_sel[l])+'</sup>'

              printf, lup,text
              printf, lup,' cm<sup>-3</sup>; '

              printf, lup, 'together with: <BR> Configuration; Term; '+$
                      ' Index of the level (1=ground); Energy (observed if available - otherwise theoretical); '
              printf, lup, ' of the lower (l) and upper levels (u) <BR> '
              printf, lup, ' Note that intensities have been calculated at log T[K]='+$
                      trim(logt,'(f5.2)')+$
                      ', where the ion fraction peaks in ionization equilibrium.'+$
                      ' Only lines brighter than '+trim(cut)+$
                      ' the strongest line (at all densities) are retained.'

              dummy = where(wvl1(i_start:i_end) LT 0, dd)
              IF dd GT 0 THEN  BEGIN 
                 printf, lup, ' Lines flagged with a "*" before their wavelengths only have theoretical energy levels.'
                 printf, lup, ' Please note that their wavelengths are not very accurate. '
              END 
              printf, lup, '<p><hr><p>'


              printf, lup, '<table BORDER CELLSPACING=5 NOSAVE >'


              printf, lup, '<tr>'
              printf, lup, ' <td><h3> Wavelength  </h3></td> '
              printf, lup, '<td><h3> gf  </h3></td>'
              printf, lup, ' <td><h3> A<sub>ul</sub>  </h3></td> '

              for l=0,nsel-1 do   printf, lup, '<td><h3> R  </h3></td>'

              printf, lup,'<td><h3> C<sub>l</sub> </h3></td> '
              printf, lup, ' <td><h3> T<sub>l</sub>  </h3></td> '
              printf, lup,'<td><h3> C<sub>u</sub> </h3></td> <td><h3> T<sub>u</sub> </h3></td>'
              printf, lup, ' <td><h3> I<sub>l</sub></h3></td> <td><h3> I<sub>u</sub></h3></td>'
              printf, lup, '<td><h3> E<sub>l</sub> </h3></td> <td><h3> E<sub>u</sub>  </h3></td>'   
              printf, lup, '</tr> '
              printf, lup, '<tr>'

              printf, lup, ' <td><h3>  (A) </h3></td>'+$
                      ' <td><h3>  </h3></td>'
              printf, lup, ' <td><h3>  (s-1) </h3></td> '

              for l=0,nsel-1 do   printf, lup, ' <td><h3>  (10<sup>' +$
                                          trim(dens_sel[l])+'</sup>) </h3></td> '
              

              printf, lup,'<td><h3>   </h3></td> <td><h3>  </h3></td>'
              printf, lup, ' <td><h3>   </h3></td> <td><h3>  </h3></td>'
              printf, lup, ' <td><h3>  </h3></td> <td><h3>  </h3></td>'
              printf, lup, '<td><h3>  (cm<sup>-1</sup>)</h3></td>'+$
                      ' <td><h3>  (cm<sup>-1</sup>) </h3></td>'   
              printf, lup, '</tr> '


;      FOR i=0, n_elements(lvl1) -1 DO BEGIN 

              FOR i=i_start, i_END DO BEGIN 

                 l1 = lvl1(i) & l2=lvl2(i)

                 itrans= where(lvl1 EQ l1 AND lvl2 EQ l2, nn)


                 IF nn EQ 1 THEN BEGIN 
                    lambda1 = wvl1(i) & lambda1=lambda1(0)
                    IF  lambda1  LT 0 THEN BEGIN
                       e1 = ecmth[l1-1] &  e2=ecmth[l2-1]
                       wave = '*'+string(abs(wvl1(i)), '(f15.4)')
                    ENDIF ELSE BEGIN 
                       e1 = ecm[l1-1] &  e2=ecm[l2-1]
                       wave = string(wvl1(i), '(f15.4)')
                    END 

; correct for problems: 
                    if l1 ne 1 and e1 eq 0 then begin 
                       print,' problem in wavelength'
                       e1=ecmth[l1-1] 
                    end 
                    if l2 ne 1 and e2 eq 0 then begin 
                       print,' problem in wavelength'
                       e2=ecmth[l2-1] 
                    end 



                    printf, lup, '<tr>'


                    text='<td> '+trim(wave)+' </td>'+$
                         '<td> '+string(gf1(i), '(e10.2)')+ ' </td>'+$
                         '<td> '+ string(a_value1(i), '(e10.2)')+' </td>'

                    for l=0,nsel-1 do text=text+$
                       '<td> '+ string(ratios[l,i], '(e10.2)')+' </td>'

                    text=text+ '<td> '+trim(config(l1-1)) +' </td>'+$
                         '<td> <sup>'+ trim(s_values(l1-1))+'</sup>'+trim(spd_values(l1-1))+$
                         '<sub>'+trim(j_str(l1-1))+' </sub></td>'+$
                         '<td> '+trim(config(l2-1)) +' </td>'+$
                         ' <td><sup>'+trim(s_values(l2-1))+'</sup>'+trim(spd_values(l2-1))+$
                         '<sub>'+trim(j_str(l2-1))+' </sub>' +' </td>'+$
                         '<td> '+string(l1, '(i3)')+' </td>'+$
                         '<td> '+  string(l2, '(i3)') +' </td>'+$
                         '<td> '+ trim(e1)+'</td> '+$
                         '<td> '+trim(e2)+' </td>'

                    printf, lup, text

                    printf, lup, '</tr> '
                 END 
              ENDFOR 
              printf, lup, '</table><p> '

              printf, lup, '<em> Created by Giulio Del Zanna on '+ systime() +'</em>'

              printf, lup, '</body></html>'

              free_lun, lup

           ENDFOR 

           if n_tables gt  1 then printf, luo, ' <p> </ul>'

           printf, luo, '<p> '
           printf, luo, '<em> Page created by Giulio Del Zanna on '+ systime() +'</em>'
           printf, luo, '</body></html>'

           free_lun, luo
           
        ENDIF    ELSE BEGIN 

           print, 'Warning  with  '+directory+'/'+dirs_elements(el)+'/'+ion_string
           print, 'MISSING complete set of file ! '

        END 
        

     ENDFOR                     ;each ion 

;filename = findfile(concat_dir(dirs(k),'*.pro'),   count = cnt)
;IF cnt GT 0 THEN
;break_file, files(i), disk_log, dir, file, ext

finish:

;stop

  ENDFOR                        ;each element


END 
