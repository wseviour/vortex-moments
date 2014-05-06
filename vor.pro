FUNCTION VOR,field,lon,lat,GP=gp,LEVEL=level,SH=sh,NOSHOW=noshow,$
  EDGE=edge,SPLIT=split,MARS=mars,ELLIPSE=ellipse

;+
; NAME : vor
; DESCRIPTION : calculate elliptical diagnostics of the stratospheric polar 
;               vortex. Default assumes the field is NH and potential vorticity 
;
; METHOD : Uses the method laid out in Mitchell et al., JAS, 2011 to calculate
;          moments of a longitude x latitude field.
;
; Current Code Owner: Dann Mitchell (mitchell@atm.ox.ac.uk)
;
; Citation: If you use this code please cite the following papers:
;           
;           Mitchell, D.M., Charlton-Perez, A.J. and Gray, L.J, 2011. Characterising 
;           the Variability and Extremes of the Stratospheric Polar Vortices 
;           Using 2D Moments, J. Atmos. Sci., 1194-1213.
;
;           Mitchell, Daniel M., Lesley J. Gray, James Anstey, Mark P. Baldwin, 
;           Andrew J. Charlton-Perez, 2013: The Influence of Stratospheric Vortex 
;           Displacements and Splits on Surface Climate. J. Climate, 26, 2668–2682. 
;
; History:
; Version Date     Comment
; ------- -------- -------
;  0.1    4/9/2008   DMM set up code
;  1.0    23/9/2008  DMM initial version
;  1.1    6/7/2012   DMM lots of iterations lost in the mists of time
;  1.2    6/7/2012   DMM script made into a user friendly IDL function
;  1.3    7/7/2012   DMM added the ability to look at split vortices
;  1.4    14/8/2012  DMM added ability to output lons/lats for ellipse
;  1.5    14/8/2012  DMM added a keywork to make code work for Mars
;  1.6    21/8/2012  DMM fixed issue issue with ELLIPSE keyword
;  1.7    3/12/2012  DMM corrected issues with the area diagnostic, also fed
;                        to ellipse diagnostic. Added in debug keyword.
;  1.8    3/1/2012   DMM,CZ Correct error in SH diagnostics, used the correct
;                        polar steriographic tansform (Lambert) (I think
;                        Waugh, 99 have this wrong, need to chat with Darryn
;                        about this).
;
; PURPOSE : To calculate elliptical diagnostics
;
; CALLING SEQUENCE : result = VOR(field,lon,lat,EDGE=edge,/SPLIT,/SH $
;                    ,/GP,LEVEL=level,/MARS,/ELLIPSE,/NOSHOW)
;
; INPUTS : field - a 2D field of dimensions longitude latitude. By default the
;                  field is assumed to be potential vorticity, but use the 
;                  keyword /GP if you want to use geopotential
;          lon - the corresponding longitudes in degrees
;          lat - the corresponding latitudes in degrees
;
; OPTIONAL INPUTS : edge - if edge is specified then code will use this value
;                          to define the vortex region, otherwise it will use 
;                          standard values for each level. Note that either 
;                          EDGE or THETA must be specified.
;                   split - if the vortex splits into 2 calculate 
;                           diagnostics for both - this slows the code
;                   sh - calculate the diagnostics for the SH (default is NH)
;                   gp - the field is geopotential (default is PV)
;                   level - height of 2D surface. For PV this must be on 
;                           theta surfaces, for GP this must be on pressure 
;                           levels. Be sure to make sure you read in fields 
;                           with standard units. Only options are 850K (for PV)
;                           or 10hPa (for geopotential). NOTE THAT SPECIFYING 
;                           EDGE IS FAR MORE EFFECTIVE.
;                   noshow - do now output summary to screen
;                   mars - to calculate diagnostics on Martian polar vortices.
;                   ellipse - also outputs an array of 150 points for longitude 
;                             and latitude of ellipse. Array size= [150,2], where 
;                             the (*,0)=longitudes and (*,1)=latitudes.
;                   debug - purely for debugging the code if it is used as function.
;
; OUTPUTS : area - area of the vortex in km^2
;           obj_area - objective area of the vortex
;           latcent - latitude of the centre of the vortex
;           loncent - longitude of the centre of the vortex
;           theta - angle of the major axis of the vortex
;           ar - aspect ratio of the vortex
;           kurtosis - kurtosis of the vortex
;
;           **note if the names have a 1 or 2 after they are for each 
;           individual vortex for a split. This will only happen if the
;           /SPLIT keyword is used**
;
; RETURN VALUE : Data structure containing elliptical diagnostics.
;
; EXAMPLES : 
;
; ;minimal calling sequence, calculate diagnostics of the 850K PV field
; ;note it is normally better to specify EDGE rather than LEVEL
; d=VOR(field,lon,lat,LEVEL=850)
; help, d
;
;
; ;alternate minimal calling sequence, specifying edge value yourself 
; ; (this is safest option)
; d=VOR(field,lon,lat,EDGE=0.00046)
; help,max(field),0.00046
; help,d
;
;
; ;calculate diagnostics for the SH using geopotential at 10hPa and not 
; ;caring about splits vortices. 
;
; d=VOR(field,lon,lat,EDGE=2.9e5,/SH,/GP,/SPLIT)
; print,d.ar,d.ar1,d.ar2
; 
;-

; PRINT, '>>> eldi_vor: .begin...'

; >>> 17 Dec 2012 CZ >>>
; ; --- test ---
;
; HELP, field, lon, lat
; field( *, * ) = edge
; take = WHERE( lat GE 60. )
; field( *, take ) = 2. * edge
; <<< 17 Dec 2012 CZ <<<


; ----- preparations ---

longi=lon
latit=lat


; --- edge defaults ---

; If user has used LEVEL option rather than defining EDGE then standard EDGE values for that value
; are used, these can be found in Waugh et al, 1999 and Mitchell et al, 2012

IF (KEYWORD_SET(level) and KEYWORD_SET(edge)) THEN BEGIN
    print,'>>> eldi_vor: You can not specify both LEVEL and EDGE'
    stop
ENDIF

IF (KEYWORD_SET(level)) THEN BEGIN
                                ; For now only one level has been defined 850K or 10hPa
    IF (KEYWORD_SET(gp)) THEN edge=30400.0
    IF (NOT KEYWORD_SET(gp)) THEN edge=0.0005
ENDIF
; PRINT, '>>> eldi_vor: edge = ', edge


; --- show field/level information ---

;Display to user what they will be calculating, unless /noshow keyword is used

IF (NOT KEYWORD_SET(noshow)) THEN BEGIN

    IF (NOT KEYWORD_SET(gp)) THEN BEGIN
        print,'>>> eldi_vor: FIELD = POTENTIAL VORTICITY'
        IF (NOT KEYWORD_SET(level)) THEN print,'>>> eldi_vor: LEVEL = UNKNOWN'
        IF (KEYWORD_SET(level)) THEN print,'>>> eldi_vor: LEVEL = '+STRTRIM(level, 2)+'K'
        IF (NOT KEYWORD_SET(sh)) THEN print,'>>> eldi_vor: HEMISPHERE = NORTHERN'
        IF (KEYWORD_SET(sh)) THEN print,'>>> eldi_vor: HEMISPHERE = SOUTHERN'
    ENDIF

    IF (KEYWORD_SET(gp)) THEN BEGIN
        print,'>>> eldi_vor: FIELD = GEOPOTENTIAL'
        IF (NOT KEYWORD_SET(level)) THEN print,'>>> eldi_vor: LEVEL = UNKNOWN'
        IF (KEYWORD_SET(level)) THEN print,'>>> eldi_vor: LEVEL = '+STRTRIM(level, 2)+ 'hPa'
        IF (NOT KEYWORD_SET(sh)) THEN print,'>>> eldi_vor: HEMISPHERE = NORTHERN'
        IF (KEYWORD_SET(sh)) THEN print,'>>> eldi_vor: HEMISPHERE = SOUTHERN'
    ENDIF

ENDIF

; Define hemisphere
IF (NOT KEYWORD_SET(sh)) THEN latpos=WHERE(latit gt 0)
IF (KEYWORD_SET(sh)) THEN latpos=WHERE(latit lt 0)

; Conversions
latit=latit(latpos)/!radeg
longi=longi/!radeg

;Isolate only the correct hemisphere
data=field(*,latpos)


; - Calculate Vortex Area -

; Earth radius [ meters ]
a = 6374.0e3
IF (KEYWORD_SET(mars)) THEN a = 3397.0e3
; >>> 06 Dec 2012 CZ >>>
; save as "radius" because "a" becomes an index, later
radius = a
; <<< 06 Dec 2012 CZ <<<

;the length in km of each latitude change
; >>> 17 Dec 2012 CZ >>>
; latlen=111.2*ABS(latit(0)*!radeg-latit(1)*!radeg)
latlen=2*!PI*a/360.e3*ABS(latit(0)*!radeg-latit(1)*!radeg)
; <<< 17 Dec 2012 CZ <<<

;the length in km of each longitude change
diam=(2*!PI*a*COS(latit)/n_elements(longi))/1000.

;the area of each sector
sector=latlen*diam

world=fltarr(n_elements(longi),n_elements(latit))
mask=fltarr(n_elements(longi),n_elements(latit))
mask(*,*)=0

FOR count=0,n_elements(longi)-1  DO world(count,*)=sector


; - define vortex -

; define vortex indices
IF (NOT KEYWORD_SET(gp)) THEN vortex=WHERE(data gt edge)
IF (KEYWORD_SET(gp)) THEN vortex=WHERE(data lt edge)

; mask vortex
mask(vortex)=1

; vortex area
varea=mask*world
area=TOTAL(varea)


; - convert from spherical polar to cartesian -

; >>> 06 Dec 2012 CZ >>> done above!
; Earth radius
; a = 6374.0e3
;  IF (KEYWORD_SET(mars)) THEN a = 3397.0e3
; <<< 06 Dec 2012 CZ <<<

;ydir=findgen(9200)

; >>> 17 Dec 2012 CZ >>>
; R = SQRT(2*(1-sin(latit)))
R = SQRT(2*(1-abs(sin(latit))))
; <<< 17 Dec 2012 CZ <<<

X=FLTARR(n_elements(longi),n_elements(latit))
Y=FLTARR(n_elements(longi),n_elements(latit))

FOR lo=0,n_elements(longi)-1 DO BEGIN
    FOR la=0,n_elements(latit)-1 DO BEGIN
; >>> 17 Dec 2012 CZ >>>
; old:
;                                 ;Conversion for the NH
;         IF (NOT KEYWORD_SET(sh)) THEN BEGIN
;             x(lo,la)=((cos(longi(lo)))*(cos(latit(la))))/(1+(sin(latit(la))))
;             y(lo,la)=(((sin(longi(lo)))*(cos(latit(la)))))/(1+(sin(latit(la))))
;         ENDIF
;                                 ;Conversion for the SH
;         IF (KEYWORD_SET(sh)) THEN BEGIN
;             x(lo,la)=((cos(longi(lo)))*(cos(latit(la))))/(1-(sin(latit(la))))
;             y(lo,la)=(-1)*(((sin(longi(lo)))*(cos(latit(la)))))/(1-(sin(latit(la))))
;         ENDIF
; new:
                                ;Conversion for the NH
        IF (NOT KEYWORD_SET(sh)) THEN BEGIN
            x(lo,la)=R(la)*cos(longi(lo))
            y(lo,la)=R(la)*sin(longi(lo))
        ENDIF
                                ;Conversion for the SH
        IF (KEYWORD_SET(sh)) THEN BEGIN
            x(lo,la)=R(la)*cos(longi(lo))
            y(lo,la)=-R(la)*sin(longi(lo))
        ENDIF
; <<< 17 Dec 2012 CZ <<<
    ENDFOR
ENDFOR

triangulate,x,y,tr

xgrid=-1+findgen(n_elements(longi))/(.5*n_elements(longi))
ygrid=-1+findgen(n_elements(longi))/(.5*n_elements(longi))

;This is now projected onto a cartesian grid
PVXY=trigrid(x,y,data,tr,xout=xgrid,yout=ygrid)

; Replace PV into two fields, an area = to pvsurf, an area above pvsurf
; The vortex in PV is greater than surrounding pv, but the vortex is GP is less than surrounding PV


FOR a=0,n_elements(longi)-1 DO BEGIN
    FOR b=0,n_elements(longi)-1 DO BEGIN
        
        IF (KEYWORD_SET(gp)) THEN BEGIN
            
            IF (PVXY(a,b) GT edge) THEN PVXY(a,b)=edge
            
        ENDIF
        
        IF (NOT KEYWORD_SET(gp)) THEN BEGIN
            
            IF (NOT KEYWORD_SET(sh)) THEN BEGIN
                IF (PVXY(a,b) LT edge) THEN PVXY(a,b)=edge
            ENDIF
        
            IF (KEYWORD_SET(sh)) THEN BEGIN
                IF (PVXY(a,b) GT edge) THEN PVXY(a,b)=edge
            ENDIF
            
        ENDIF    

    ENDFOR
ENDFOR



;Ensure low latitudes are zero for GP
IF (KEYWORD_SET(gp)) THEN BEGIN

    bad=WHERE(pvxy eq 0)
    pvxy(bad)=edge

ENDIF

; >>> 06 Dec 2012 CZ >>>
; ;Setup arrays to fit ellipse
;
; rellip=fltarr(150)
; thiel=fltarr(150)
; lamel=fltarr(150)
;
; rellipR1=fltarr(150)
; thielR1=fltarr(150)
; lamelR1=fltarr(150)
;
; rellipR2=fltarr(150)
; thielR2=fltarr(150)
; lamelR2=fltarr(150)
; <<< 06 Dec 2012 CZ <<<

;========================== Calculate the moments for a displacement========================

; Find the area of each box in cartesian coord system

; >>> 06 Dec 2012 CZ >>>
; radius=6371.0e3
; <<< 06 Dec 2012 CZ <<<
blength=(2*radius)/n_elements(longi)
barea=blength^2
; >>> 17 Dec 2012 CZ >>>
;;; bfact = !PI / 1.5e6
bfact = 1. / 1.E6
efact = 1. / 1.E6
; <<< 17 Dec 2012 CZ <<<

; >>> 06 Dec 2012 CZ >>>
; old:
; ; Calculate moments
;
; M00=0
; M10=0
; M01=0
; M00area=0
;
; FOR a=0,n_elements(longi)-1 DO BEGIN
;     FOR b=0,n_elements(longi)-1 DO BEGIN
;
;         M1 = ABS(pvxy(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^0)
;         M00 = M00+(M1)
;         M2 = ABS(pvxy(a,b)-edge)*(xgrid(a)^1)*(ygrid(b)^0)
;         M10 = M10+(M2)
;         M3 = ABS(pvxy(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^1)
;         M01 = M01+M3
;         MAREA = ABS(pvxy(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^0)*barea
;         M00area=M00area+MAREA
;
;     ENDFOR
; ENDFOR
;
; ; Calculate the centroid
;
; centx=M10/M00
; centy=M01/M00
;
; ;convert back to polar coords
;
; R=centx^2+centy^2
; loncent=90.-atan(centx,centy)*!RADEG ;centroid longitude
; latcent=asin((1-R)/(1+R))*!RADEG ;centroid latitude
;
;
; ;=============== Calculate the relative moments of the vortex displacement========================
;
;  J11=0
;  J20=0
;  J02=0
; J22=0
; J40=0
; J04=0
;
; FOR a=0, n_elements(longi)-1 DO BEGIN
;     FOR b=0,n_elements(longi)-1 DO BEGIN
;
;         J1 = ABS(pvxy(a,b)-edge)*((xgrid(a)-centx)^1)*((ygrid(b)-centy)^1)
;         J11 = J11+J1
;         J2 = ABS(pvxy(a,b)-edge)*((xgrid(a)-centx)^2)*((ygrid(b)-centy)^0)
;         J20 = J20+J2
;         J3 = ABS(pvxy(a,b)-edge)*((xgrid(a)-centx)^0)*((ygrid(b)-centy)^2)
;         J02 = J02+J3
;         J4 = ABS(pvxy(a,b)-edge)*((xgrid(a)-centx)^2)*((ygrid(b)-centy)^2)
;         J22 = J22+J4
;         J5 = ABS(pvxy(a,b)-edge)*((xgrid(a)-centx)^4)*((ygrid(b)-centy)^0)
;         J40 = J40+J5
;         J6 = ABS(pvxy(a,b)-edge)*((xgrid(a)-centx)^0)*((ygrid(b)-centy)^4)
;         J04 = J04+J6
;
;     ENDFOR
; ENDFOR
;
; ; Angle between the x-axis and major axis of ellipse
; ; note that atan is calculated so the ration (2*J11)/(J20-J02) is considered
;
; ang = 0.5*atan((2*J11),(J20-J02))
;
; ; Aspect ratio
;
; AR = SQRT(abs(((J20+J02)+SQRT(4*(J11^2)+(J20-J02)^2))/((J20+J02)-SQRT(4*(J11^2)+(J20-J02)^2))))
;
; ; Equivalent area
;
; EA = M00area/(edge)
;
;
; ; major/minor axes
;
; minor=SQRT((area*2e5)/(!pi*AR))*2
; major=AR*minor
;
; ; Excess kurtosis
;
; KUR=M00*((J40+(2.0*J22)+J04)/(J20+J02)^2.0) - (2.0/3.0)*(((3.0*(ar^4.0))+(2.0*(ar^2.0))+3.0) $
;                                                          /((ar^2.0)+1.0)^2.0)
;
; ;str={area:area, obj_area:ea, ar:ar, latcent:latcent, loncent:loncent, theta:ang, kurtosis:kur}
; new:
; PRINT, '>>> eldi_vor: .str0 = CALC_VOR(PVXY)...'
; >>> 17 Dec 2012 CZ >>>
; str0 = CALC_VOR( longi, pvxy, edge, xgrid, ygrid, radius, area, barea )

; transfer sh key to new function

IF (KEYWORD_SET(sh)) THEN str0 = CALC_VOR( longi, pvxy, edge, xgrid, ygrid, radius, barea, /SH)

IF NOT (KEYWORD_SET(sh)) THEN str0 = CALC_VOR( longi, pvxy, edge, xgrid, ygrid, radius, barea)

; PRINT, '>>> vor:   area old = ', area, ' vs new ', str0.area * bfact
; STOP
area = str0.area * bfact
; ea = str0.obj_area
ea = str0.obj_area * efact
; <<< 17 Dec 2012 CZ <<<
ar = str0.ar
latcent = str0.latcent
loncent = str0.loncent
ang = str0.theta
kur = str0.kurtosis
ell = str0.ellipse
centx = str0.centx
centy = str0.centy
; STOP
; <<< 06 Dec 2012 <<<


IF (KEYWORD_SET(split)) THEN BEGIN


; ----- Characterising warmings into wave 1/wave 2 -----

; Define a splitting event to be when the excess kurtosis is less than -0.1

    R1=fltarr(n_elements(longi),n_elements(longi))
    R2=fltarr(n_elements(longi),n_elements(longi))


; --- SPLITTING SECTION ---

    IF (kur LT -0.1) THEN BEGIN

        FOR xp=0,n_elements(longi)-1 DO BEGIN
            FOR yp=0,n_elements(longi)-1 DO BEGIN

                IF (ygrid(yp) GT ((-1.0/tan(ang))*(xgrid(xp)-centx)+centy)) THEN BEGIN
                    R1(xp,yp)=0
                ENDIF ELSE BEGIN
                    R1(xp,yp)=pvxy(xp,yp)

                ENDELSE

                IF (ygrid(yp) LT ((-1.0/tan(ang))*(xgrid(xp)-centx)+centy)) THEN BEGIN
                    R2(xp,yp)=0
                ENDIF ELSE BEGIN
                    R2(xp,yp)=pvxy(xp,yp)

                ENDELSE
            ENDFOR
        ENDFOR


; - Make everything outside the vortex region equal to the edge value -

        FOR a=0,n_elements(longi)-1 DO BEGIN
            FOR b=0,n_elements(longi)-1 DO BEGIN

                IF (KEYWORD_SET(gp)) THEN BEGIN
                    IF (R1(a,b) GT edge) THEN R1(a,b)=edge
                    IF (R2(a,b) GT edge) THEN R2(a,b)=edge
                ENDIF

                IF (NOT KEYWORD_SET(gp)) THEN BEGIN
                    IF (R1(a,b) LT edge) THEN R1(a,b)=edge
                    IF (R2(a,b) LT edge) THEN R2(a,b)=edge
                ENDIF

            ENDFOR	; b
        ENDFOR ; a


; - Ensure low latitudes are zero for GP -

        IF (KEYWORD_SET(gp)) THEN BEGIN

            bad=WHERE(R1 eq 0)
            R1(bad)=edge

            bad=WHERE(R2 eq 0)
            R2(bad)=edge

        ENDIF ; KEYWORD_SET(gp)

; >>> 06 Dec 2012 CZ >>>
; old:
; ; Calculate moments for each of the split vortices
;
;         M00R1=0
;         M10R1=0
;         M01R1=0
;         M00areaR1=0
;         M00R2=0
;         M10R2=0
;         M01R2=0
;         M00areaR2=0
;
;         FOR a=0,n_elements(longi)-1 DO BEGIN
;             FOR b=0,n_elements(longi)-1 DO BEGIN
;
;                 M1 = ABS(R1(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^0)
;                 M00R1 = M00R1+M1
;                 M2 = ABS(R1(a,b)-edge)*(xgrid(a)^1)*(ygrid(b)^0)
;                 M10R1 = M10R1+M2
;                 M3 = ABS(R1(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^1)
;                 M01R1 = M01R1+M3
;                 MAREA = ABS(R1(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^0)*barea
;                 M00areaR1=M00areaR1+MAREA
;                 M1 = ABS(R2(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^0)
;                 M00R2 = M00R2+M1
;                 M2 = ABS(R2(a,b)-edge)*(xgrid(a)^1)*(ygrid(b)^0)
;                 M10R2 = M10R2+M2
;                 M3 = ABS(R2(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^1)
;                 M01R2 = M01R2+M3
;                 MAREA = ABS(R2(a,b)-edge)*(xgrid(a)^0)*(ygrid(b)^0)*barea
;                 M00areaR2=M00areaR2+MAREA
;
;             ENDFOR
;         ENDFOR
;
; ; Calculate the centroid
;
;         centxR1=M10R1/M00R1
;         centyR1=M01R1/M00R1
;         centroidR1=[centxR1,centyR1]
;         centxR2=M10R2/M00R2
;         centyR2=M01R2/M00R2
;         centroidR2=[centxR2,centyR2]
;
; ;=============== Calculate the relative moments of the vortex split========================
;
;         J11R1=0
;         J20R1=0
;         J02R1=0
;         J22R1=0
;         J40R1=0
;         J04R1=0
;         J11R2=0
;         J20R2=0
;         J02R2=0
;         J22R2=0
;         J40R2=0
;         J04R2=0
;
;         FOR a=0,n_elements(longi)-1 DO BEGIN
;             FOR b=0,n_elements(longi)-1 DO BEGIN
;
;                 J1 = ABS(R1(a,b)-edge)*((xgrid(a)-centxR1)^1)*((ygrid(b)-centyR1)^1)
;                 J11R1 = J11R1+J1
;                 J2 = ABS(R1(a,b)-edge)*((xgrid(a)-centxR1)^2)*((ygrid(b)-centyR1)^0)
;                 J20R1 = J20R1+J2
;                 J3 = ABS(R1(a,b)-edge)*((xgrid(a)-centxR1)^0)*((ygrid(b)-centyR1)^2)
;                 J02R1 = J02R1+J3
;                 J4 = ABS(R1(a,b)-edge)*((xgrid(a)-centxR1)^2)*((ygrid(b)-centyR1)^2)
;                 J22R1 = J22R1+J4
;                 J5 = ABS(R1(a,b)-edge)*((xgrid(a)-centxR1)^4)*((ygrid(b)-centyR1)^0)
;                 J40R1 = J40R1+J5
;                 J6 = ABS(R1(a,b)-edge)*((xgrid(a)-centxR1)^0)*((ygrid(b)-centyR1)^4)
;                 J04R1 = J04R1+J6
;                 J1 = ABS(R2(a,b)-edge)*((xgrid(a)-centxR2)^1)*((ygrid(b)-centyR2)^1)
;                 J11R2 = J11R2+J1
;                 J2 = ABS(R2(a,b)-edge)*((xgrid(a)-centxR2)^2)*((ygrid(b)-centyR2)^0)
;                 J20R2 = J20R2+J2
;                 J3 = ABS(R2(a,b)-edge)*((xgrid(a)-centxR2)^0)*((ygrid(b)-centyR2)^2)
;                 J02R2 = J02R2+J3
;                 J4 = ABS(R2(a,b)-edge)*((xgrid(a)-centxR2)^2)*((ygrid(b)-centyR2)^2)
;                 J22R2 = J22R2+J4
;                 J5 = ABS(R2(a,b)-edge)*((xgrid(a)-centxR2)^4)*((ygrid(b)-centyR2)^0)
;                 J40R2 = J40R2+J5
;                 J6 = ABS(R2(a,b)-edge)*((xgrid(a)-centxR2)^0)*((ygrid(b)-centyR2)^4)
;                 J04R2 = J04R2+J6
;
;             ENDFOR
;         ENDFOR
;
; ; Angle between the x-axis and major axis of ellipse
; ; note that atan is calculated so the ration (2*J11)/(J20-J02) is considered
;
;         angR1 = 0.5*atan((2*J11R1),(J20R1-J02R1))
;         angR2 = 0.5*atan((2*J11R2),(J20R2-J02R2))
;
; ; Aspect ratio
;
;         ARR1 = SQRT(abs(((J20R1+J02R1)+SQRT(4*(J11R1^2)+$
;                                             (J20R1-J02R1)^2))/((J20R1+J02R1)-SQRT(4*(J11R1^2)+(J20R1-J02R1)^2))))
;         ARR2 = SQRT(abs(((J20R2+J02R2)+SQRT(4*(J11R2^2)+$
;                                             (J20R2-J02R2)^2))/((J20R2+J02R2)-SQRT(4*(J11R2^2)+(J20R2-J02R2)^2))))
;
; ; Equivalent area
;
;         EAR1 = M00areaR1/edge
;         EAR2 = M00areaR2/edge
;
; ; major/minor axes
;
;         minorR1=SQRT(EAR1/(!pi*ARR1))*2
;         majorR1=ARR1*minorR1
;
;         axesR1=[majorR1,minorR1]
;
;         minorR2=SQRT(EAR2/(!pi*ARR2))*2
;         majorR2=ARR2*minorR2
;
;         axesR2=[majorR2,minorR2]
;
; ;==========================Plotting the ellipse for splitting events===========================
;
;         points=150
;
; ; Divide a circle into sections
;
;         phi = 2 * !PI * (Findgen(points)/(points-1))
;
;                                 ; Parameterized equation of ellipse.
;
;         xR1 =  (majorR1/2)*Cos(phi)
;         yR1 =  (minorR1/2)*Sin(phi)
;         xR2 =  (majorR2/2)*Cos(phi)
;         yR2 =  (minorR2/2)*Sin(phi)
;
;                                 ; Rotate to desired position angle.
;
;         xprimeR1 = (centxR1*radius) + (xR1*cos(angR1)) - (yR1*sin(angR1))
;         yprimeR1 = (centyR1*radius) + (xR1*sin(angR1)) + (yR1*cos(angR1))
;         xprimeR2 = (centxR2*radius) + (xR2*cos(angR2)) - (yR2*sin(angR2))
;         yprimeR2 = (centyR2*radius) + (xR2*sin(angR2)) + (yR2*cos(angR2))
;
;                                 ; Extract the points to return.
;
;         ptsR1 = FltArr(2, N_Elements(xprimeR1))
;
;         ptsR1[0,*] = xprimeR1
;         ptsR1[1,*] = yprimeR1
;
;         ptsR2 = FltArr(2, N_Elements(xprimeR2))
;
;         ptsR2[0,*] = xprimeR2
;         ptsR2[1,*] = yprimeR2
;
; ; convert back to polars
; ; Polar radius
;         RR1=centxR1^2+centyR1^2
;         RR2=centxR2^2+centyR2^2
;
;         FOR j=0,149 DO BEGIN
;
;             rellipR1(j)=(ptsR1[0,j]/radius)^2+(ptsR1[1,j]/radius)^2
;             thielR1(j)=asin((1-rellipR1(j))/(1+rellipR1(j)))*!radeg
;
;
;             rellipR2(j)=(ptsR2[0,j]/radius)^2+(ptsR2[1,j]/radius)^2
;             thielR2(j)=asin((1-rellipR2(j))/(1+rellipR2(j)))*!radeg
;
;
;         ENDFOR
;
; ; Centroid
; ; >>> 05 Dec 2012 CZ >>>
; ; old:
; ;         lamR1=90.-atan(centxR1,centyR1)
; ;         thiR1=asin((1-RR1)/(1+RR1))
; ;
; ;         lamR2=90.-atan(centxR2,centyR2)
; ;         thiR2=asin((1-RR2)/(1+RR2))
; ; new:
;         lonCentR1=90.-atan(centxR1,centyR1)*!RADEG
;         latCentR1=asin((1-RR1)/(1+RR1))*!RADEG
;
;         lonCentR2=90.-atan(centxR2,centyR2)*!RADEG
;         latCentR2=asin((1-RR2)/(1+RR2))*!RADEG
; ; <<< 05 Dec 2012 CZ <<<
;
; ; Convert to unit coords
;         uptsR1=ptsR1/radius
;         uptsR2=ptsR2/radius
;
; ; Ellipse
;
;         FOR j=0,149 DO BEGIN
;
; ; >>> 05 Dec 2012 CZ >>>
; ;             lamelR1(j)=90.-atan(uptsR1[0,j],uptsR1[1,j])
;             lamelR1(j)=90.-atan(uptsR1[0,j],uptsR1[1,j])*!RADEG
; ; <<< 05 Dec 2012 CZ <<<
;
;         ENDFOR
;
;         FOR j=0,149 DO BEGIN
;
; ; >>> 05 Dec 2012 CZ >>>
; ;             lamelR2(j)=90.-atan(uptsR2[0,j],uptsR2[1,j])
;             lamelR2(j)=90.-atan(uptsR2[0,j],uptsR2[1,j])*!RADEG
; ; <<< 05 Dec 2012 CZ <<<
;
;         ENDFOR
;
;         ellR1=fltarr(150,2)
;         ellR1(*,0)=lamelR1
;         ellR1(*,1)=thielR1
;
;         ellR2=fltarr(150,2)
;         ellR2(*,0)=lamelR2
;         ellR2(*,1)=thielR2
;
; new:

        ; first split vortex
        ; PRINT, '>>> eldi_vor: .str1 = CALC_VOR(R1)...'
; >>> 17 Dec 2012 CZ >>>
; 		str1 = CALC_VOR( longi, R1, edge, xgrid, ygrid, radius, area, barea )
		str1 = CALC_VOR( longi, R1, edge, xgrid, ygrid, radius, barea )
		areaR1 = str1.area * bfact
; 		eaR1 = str1.obj_area
       	eaR1 = str1.obj_area * efact
; <<< 17 Dec 2012 CZ <<<
		arR1 = str1.ar
		latcentR1 = str1.latcent
		loncentR1 = str1.loncent
		angR1 = str1.theta
		kurR1 = str1.kurtosis
		ellR1 = str1.ellipse
		centxR1 = str1.centx
		centyR1 = str1.centy

		; second split vortex
		; PRINT, '>>> eldi_vor: .str2 = CALC_VOR(R2)...'
; >>> 17 Dec 2012 CZ >>>
; 		str2 = CALC_VOR( longi, R2, edge, xgrid, ygrid, radius, area, barea )
		str2 = CALC_VOR( longi, R2, edge, xgrid, ygrid, radius, barea )
		areaR2 = str2.area * bfact
; 		eaR2 = str2.obj_area
		eaR2 = str2.obj_area * efact
; <<< 17 Dec 2012 CZ <<<
		arR2 = str2.ar
		latcentR2 = str2.latcent
		loncentR2 = str2.loncent
		angR2 = str2.theta
		kurR2 = str2.kurtosis
		ellR2 = str2.ellipse
		centxR2 = str2.centx
		centyR2 = str2.centy

        ; PRINT, '>>> eldi_vor:   ellR1 = ', ellR1( 0, 0 ), ellR1( 0, 1 )
; <<< 06 Dec 2012 CZ <<<
                                ;this ends the split section
	ENDIF ; (kur LT -0.1)

ENDIF ; KEYWORD_SET(split)


; >>> 06 Dec 2012 CZ >>>
; old:
; ; IF (kur GE -0.1) THEN BEGIN
;
;
; ;==========================Plotting the ellipse for events where kur GE -0.1 ==============================
;
;     points=150
;
; ; Divide a circle into sections
;
;     phi = 2 * !PI * (Findgen(points)/(points-1))
;
;                                 ; Parameterized equation of ellipse.
;
;     x =  (major/2)*Cos(phi)
;     y =  (minor/2)*Sin(phi)
;
;                                 ; Rotate to desired position angle.
;
;     xprime = (centx*radius) + (x*cos(ang)) - (y*sin(ang))
;     yprime = (centy*radius) + (x*sin(ang)) + (y*cos(ang))
;
;                                 ; Extract the points to return.
;
;     pts = FltArr(2, N_Elements(xprime))
;
;     pts[0,*] = xprime
;     pts[1,*] = yprime
;
; ;RETURN, pts
;
;     pi=3.1414926
;     a=6374.0e3
;
;     IF (KEYWORD_SET(mars)) THEN a = 3397.0e3
;
; ; convert back to polars
; ; Polar radius
;     R=centx^2+centy^2
;
;     FOR j=0,149 DO BEGIN
;
;         rellip(j)=(pts[0,j]/radius)^2+(pts[1,j]/radius)^2
;         thiel(j)=asin((1-rellip(j))/(1+rellip(j)))*!radeg
;
;     ENDFOR
;
; ; Convert to unit coords
;     upts=pts/radius
;
;  ; Ellipse
;
;     FOR j=0,149 DO BEGIN
;                                 ;by using 2 arguments in atan is takes account range larger than -pi to pi.
;         lamel(j)=90.-atan(upts[0,j],upts[1,j])*!RADEG
;
;     ENDFOR
;
; ;store the lons and lats of ellipse in a single array
; ell=fltarr(150,2)
; ell(*,0)=lamel
; ell(*,1)=thiel
; <<< 06 Dec 2012 <<<

; >>> 17 Dec 2012 CZ >>>
; ;set areas of individual vortices to NAN (this was done with at lat/lon grid with data > edge
; AREAR1=!Values.F_NAN
; AREAR2=!Values.F_NAN
; <<< 17 Dec 2012 CZ <<<

IF (NOT KEYWORD_SET(split)) THEN BEGIN
    str={area:area, obj_area:ea, ar:ar, latcent:latcent, loncent:loncent, theta:ang, kurtosis:kur}
ENDIF

IF (NOT KEYWORD_SET(split)) and (KEYWORD_SET(ellipse)) THEN BEGIN
    str={area:area, obj_area:ea, ar:ar, latcent:latcent, loncent:loncent, theta:ang, kurtosis:kur, $
    	ellipse:ell}
ENDIF

IF (KEYWORD_SET(split)) THEN BEGIN
    IF (kur GE -0.1) or (FINITE(kur) eq 0) THEN BEGIN
        EAR1=!Values.F_NAN
        EAR2=!Values.F_NAN
        ARR1=!Values.F_NAN
        ARR2=!Values.F_NAN
        ANGR1=!Values.F_NAN
        ANGR2=!Values.F_NAN
; >>> 05 Dec 2012 CZ >>>
; old:
;         LAMR1=!Values.F_NAN
;         LAMR2=!Values.F_NAN
;         THIR1=!Values.F_NAN
;         THIR2=!Values.F_NAN
; new:
        lonCentR1=!Values.F_NAN
        lonCentR2=!Values.F_NAN
        latCentR1=!Values.F_NAN
        latCentR2=!Values.F_NAN
; <<< 05 Dec 2012 CZ <<<
        AREAR1=!Values.F_NAN
        AREAR2=!Values.F_NAN
        ellR1=!Values.F_NAN
        ellR2=!Values.F_NAN
    ENDIF

; >>> 05 Dec 2012 CZ >>>
; old:
;     str={area:area, obj_area:ea, ar:ar, latcent:latcent, loncent:loncent, theta:ang, kurtosis:kur, ellipse:ell, $
;          area1:arear1, obj_area1:ear1, ar1:arr1, latcent1:thir1*!RADEG, loncent1:lamr1*!RADEG, theta1:angr1, ellipse1:ellR1,$
;          area2:arear2, obj_area2:ear2, ar2:arr2, latcent2:thir2*!RADEG, loncent2:lamr2*!RADEG, theta2:angr2, ellipse2:ellR2}
; new:
    str={area:area, obj_area:ea, ar:ar, latcent:latcent, loncent:loncent, theta:ang, kurtosis:kur, ellipse:ell, $
         area1:arear1, obj_area1:ear1, ar1:arr1, latcent1:latCentR1, loncent1:lonCentR1, theta1:angr1, ellipse1:ellR1,$
         area2:arear2, obj_area2:ear2, ar2:arr2, latcent2:latCentR2, loncent2:lonCentR2, theta2:angr2, ellipse2:ellR2}
; <<< 05 Dec 2012 CZ <<<

ENDIF

RETURN, str

END
