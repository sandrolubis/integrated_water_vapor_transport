; Calculate mass-weighted column integrated moisture budget
; Author: Dr. Sandro Lubis, Feb 2022 - PNNL
; Contact: slubis.geomar@gmail.com
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Load necessary libraries
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  ; Define constants
  g = 9.81  ; Gravity (m/s²)
  p_top = 100   ; Upper pressure level (input must be in Pa) 
  p_bot = 100000  ; Lower pressure level (Pa)


  ; Loop over years
  do year = 1979, 1980
    print("Processing year: " + year)

    ; Open input files for each year
    file_q = addfile("/pscratch/sd/s/slubis/ERA5_Data/SH/daily/q." + year + ".nc", "r")
    file_u = addfile("/pscratch/sd/s/slubis/ERA5_Data/U/daily/u." + year + ".nc", "r")
    file_v = addfile("/pscratch/sd/s/slubis/ERA5_Data/V/daily/v." + year + ".nc", "r")

    ; Read variables
    q = file_q->var133  ; (time, lev, lat, lon)
    u = file_u->var131     ; (time, lev, lat, lon)
    v = file_v->var132    ; (time, lev, lat, lon)
    
    printVarSummary(q)
    printVarSummary(u)
    printVarSummary(v)

    ; Get pressure levels from coordinate attribute
    p = u&plev    ; level must be in Pa

    ; Convert pressure to Pascals for integration
    ;p = p * 100.0

    ; Compute vertical weight dp/g
    dp = dpres_plevel_Wrap(p, p_bot, p_top, 0)  
    ;print(dp)
    dpm = conform_dims(dimsizes(q), dp, 1)  
    wgt = dpm / g  ; Compute weight for vertical integration


    ; Perform vertical integration using wgt_vertical_n
    IVTx = wgt_vertical_n(q*u, wgt, 1, 1)  ; Integrate over vertical axis (1)
    IVTy = wgt_vertical_n(q*v, wgt, 1, 1)  ; Integrate over vertical axis (1)

    ; Compute total IVT
    IVT = sqrt(IVTx^2 + IVTy^2)
    IVT@units="kg/m/s"
    
    copy_VarCoords(u(:,0,:,:),IVTx)
    copy_VarCoords(u(:,0,:,:),IVTy)
    copy_VarCoords(u(:,0,:,:),IVT)
    
    printVarSummary(IVTx)
    printVarSummary(IVTy)
    printVarSummary(IVT)

    ; Define output file for each year
    outFile = addfile("IVT_" + year + ".nc", "c")


    ; Write data to output file
    outFile->IVT = IVT
    outFile->IVTx = IVTx
    outFile->IVTy = IVTy

    print("IVT calculation complete for year: " + year + ". Output saved in IVT_" + year + ".nc")

   delete(IVT)
   delete(IVTx)
   delete(IVTy)
   delete(q)
   delete(u)
   delete(v)
   delete(dpm)
   delete(wgt)
   
  end do
end

