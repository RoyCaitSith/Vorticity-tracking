load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
    
    file_dir      = "/raid63/cffeng/14_vorticity_track_hodges_all_year"
    file_vor      = file_dir+"/pf_smooth_data"
    file_pic      = file_dir+"/picture"
    file_tra      = file_dir+"/best_track"

    file_vor_cho       = new(2,"string")
    file_tra_cho       = new(1,"string")

    res = True
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = True
    res@gsnAddCyclic         = False
    res@cnFillOn             = True
    res@cnLinesOn            = False
;    res@lbLabelBarOn         = False
    res@lbLabelBarOn         = True
    res@cnLevelSelectionMode = "ExplicitLevels"
;    res@cnLevels             = (/0,1,2,3,5,7,9,11,13,15,17,19,21,23,25/)*10^-6
;    res@cnLevels             = (/5,10,15,20,25,30,35,40,45,50/)*10^-6
    res@mpCenterLonF         = 180
;    res@mpMinLonF            = 60
;    res@mpMaxLonF            = 230
;    res@mpMinLatF            = -40
;    res@mpMaxLatF            = 40
    res@mpMinLonF            = 60
    res@mpMaxLonF            = 230
    res@mpMinLatF            = 0
    res@mpMaxLatF            = 40

    resP = True
    resP@gsnMaximize         = True
;    resP@gsnPanelLabelBar    = True 
;    resP@lbOrientation       = "vertical"

    pres = True
    pres@gsLineThicknessF    = 2.0

    mres = True
    mres@gsMarkerIndex       = 16
    mres@gsMarkerSizeF       = 4.0
    mres@gsMarkerColor       = "black"

    res_now = True
    res_now@gsMarkerIndex    = 16
    res_now@gsMarkerSizeF    = 4.0
    res_now@gsMarkerColor    = "green"

    vcres = True
    vcres@gsnDraw            = False
    vcres@gsnFrame           = False
    vcres@gsnMaximize        = True
    vcres@vcRefMagnitudeF    = 10.0
    vcres@vcRefLengthF       = 0.025
    vcres@vcGlyphStyle       = "CurlyVector"
    vcres@gsnAddCyclic       = False
;    vcres@mpCenterLonF       = 180
;    vcres@mpMinLonF          = 60
;    vcres@mpMaxLonF          = 230
;    vcres@mpMinLatF          = 0
;    vcres@mpMaxLatF          = 40

    plots = new(2,"graphic")

    do i_year = 1979,1979

        file_vor_cho(0) = file_vor+"/vor_wnp_1p000_19000101-19001231.nc"
        file_vor_cho(1) = file_vor+"/vor_wnp_1p000_npf_2p00_gaus_19000101-19001231.nc" 
        file_tra_cho(0) = file_tra+"/ff_trs_pos_low_select_1900.new"
;        file_tra_cho(0) = file_tra+"/ff_trs_pos_low_1900.new"
        file_vor_cho(0) = str_sub_str(file_vor_cho(0),"1900",tostring(i_year))
        file_vor_cho(1) = str_sub_str(file_vor_cho(1),"1900",tostring(i_year))
        file_tra_cho(0) = str_sub_str(file_tra_cho(0),"1900",tostring(i_year))

        system("mkdir "+file_pic+"/"+tostring(i_year)+"/")

        a = addfile(file_vor_cho(0),"r") 
        low_1p000_raw = a->vo
;        printVarSummary(low_1p000_raw)
        delete(a)

        a = addfile(file_vor_cho(1),"r") 
        low_1p000_npf_spectral_filter = a->vo
;        printVarSummary(low_1p000_npf_spectral_filter)
        delete(a)

        time = low_1p000_raw&time
        n_time = dimsizes(time)

        nrow_low_1p000_spectral_filter_fnl = numAsciiRow(file_tra_cho(0))
        besttrack_low_1p000_spectral_filter_fnl = asciiread(file_tra_cho(0),nrow_low_1p000_spectral_filter_fnl,"string")

        do i = 0,dimsizes(time)-1
;        do i = 0,100

            time_now = cd_calendar(time(i),3)
            print(time_now)
            wks = gsn_open_wks("png",file_pic+"/"+tostring(i_year)+"/"+time_now)

            gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")

            res@cnLevelSelectionMode = "ExplicitLevels"
            res@cnLevels             = (/5,10,15,20,25,50,100,150,200/)*10^-6
            plots(0) = gsn_csm_contour_map(wks,low_1p000_raw(i,:,:),res)
            plots(1) = gsn_csm_contour_map(wks,low_1p000_npf_spectral_filter(i,:,:),res)

            dums_low_1p000_spectral_filter_fnl = new(600,"graphic")

            i_draw = 0
;            count = 4
            count = 1
            do while (count.lt.nrow_low_1p000_spectral_filter_fnl)
                lines = stringtointeger(str_get_field(besttrack_low_1p000_spectral_filter_fnl(count),2," "))
                count_next = count+lines+2
                flag = False
                do j = count+1,count_next-2
                    if(stringtointeger(str_get_field(besttrack_low_1p000_spectral_filter_fnl(j),1," ")).eq.time_now)
                        flag = True
                        find_j = j-count-1
                    end if
                end do
                if(flag)
                    xpt = new(lines,"double")
                    ypt = new(lines,"double")
                    do j = count+1, count_next-2
                        xpt(j-count-1) = stringtodouble(str_get_field(besttrack_low_1p000_spectral_filter_fnl(j),2," "))
                        ypt(j-count-1) = stringtodouble(str_get_field(besttrack_low_1p000_spectral_filter_fnl(j),3," "))
                    end do
                    dums_low_1p000_spectral_filter_fnl(i_draw) = gsn_add_polyline(wks,plots(1),xpt,ypt,pres)
                    i_draw = i_draw+1
                    dums_low_1p000_spectral_filter_fnl(i_draw) = gsn_add_polymarker(wks,plots(1),xpt,ypt,mres)
                    i_draw = i_draw+1
                    dums_low_1p000_spectral_filter_fnl(i_draw) = gsn_add_polymarker(wks,plots(1),xpt(find_j),ypt(find_j),res_now)
                    i_draw = i_draw+1
                    delete(xpt)
                    delete(ypt)
                end if
                count = count_next
            end do

            gsn_panel(wks,plots,(/2,1/),resP)
            delete(dums_low_1p000_spectral_filter_fnl)

        end do

        delete([/low_1p000_raw,low_1p000_npf_spectral_filter/])
        delete([/time,n_time/])
        delete([/nrow_low_1p000_spectral_filter_fnl,besttrack_low_1p000_spectral_filter_fnl/])

    end do

end
