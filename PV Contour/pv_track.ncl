;编辑软件Notepad++
;该程序是用来寻找PV的等值线
;使用该程序的时候，有什么问题请联系冯城烽，lzfchfthdl@163.com

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
    
	;设置读入文件的路径
    file_dir = "/nfs3/jfang/fung/pv_track/"
    file_name = systemfunc("ls "+file_dir+"*.grib2")
	i_nc_file = 0
	
	;设置读取的等高面的层次
    lev = 20000
	
	;设置寻找的pv的等值线的值
    contour_pv = 0.75
    

    wks_type = "png"
    wks_type@wkHeight = 1024
    wks_type@wkWidth = 1024
    wks = gsn_open_wks(wks_type,"test")
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame = True
    resxy@xyLineThicknessF = 4.0
    res = True
    res@gsnDraw = False
    res@gsnFrame = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnAddCyclic = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/0.6,0.75,1.0/)

    nc_file = file_name(i_nc_file) 
    a = addfile(nc_file,"r")
    lv = a->lv_ISBL0
    n_lv = dimsizes(lv)
    do i = 0,n_lv-1
        if(lev.eq.lv(i))
            i_lv = i
        end if
    end do

	;计算pv
    absv_save = a->ABSV_P0_L100_GLL0(i_lv,::-1,:)
    tmp_save = a->TMP_P0_L100_GLL0(i_lv-1:i_lv+1:2,::-1,:)
    lon_save = a->lon_0
    lat = a->lat_0(::-1)
    n_lon = dimsizes(lon_save)+1
    n_lat = dimsizes(lat)
    absv = new((/n_lat,n_lon/),"double")
    tmp = new((/2,n_lat,n_lon/),"double")
    lon = new((/n_lon/),"double")
    absv(:,0:n_lon-2) = absv_save
    absv(:,n_lon-1) = absv_save(:,0)
    tmp(:,:,0:n_lon-2) = tmp_save
    tmp(:,:,n_lon-1) = tmp_save(:,:,0)
    lon(0:n_lon-2) = lon_save
    lon(n_lon-1) = 360
    theta = tmp*(100000/lev)^0.286
    pv = -10.0^6*9.8*absv*(theta(1,:,:)-theta(0,:,:))/(lv(i_lv+1)-lv(i_lv-1))
    pv!0 = "lat"
    pv!1 = "lon"
    pv@lat = lat
    pv@lon = lon

    edge = new((/n_lat-1,n_lon-1/),"logical") 
    edge = False
    edge(0,:) = True
    edge(n_lat-2,:) = True
    edge(:,0) = True
    edge(:,n_lon-2) = True

    point = new((/4,n_lat-1,n_lon-1/),"double")
    finished = new((/4,n_lat-1,n_lon-1/),"logical")
    finished = True
    do i = 0,n_lat-2
        do j = 0,n_lon-2 
            if((pv(i,j)-contour_pv)*(pv(i,j+1)-contour_pv).lt.0)
                point(0,i,j) = (contour_pv+pv(i,j+1)*j-pv(i,j)*(j+1))/(pv(i,j+1)-pv(i,j)) 
                finished(0,i,j) = False
            end if
            if((pv(i,j+1)-contour_pv)*(pv(i+1,j+1)-contour_pv).lt.0)
                point(1,i,j) = (contour_pv+pv(i+1,j+1)*i-pv(i,j+1)*(i+1))/(pv(i+1,j+1)-pv(i,j+1)) 
                finished(1,i,j) = False
            end if
            if((pv(i+1,j+1)-contour_pv)*(pv(i+1,j)-contour_pv).lt.0)
                point(2,i,j) = (contour_pv*(-1)+pv(i+1,j)*(j+1)-pv(i+1,j+1)*j)/(pv(i+1,j)-pv(i+1,j+1)) 
                finished(2,i,j) = False
            end if
            if((pv(i+1,j)-contour_pv)*(pv(i,j)-contour_pv).lt.0)
                point(3,i,j) = (contour_pv*(-1)+pv(i,j)*(i+1)-pv(i+1,j)*i)/(pv(i,j)-pv(i+1,j)) 
                finished(3,i,j) = False
            end if
        end do
    end do

    seq = new((/2,(n_lat-1)*(n_lon-1)/),"integer")
    seq_count = 0
    do i = 0,n_lon-2
        if(.not.finished(0,0,i))
            seq(0,seq_count) = i
            seq(1,seq_count) = 0
            seq_count = seq_count+1
        end if
    end do
    do i = 0,n_lon-2
        if(.not.finished(2,n_lat-2,i))
            seq(0,seq_count) = i
            seq(1,seq_count) = n_lat-2 
            seq_count = seq_count+1
        end if
    end do
    do i = 1,n_lat-3
        if(.not.finished(3,i,0))
            seq(0,seq_count) = 0  
            seq(1,seq_count) = i 
            seq_count = seq_count+1
        end if
    end do
    do i = 1,n_lat-3
        if(.not.finished(1,i,n_lon-2))
            seq(0,seq_count) = n_lon-2
            seq(1,seq_count) = i
            seq_count = seq_count+1
        end if
    end do
    do i = 0,n_lon-2
        if(finished(0,0,i))
            seq(0,seq_count) = i
            seq(1,seq_count) = 0
            seq_count = seq_count+1
        end if
    end do
    do i = 0,n_lon-2
        if(finished(2,n_lat-2,i))
            seq(0,seq_count) = i
            seq(1,seq_count) = n_lat-2 
            seq_count = seq_count+1
        end if
    end do
    do i = 1,n_lat-3
        if(finished(3,i,0))
            seq(0,seq_count) = 0  
            seq(1,seq_count) = i 
            seq_count = seq_count+1
        end if
    end do
    do i = 1,n_lat-3
        if(finished(1,i,n_lon-2))
            seq(0,seq_count) = n_lon-2
            seq(1,seq_count) = i
            seq_count = seq_count+1
        end if
    end do
    do i = 1,n_lat-3
        do j = 1,n_lon-3
            seq(0,seq_count) = j
            seq(1,seq_count) = i
            seq_count = seq_count+1
        end do
    end do
    n_seq = dimsizes(seq(0,:))

    contour = new((/2,60000/),"double") 
    i_cont = 0 
    seq_count = 0
    do while(seq_count.le.n_seq-1)
        x = seq(0,seq_count)
        y = seq(1,seq_count)
        i_s_start = 0 
        if((x.eq.0).and.(.not.finished(3,y,x)))
            i_s_start = 3
        end if
        if((x.eq.n_lon-2).and.(.not.finished(1,y,x)))
            i_s_start = 1
        end if
        if((y.eq.0).and.(.not.finished(0,y,x)))
            i_s_start = 0
        end if
        if((y.eq.n_lat-2).and.(.not.finished(2,y,x)))
            i_s_start = 2
        end if
        do i_s = i_s_start,3
            in_flag = False
            if(.not.finished(i_s,y,x))
                in_flag = True
                if(i_s.eq.0)
                    contour(0,i_cont) = point(i_s,y,x)   
                    contour(1,i_cont) = y  
                    start_cont = i_cont
                    i_cont = i_cont+1
                end if
                if(i_s.eq.1)
                    contour(0,i_cont) = x+1
                    contour(1,i_cont) = point(i_s,y,x)
                    start_cont = i_cont
                    i_cont = i_cont+1
                end if
                if(i_s.eq.2)
                    contour(0,i_cont) = point(i_s,y,x)
                    contour(1,i_cont) = y+1
                    start_cont = i_cont
                    i_cont = i_cont+1
                end if
                if(i_s.eq.3)
                    contour(0,i_cont) = x
                    contour(1,i_cont) = point(i_s,y,x)
                    start_cont = i_cont
                    i_cont = i_cont+1
                end if 
                finished(i_s,y,x) = True 
                flag = 0 
                do k = 0,3   
                    if((.not.finished(k,y,x)))
                        flag = flag+1
                    end if
                end do
                if(flag.eq.1)
                    do k = 0,3
                        if((.not.finished(k,y,x)).and.(k.ne.i_s))    
                            if(k.eq.0)
                                x_next = x
                                y_next = y-1
                                p_next = 2
                                finished(0,y,x) = True   
                                contour(0,i_cont) = point(k,y,x)
                                contour(1,i_cont) = y
                                i_cont = i_cont+1
                            end if
                            if(k.eq.1)
                                x_next = x+1
                                y_next = y
                                p_next = 3
                                finished(1,y,x) = True
                                contour(0,i_cont) = x+1
                                contour(1,i_cont) = point(k,y,x)
                                i_cont = i_cont+1
                            end if
                            if(k.eq.2)
                                x_next = x
                                y_next = y+1
                                p_next = 0
                                finished(2,y,x) = True
                                contour(0,i_cont) = point(k,y,x)
                                contour(1,i_cont) = y+1
                                i_cont = i_cont+1
                            end if
                            if(k.eq.3)
                                x_next = x-1
                                y_next = y
                                p_next = 1
                                finished(3,y,x) = True
                                contour(0,i_cont) = x
                                contour(1,i_cont) = point(k,y,x)
                                i_cont = i_cont+1
                            end if 
                        end if 
                    end do
                else
                    p = i_s
;                    print(pv(y,x))
;                    print(pv(y,x+1))
;                    print(pv(y+1,x+1))
;                    print(pv(y+1,x))
;                    cent_pv = 0.275*(pv(y,x)+pv(y,x+1)+pv(y+1,x+1)+pv(y+1,x))
                    cent_pv = 1.1*contour_pv
;                    print(cent_pv)
                    if(p.eq.0)
                        if((cent_pv-contour_pv)*(contour_pv-pv(y,x)).gt.0)
                            p_in_next = 3
                        else
                            p_in_next = 1
                        end if
                    end if
                    if(p.eq.1)
                        if((cent_pv-contour_pv)*(contour_pv-pv(y,x+1)).gt.0)
                            p_in_next = 0
                        else
                            p_in_next = 2
                        end if
                    end if
                    if(p.eq.2)
                        if((cent_pv-contour_pv)*(contour_pv-pv(y+1,x+1)).gt.0)
                            p_in_next = 1
                        else
                            p_in_next = 3
                        end if
                    end if
                    if(p.eq.3)
                        if((cent_pv-contour_pv)*(contour_pv-pv(y+1,x)).gt.0)
                            p_in_next = 2
                        else
                            p_in_next = 0
                        end if
                    end if
;                    print(p_in_next)
;                    if(p.eq.0)
;                        if(point(1,y,x).lt.point(3,y,x))
;                            p_in_next = 1
;                        else
;                            if(point(1,y,x).gt.point(3,y,x))
;                                p_in_next = 3
;                            else
;                                dist_1 = sqrt((x+1-point(0,y,x))^2+(point(1,y,x)-y)^2)
;                                dist_2 = sqrt((x-point(0,y,x))^2+(point(3,y,x)-y)^2)
;                                if(dist_1.lt.dist_2)
;                                    p_in_next = 1
;                                else
;                                    p_in_next = 3
;                                end if
;                            end if
;                        end if
;                    end if
;                    if(p.eq.1)
;                        if(point(2,y,x).gt.point(0,y,x))
;                            p_in_next = 2
;                        else
;                            if(point(2,y,x).lt.point(0,y,x))
;                                p_in_next = 0
;                            else
;                                dist_1 = sqrt((point(2,y,x)-(x+1))^2+(y+1-point(1,y,x))^2)
;                                dist_2 = sqrt((point(0,y,x)-(x+1))^2+(y-point(1,y,x))^2)
;                                if(dist_1.lt.dist_2)
;                                    p_in_next = 2
;                                else
;                                    p_in_next = 0
;                                end if
;                            end if
;                        end if
;                    end if
;                    if(p.eq.2)
;                        if(point(3,y,x).gt.point(1,y,x))
;                            p_in_next = 3
;                        else
;                            if(point(3,y,x).lt.point(1,y,x))
;                                p_in_next = 1
;                            else
;                                dist_1 = sqrt((x-point(2,y,x))^2+(point(3,y,x)-(y+1))^2)
;                                dist_2 = sqrt((x+1-point(2,y,x))^2+(point(1,y,x)-(y+1))^2)
;                                if(dist_1.lt.dist_2)
;                                    p_in_next = 3
;                                else
;                                    p_in_next = 1
;                                end if
;                            end if
;                        end if
;                    end if
;                    if(p.eq.3)
;                        if(point(0,y,x).lt.point(2,y,x))
;                            p_in_next = 0
;                        else
;                            if(point(0,y,x).gt.point(2,y,x))
;                                p_in_next = 2
;                            else
;                                dist_1 = sqrt((point(0,y,x)-x)^2+(y-point(3,y,x))^2)
;                                dist_2 = sqrt((point(2,y,x)-x)^2+(y+1-point(3,y,x))^2)
;                                if(dist_1.lt.dist_2)
;                                    p_in_next = 0
;                                else
;                                    p_in_next = 2
;                                end if
;                            end if
;                        end if
;                    end if
                    if(p_in_next.eq.0)    
                        x_next = x
                        y_next = y-1
                        p_next = 2
                        finished(0,y,x) = True
                        contour(0,i_cont) = point(0,y,x)   
                        contour(1,i_cont) = y
                        i_cont = i_cont+1
                    end if
                    if(p_in_next.eq.1)
                        x_next = x+1
                        y_next = y
                        p_next = 3
                        finished(1,y,x) = True
                        contour(0,i_cont) = x+1  
                        contour(1,i_cont) = point(1,y,x)
                        i_cont = i_cont+1
                    end if
                    if(p_in_next.eq.2)
                        x_next = x 
                        y_next = y+1
                        p_next = 0
                        finished(2,y,x) = True
                        contour(0,i_cont) = point(2,y,x) 
                        contour(1,i_cont) = y+1
                        i_cont = i_cont+1
                    end if
                    if(p_in_next.eq.3)
                        x_next = x-1
                        y_next = y
                        p_next = 1
                        finished(3,y,x) = True
                        contour(0,i_cont) = x 
                        contour(1,i_cont) = point(3,y,x)
                        i_cont = i_cont+1
                    end if
                end if    
;                print(x)
;                print(y)
;                print(point(:,y,x))
;                print(x_next)
;                print(y_next)
;                print(p_next)
                next_flag = True
                if((x_next.ge.0).and.(x_next.le.n_lon-2).and.\
                   (y_next.ge.0).and.(y_next.le.n_lat-2))
                    next_flag = .not.finished(p_next,y_next,x_next)
;                    print(point(:,y_next,x_next))
;                    print(next_flag)
                else
                    next_flag = False    
;                    print(next_flag)
                end if
                do while(next_flag)
                    x = x_next
                    y = y_next
                    p = p_next
                    flag = 0
                    do k = 0,3  
                        if(.not.finished(k,y,x))    
                            flag = flag+1
                            if(k.ne.p)
                                p_in_next = k
                            end if
                        end if
                    end do
                    finished(p,y,x) = True
                    if(flag.eq.2)
                        if(p_in_next.eq.0)    
                            x_next = x
                            y_next = y-1
                            p_next = 2
                            finished(0,y,x) = True
                            contour(0,i_cont) = point(0,y,x)   
                            contour(1,i_cont) = y
                            i_cont = i_cont+1
                        end if
                        if(p_in_next.eq.1)
                            x_next = x+1
                            y_next = y
                            p_next = 3
                            finished(1,y,x) = True
                            contour(0,i_cont) = x+1  
                            contour(1,i_cont) = point(1,y,x)
                            i_cont = i_cont+1
                        end if
                        if(p_in_next.eq.2)
                            x_next = x 
                            y_next = y+1
                            p_next = 0
                            finished(2,y,x) = True
                            contour(0,i_cont) = point(2,y,x) 
                            contour(1,i_cont) = y+1
                            i_cont = i_cont+1
                        end if
                        if(p_in_next.eq.3)
                            x_next = x-1
                            y_next = y
                            p_next = 1
                            finished(3,y,x) = True
                            contour(0,i_cont) = x 
                            contour(1,i_cont) = point(3,y,x)
                            i_cont = i_cont+1
                        end if
                    else
;                        print(pv(y,x))
;                        print(pv(y,x+1))
;                        print(pv(y+1,x+1))
;                        print(pv(y+1,x))
                        cent_pv = 0.275*(pv(y,x)+pv(y,x+1)+pv(y+1,x+1)+pv(y+1,x))
                        cent_pv = 1.1*contour_pv
;                        print(cent_pv)
                        if(p.eq.0)
                            if((cent_pv-contour_pv)*(contour_pv-pv(y,x)).gt.0)
                                p_in_next = 3
                            else
                                p_in_next = 1
                            end if
                        end if
                        if(p.eq.1)
                            if((cent_pv-contour_pv)*(contour_pv-pv(y,x+1)).gt.0)
                                p_in_next = 0
                            else
                                p_in_next = 2
                            end if
                        end if
                        if(p.eq.2)
                            if((cent_pv-contour_pv)*(contour_pv-pv(y+1,x+1)).gt.0)
                                p_in_next = 1
                            else
                                p_in_next = 3
                            end if
                        end if
                        if(p.eq.3)
                            if((cent_pv-contour_pv)*(contour_pv-pv(y+1,x)).gt.0)
                                p_in_next = 2
                            else
                                p_in_next = 0
                            end if
                        end if
;                        print(p_in_next)
;                        if(p.eq.0)
;                            if(point(1,y,x).lt.point(3,y,x))
;                                p_in_next = 1
;                            else
;                                if(point(1,y,x).gt.point(3,y,x))
;                                    p_in_next = 3
;                                else
;                                    dist_1 = sqrt((x+1-point(0,y,x))^2+(point(1,y,x)-y)^2)
;                                    dist_2 = sqrt((x-point(0,y,x))^2+(point(3,y,x)-y)^2)
;                                    if(dist_1.lt.dist_2)
;                                        p_in_next = 1
;                                    else
;                                        p_in_next = 3
;                                    end if
;                                end if
;                            end if
;                        end if
;                        if(p.eq.1)
;                            if(point(2,y,x).gt.point(0,y,x))
;                                p_in_next = 2
;                            else
;                                if(point(2,y,x).lt.point(0,y,x))
;                                    p_in_next = 0
;                                else
;                                    dist_1 = sqrt((point(2,y,x)-(x+1))^2+(y+1-point(1,y,x))^2)
;                                    dist_2 = sqrt((point(0,y,x)-(x+1))^2+(y-point(1,y,x))^2)
;                                    if(dist_1.lt.dist_2)
;                                        p_in_next = 2
;                                    else
;                                        p_in_next = 0
;                                    end if
;                                end if
;                            end if
;                        end if
;                        if(p.eq.2)
;                            if(point(3,y,x).gt.point(1,y,x))
;                                p_in_next = 3
;                            else
;                                if(point(3,y,x).lt.point(1,y,x))
;                                    p_in_next = 1
;                                else
;                                    dist_1 = sqrt((x-point(2,y,x))^2+(point(3,y,x)-(y+1))^2)
;                                    dist_2 = sqrt((x+1-point(2,y,x))^2+(point(1,y,x)-(y+1))^2)
;                                    if(dist_1.lt.dist_2)
;                                        p_in_next = 3
;                                    else
;                                        p_in_next = 1
;                                    end if
;                                end if
;                            end if
;                        end if
;                        if(p.eq.3)
;                            if(point(0,y,x).lt.point(2,y,x))
;                                p_in_next = 0
;                            else
;                                if(point(0,y,x).gt.point(2,y,x))
;                                    p_in_next = 2
;                                else
;                                    dist_1 = sqrt((point(0,y,x)-x)^2+(y-point(3,y,x))^2)
;                                    dist_2 = sqrt((point(2,y,x)-x)^2+(y+1-point(3,y,x))^2)
;                                    if(dist_1.lt.dist_2)
;                                        p_in_next = 0
;                                    else
;                                        p_in_next = 2
;                                    end if
;                                end if
;                            end if
;                        end if
;                        print(p_in_next)
                        if(p_in_next.eq.0)    
                            x_next = x
                            y_next = y-1
                            p_next = 2
                            finished(0,y,x) = True
                            contour(0,i_cont) = point(0,y,x)   
                            contour(1,i_cont) = y
                            i_cont = i_cont+1
                        end if
                        if(p_in_next.eq.1)
                            x_next = x+1
                            y_next = y
                            p_next = 3
                            finished(1,y,x) = True
                            contour(0,i_cont) = x+1  
                            contour(1,i_cont) = point(1,y,x)
                            i_cont = i_cont+1
                        end if
                        if(p_in_next.eq.2)
                            x_next = x 
                            y_next = y+1
                            p_next = 0
                            finished(2,y,x) = True
                            contour(0,i_cont) = point(2,y,x) 
                            contour(1,i_cont) = y+1
                            i_cont = i_cont+1
                        end if
                        if(p_in_next.eq.3)
                            x_next = x-1
                            y_next = y
                            p_next = 1
                            finished(3,y,x) = True
                            contour(0,i_cont) = x 
                            contour(1,i_cont) = point(3,y,x)
                            i_cont = i_cont+1
                        end if
                    end if
;                    print(x_next)
;                    print(y_next)
;                    print(p_next)
                    next_flag = True
                    if((x_next.ge.0).and.(x_next.le.n_lon-2).and.\
                       (y_next.ge.0).and.(y_next.le.n_lat-2))
                        next_flag = .not.finished(p_next,y_next,x_next)
;                        print(point(:,y_next,x_next))
;                        print(next_flag)
                    else
                        next_flag = False    
;                        print(next_flag)
                    end if
                end do
                i_cont = i_cont+1
            end if
            if(in_flag)
                sx = doubletointeger(floor(min(contour(0,start_cont:i_cont-1))))-2
                ex = doubletointeger(ceil(max(contour(0,start_cont:i_cont-1))))+2
                sy = doubletointeger(floor(min(contour(1,start_cont:i_cont-1))))-2
                ey = doubletointeger(ceil(max(contour(1,start_cont:i_cont-1))))+2
                if(sx.le.0)
                    sx = 0
                end if
                if(ex.ge.n_lon-1)
                    ex = n_lon-1
                end if
                if(sy.le.0)
                    sy = 0
                end if
                if(ey.ge.n_lat-1)
                    ey = n_lat-1
                end if
;                plot_1 = gsn_csm_xy(wks,contour(0,start_cont:i_cont-1)-sx,contour(1,start_cont:i_cont-1)-sy,resxy) 
;                plot_2 = gsn_csm_contour(wks,pv(sy:ey,sx:ex),res)
;                overlay(plot_2,plot_1)
;                draw(plot_2)
            end if
        end do
        seq_count = seq_count+1
    end do

	;将程序找到的PV等值线与NCL画出来的等值线进行对比
;    plot_1 = gsn_csm_xy(wks,contour(0,0:i_cont-1)-140,contour(1,0:i_cont-1)-140,resxy) 
;    plot_2 = gsn_csm_contour(wks,pv(140:320,140:320),res)
    plot_1 = gsn_csm_xy(wks,contour(0,0:i_cont-1),contour(1,0:i_cont-1),resxy) 
    plot_2 = gsn_csm_contour(wks,pv,res)
    overlay(plot_2,plot_1)
    draw(plot_2)
    plot_1 = gsn_csm_xy(wks,contour(0,0:i_cont-1),contour(1,0:i_cont-1),resxy)

end
