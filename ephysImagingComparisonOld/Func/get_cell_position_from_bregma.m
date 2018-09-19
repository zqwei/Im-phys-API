function postion_from_bregma=get_cell_position_from_bregma(RoiCenterPos,fov_origin,angle)
% position_from_bregma : 1x2 vector [ML, AP] coordinate, in um
% fov_origin           : 1x2 vector [ML,AP] cooridnate of the upper left corner of the field of view

um_per_pixel=600/512;
pixx=RoiCenterPos(1)*cos(angle/180*pi)-RoiCenterPos(2)*sin(angle/180*pi);
pixy=RoiCenterPos(1)*sin(angle/180*pi)+RoiCenterPos(2)*cos(angle/180*pi);
postion_from_bregma=fov_origin+[pixx,pixy]*um_per_pixel;

    
    
