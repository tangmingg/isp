clear all;
close all;
% image_left  = imread('picture_lib/sample00.png');
% image_right = imread('picture_lib/sample01.png');
image_left  = imread('picture_lib/im2.png');
image_right = imread('picture_lib/im6.png');
% image_left = imread('picture_lib/example.png');


%%test: plot rgb channel 
image_left_r=uint16(image_left(:,:,1));
image_left_g=uint16(image_left(:,:,2));
image_left_b=uint16(image_left(:,:,3));
image_right_r=uint16(image_right(:,:,1));
image_right_g=uint16(image_right(:,:,2));
image_right_b=uint16(image_right(:,:,3));
image_width  = length(image_left_r(1,:));
image_height = length(image_left_r(:,1));

% figure,
% subplot(2,2,1),imshow(image_left_r),title('Red component');
% subplot(2,2,2),imshow(image_left_g),title('green component');
% subplot(2,2,3),imshow(image_left_b),title('blue component');
% subplot(2,2,4),imagesc(image_left),title('original image');
% figure,
% subplot(2,1,1),imagesc(image_left),title('original image: left');
% subplot(2,1,2),imagesc(image_right),title('original image: right');
dmax_num = floor(image_width/8);
census_ad = zeros(image_height,image_width,dmax_num);

%%Cad
for dmax_n = 1 : dmax_num
    for width_n = 1: image_width
        for height_n = 1: image_height
            dreal = dmax_n - 1;
            if(width_n - dreal <= 0)
                census_ad(height_n,width_n,dmax_n) = 0;
            else
                census_ad(height_n,width_n,dmax_n) = (abs(image_left_r(height_n,width_n) - image_right_r(height_n,width_n - dreal)) ...
                                                   +  abs(image_left_g(height_n,width_n) - image_right_g(height_n,width_n - dreal)) ...
                                                   +  abs(image_left_b(height_n,width_n) - image_right_b(height_n,width_n - dreal)))/3;
            end
        end
    end
end

%%Gray=(R*30+G*59+B*11)/100
image_left_grey  = uint8((uint16(image_left_r).*30 + uint16(image_left_g).*59 + uint16(image_left_b).*11)./100);
image_right_grey = uint8((uint16(image_right_r).*30 + uint16(image_right_g).*59 + uint16(image_right_b).*11)./100);
image_left_grey_mirror_col  = [image_left_grey(:,1)   image_left_grey(:,1)   image_left_grey(:,1)    image_left_grey...
                              image_left_grey(:,end) image_left_grey(:,end) image_left_grey(:,end)];
image_left_grey_mirror_row  = [image_left_grey_mirror_col(1,:)  ;image_left_grey_mirror_col(1,:)   ;image_left_grey_mirror_col;...
                              image_left_grey_mirror_col(end,:);image_left_grey_mirror_col(end,:)];
image_right_grey_mirror_col = [image_right_grey(:,1)   image_right_grey(:,1)   image_right_grey(:,1)    image_right_grey...
                              image_right_grey(:,end) image_right_grey(:,end) image_right_grey(:,end)];
image_right_grey_mirror_row = [image_right_grey_mirror_col(1,:) ;image_right_grey_mirror_col(1,:)   ;image_right_grey_mirror_col;...
                              image_right_grey_mirror_col(end,:);image_right_grey_mirror_col(end,:)];
                     
                                
                          
% figure,
% subplot(2,2,1),imagesc(image_left_grey),title('grey image: left');
% subplot(2,2,2),imagesc(image_right_grey),title('grey image: right');
% subplot(2,2,3),imagesc(image_left_grey_mirror_row),title('grey image: left');
% subplot(2,2,4),imagesc(image_right_grey_mirror_row),title('grey image: right');



i_ca_aver_left = uint8(  (...
                    uint16(image_left_grey_mirror_row(2:1+image_height,4:3+image_width)) + ...
                    uint16(image_left_grey_mirror_row(3:2+image_height,4:3+image_width)) + ...
                    uint16(image_left_grey_mirror_row(4:3+image_height,4:3+image_width)) + ...
                    uint16(image_left_grey_mirror_row(3:2+image_height,3:2+image_width)) + ...
                    uint16(image_left_grey_mirror_row(3:2+image_height,4:3+image_width)) + ...
                    uint16(image_left_grey_mirror_row(3:2+image_height,5:4+image_width)) ...
                    )/6);
i_ca_aver_right = uint8(  (...
                    uint16(image_right_grey_mirror_row(2:1+image_height,4:3+image_width)) + ...
                    uint16(image_right_grey_mirror_row(3:2+image_height,4:3+image_width)) + ...
                    uint16(image_right_grey_mirror_row(4:3+image_height,4:3+image_width)) + ...
                    uint16(image_right_grey_mirror_row(3:2+image_height,3:2+image_width)) + ...
                    uint16(image_right_grey_mirror_row(3:2+image_height,4:3+image_width)) + ...
                    uint16(image_right_grey_mirror_row(3:2+image_height,5:4+image_width)) ...
                    )/6);     
           
window_n = 0;
for window_row_n = 1 : 5
    for window_col_n = 1 : 7
%         window_n = (window_row_n - 1) * 7 + window_col_n;
        if(window_row_n == 1 && window_col_n == 1 || window_row_n == 1 && window_col_n == 7 || ...
                window_row_n == 5 && window_col_n == 1 || window_row_n == 5 && window_col_n == 7)
             window_n =  window_n;
        else
            window_n = window_n + 1;
            census_bin_left(:,:,window_n) = sign(sign(sign(double(i_ca_aver_left)-...
                                            double(image_left_grey_mirror_row(window_row_n:window_row_n+image_height-1,window_col_n:window_col_n+image_width-1)))...
                                            -0.5)+1);
            census_bin_right(:,:,window_n) = sign(sign(sign(double(i_ca_aver_right)-...
                                            double(image_right_grey_mirror_row(window_row_n:window_row_n+image_height-1,window_col_n:window_col_n+image_width-1)))...
                                            - 0.5)+1);
        end
    end
end

census_ca = zeros(image_height,image_width,dmax_num,'uint8');
for dmax_n = 1 : dmax_num
    for width_n = 1: image_width
        for height_n = 1: image_height
            dreal = dmax_n - 1;
            if(width_n - dreal <= 0)
                census_ca(height_n,width_n,dmax_n) = 0;
            else
                for window_n = 1 : 31
%                     census_bin_left(height_n,width_n,window_n)
%                     census_bin_right(height_n,width_n-dreal,window_n)
                    xor_lr = xor(uint8(census_bin_left(height_n,width_n,window_n)),uint8(census_bin_right(height_n,width_n-dreal,window_n)));
                    xor_lr_uint8 = uint8(xor_lr) * 8;
                    census_ca(height_n,width_n,dmax_n) = census_ca(height_n,width_n,dmax_n) + xor_lr_uint8;
                
                end
            end
        end
    end
end

c_fusion = (uint16(census_ca) + uint16(census_ad))/2;
% for nn = 1:dmax_num
%     figure,imagesc(census_ca(:,:,nn))
% end

row_ff = ones(1,image_width+2,'uint16') * 10000;
col_ff = ones(image_height,1,'uint16') * 10000;
inter  = zeros(image_height,image_width);
whole  = [col_ff inter col_ff];
whole  = [row_ff;whole;row_ff];

lr_00 = zeros(image_height+2,image_width+2,dmax_num+2,'uint16');
lr_01 = zeros(image_height+2,image_width+2,dmax_num+2,'uint16');
lr_02 = zeros(image_height+2,image_width+2,dmax_num+2,'uint16');
lr_03 = zeros(image_height+2,image_width+2,dmax_num+2,'uint16');
for nn = 1:dmax_num
    lr_00(:,:,nn) = whole;
    lr_01(:,:,nn) = whole;
    lr_02(:,:,nn) = whole;
    lr_03(:,:,nn) = whole;
end


% for dmax_n = 1 : dmax_num + 1
%     for width_n = 2: image_width + 1
%         for height_n = 2: image_height + 1
%         end
%     end
% end

p1 = 5;
p2 = 30;%4*p1;
e1 = 0.25;%0.4;
e2 = 0.125;%.8;
texture_p = zeros(image_height,image_width,'uint16'); 
image_diff = uint16(abs(image_left_grey_mirror_row(:,2:image_width+6) - image_left_grey_mirror_row(:,1:image_width+5)));
for window_row_n = 1 : 5
    for window_col_n = 1 : 6
        texture_p = texture_p + image_diff(window_row_n:window_row_n+image_height-1,window_col_n:window_col_n+image_width-1);
    end
end
real_p1 = p1 + e1 * (255 - texture_p);
real_p2 = p2 + e2 * (255 - texture_p);
%%00
for dmax_n = 2 : dmax_num + 1
    for width_n = 2: image_width + 1
        for height_n = 2: image_height + 1
            dreal = dmax_n - 2;
            if(dreal == 0)
                lr_disp_2_min = min(lr_00(height_n,width_n-1,3:dmax_num+1));
            elseif(dreal == 1)
                lr_disp_2_min = min(lr_00(height_n,width_n-1,4:dmax_num+1));
            elseif(dreal == 3)
                lr_disp_2_min = min(min(lr_00(height_n,width_n-1,2)),min( lr_00(height_n,width_n-1,5:dmax_num+1)));
            elseif(dreal == dmax_n-1)
                lr_disp_2_min = min(lr_00(height_n,width_n-1,2:dmax_n-1));
            elseif(dreal == dmax_n-2)
                lr_disp_2_min = min(lr_00(height_n,width_n-1,2:dmax_n-2));
            else
                lr_disp_2_min = min(min(lr_00(height_n,width_n-1,2:dmax_n-3)),min(lr_00(height_n,width_n-1,dmax_n:dmax_num+1)));
            end
            lr_disp_equal_0  = uint16(lr_00(height_n,width_n-1,dmax_n));
            lr_disp_equal_s0 = uint16(lr_00(height_n,width_n-1,dmax_n+1)) + uint16(real_p1(height_n-1,width_n-1));
            lr_disp_equal_s1 = uint16(lr_00(height_n,width_n-1,dmax_n-1)) + uint16(real_p1(height_n-1,width_n-1));
            lr_disp_equal_2  = uint16(lr_disp_2_min) + uint16(real_p2(height_n-1,width_n-1));           
            
            lr_00(height_n,width_n,dmax_n) = uint16(c_fusion(height_n-1,width_n-1,dmax_n-1)) + ...
            uint16(min([lr_disp_equal_0 , lr_disp_equal_s0, lr_disp_equal_s1, lr_disp_equal_2]))...
            -uint16(min(lr_00(height_n,width_n-1,dmax_n:2:dmax_num+1)));
            end
        end
end
%01
for dmax_n = 2 : dmax_num + 1
    for width_n = 2: image_width + 1
        for height_n = 2: image_height + 1
            dreal = dmax_n - 2;
            if(dreal == 0)
                lr_disp_2_min = min(lr_01(height_n-1,width_n-1,3:dmax_num+1));
            elseif(dreal == 1)
                lr_disp_2_min = min(lr_01(height_n-1,width_n-1,4:dmax_num+1));
            elseif(dreal == 3)
                lr_disp_2_min = min(min(lr_01(height_n-1,width_n-1,2)),min( lr_01(height_n-1,width_n-1,5:dmax_num+1)));
            elseif(dreal == dmax_n-1)
                lr_disp_2_min = min(lr_01(height_n-1,width_n-1,2:dmax_n-1));
            elseif(dreal == dmax_n-2)
                lr_disp_2_min = min(lr_01(height_n-1,width_n-1,2:dmax_n-2));
            else
                lr_disp_2_min = min(min(lr_01(height_n-1,width_n-1,2:dmax_n-3)),min(lr_01(height_n-1,width_n-1,dmax_n:dmax_num+1)));
            end
            lr_disp_equal_0  = uint16(lr_01(height_n-1,width_n-1,dmax_n));
            lr_disp_equal_s0 = uint16(lr_01(height_n-1,width_n-1,dmax_n+1)) + uint16(p1);
            lr_disp_equal_s1 = uint16(lr_01(height_n-1,width_n-1,dmax_n-1)) + uint16(p1);
            lr_disp_equal_2  = uint16(lr_disp_2_min + p2);           
            
            lr_01(height_n,width_n,dmax_n) = uint16(c_fusion(height_n-1,width_n-1,dmax_n-1)) + ...
            uint16(min([lr_disp_equal_0 , lr_disp_equal_s0, lr_disp_equal_s1, lr_disp_equal_2]))...
            -uint16(min(lr_01(height_n-1,width_n-1,dmax_n:2:dmax_num+1)));
            end
        end
end

%02
for dmax_n = 2 : dmax_num + 1
    for width_n = 2: image_width + 1
        for height_n = 2: image_height + 1
            dreal = dmax_n - 2;
            if(dreal == 0)
                lr_disp_2_min = min(lr_02(height_n-1,width_n,3:dmax_num+1));
            elseif(dreal == 1)
                lr_disp_2_min = min(lr_02(height_n-1,width_n,4:dmax_num+1));
            elseif(dreal == 3)
                lr_disp_2_min = min(min(lr_02(height_n-1,width_n,2)),min( lr_02(height_n-1,width_n,5:dmax_num+1)));
            elseif(dreal == dmax_n-1)
                lr_disp_2_min = min(lr_02(height_n-1,width_n,2:dmax_n-1));
            elseif(dreal == dmax_n-2)
                lr_disp_2_min = min(lr_02(height_n-1,width_n,2:dmax_n-2));
            else
                lr_disp_2_min = min(min(lr_02(height_n-1,width_n,2:dmax_n-3)),min(lr_02(height_n-1,width_n,dmax_n:dmax_num+1)));
            end
            lr_disp_equal_0  = uint16(lr_02(height_n-1,width_n,dmax_n));
            lr_disp_equal_s0 = uint16(lr_02(height_n-1,width_n,dmax_n+1)) + uint16(p1);
            lr_disp_equal_s1 = uint16(lr_02(height_n-1,width_n,dmax_n-1)) + uint16(p1);
            lr_disp_equal_2  = uint16(lr_disp_2_min) + uint16(p2);           
            
            lr_02(height_n,width_n,dmax_n) = uint16(c_fusion(height_n-1,width_n-1,dmax_n-1)) + ...
            uint16(min([lr_disp_equal_0 , lr_disp_equal_s0, lr_disp_equal_s1, lr_disp_equal_2]))...
            -uint16(min(lr_02(height_n-1,width_n,dmax_n:2:dmax_num+1)));
            end
        end
end

%03
for dmax_n = 2 : dmax_num + 1
    for width_n = 2: image_width + 1
        for height_n = 2: image_height + 1
            dreal = dmax_n - 2;
            if(dreal == 0)
                lr_disp_2_min = min(lr_03(height_n-1,width_n+1,3:dmax_num+1));
            elseif(dreal == 1)
                lr_disp_2_min = min(lr_03(height_n-1,width_n+1,4:dmax_num+1));
            elseif(dreal == 3)
                lr_disp_2_min = min(min(lr_03(height_n-1,width_n+1,2)),min( lr_03(height_n-1,width_n+1,5:dmax_num+1)));
            elseif(dreal == dmax_n-1)
                lr_disp_2_min = min(lr_03(height_n-1,width_n+1,2:dmax_n-1));
            elseif(dreal == dmax_n-2)
                lr_disp_2_min = min(lr_03(height_n-1,width_n+1,2:dmax_n-2));
            else
                lr_disp_2_min = min(min(lr_03(height_n-1,width_n+1,2:dmax_n-3)),min(lr_03(height_n-1,width_n+1,dmax_n:dmax_num+1)));
            end
            lr_disp_equal_0  = uint16(lr_03(height_n-1,width_n+1,dmax_n));
            lr_disp_equal_s0 = uint16(lr_03(height_n-1,width_n+1,dmax_n+1)) + uint16(p1);
            lr_disp_equal_s1 = uint16(lr_03(height_n-1,width_n+1,dmax_n-1)) + uint16(p1);
            lr_disp_equal_2  = uint16(lr_disp_2_min) + uint16(p2);           
            
            lr_03(height_n,width_n,dmax_n) = uint16(c_fusion(height_n-1,width_n-1,dmax_n-1)) + ...
            uint16(min([lr_disp_equal_0 , lr_disp_equal_s0, lr_disp_equal_s1, lr_disp_equal_2]))...
            -uint16(min(lr_03(height_n-1,width_n+1,dmax_n:2:dmax_num+1)));
            end
        end
end

spd = lr_00(2:image_height+1,2:image_width+1,2:dmax_num+1) + lr_01(2:image_height+1,2:image_width+1,2:dmax_num+1) + ...
    lr_02(2:image_height+1,2:image_width+1,2:dmax_num+1) + lr_03(2:image_height+1,2:image_width+1,2:dmax_num+1);


coef = 0.85;
for width_n = 1: image_width
    for height_n = 1: image_height
        list_origin_reshape = reshape(spd(height_n,width_n,:),1,dmax_num);
        [dmin_0,dindex_0] = min(list_origin_reshape);
        list_second_reshape = list_origin_reshape;
        list_second_reshape(dindex_0) = 65535;
        [dmin_1,dindex_1] = min(list_second_reshape);
        if(dmin_0 > coef * dmin_1)
            disp_d(height_n,width_n) = 0;
        else
            if(dindex_0 == dmax_num)
                disp_d(height_n,width_n) = dindex_0 - 1;
            elseif(dindex_0 > 1)
                disp_d(height_n,width_n) = dindex_0 - 1 ...
                                         - (list_origin_reshape(dindex_0+1) - list_origin_reshape(dindex_0-1))...
                                         / (2*(list_origin_reshape(dindex_0+1) + list_origin_reshape(dindex_0-1) - 2* list_origin_reshape(dindex_0)));
            end
%             disp_d(height_n,width_n) = dindex_0 - 1;
        end
    end
end
        
% for nn = 1:dmax_num
%     figure,imagesc(spd(:,:,nn))
% end


figure,imagesc(disp_d)
sss = 1

