function [best_para_fit_pv best_para_fit_rms] = Zernike_Fit_Parabold(rpt_path_name, pch_path_name, surf_name, Surf_VC, Surf_c,X0,Y0,radius,termOrder,firstTermIndex,kArray,surfacetype)
% rpt_path_name:  face_nodes file name
% pch_path_name: FEA Displement output files name
% surface_name:  Output files name
% Surf_VC: vertex coordibat vector (输入顶点坐标系向量)
% Surf_c: mirror surface curvature
% vc:  vertex coordinate
% gc:  global coordinate



%%%%% Read the original data and perform coordinate transformation  %%%%% 
% ‘rpt_path_name’（初始节点坐标文件）；% ‘pch_path_name’（节点位移文件）；‘Surf_c’镜面顶点曲率(曲率半径导数)；‘surfacetype’镜面类型；
%'radius'镜面半径；(X0,Y0)归一化半径坐标，存在离心时应加上离心距离；‘termOrder’泽尼克类型；firstTermIndex-泽尼克初始阶数；‘kmax’-泽尼克最大阶数
rpt_path_name='ParaboloidD120DisPM.rpt';pch_path_name='ParaboloidD120YDisPM.pch';surf_name='ParaboloidD120YDisPM';Surf_VC=[0 0 0;0 0 100;100 0 0]';Surf_c=1/750;surfacetype='paraboloid';
radius=125/2;X0=123;Y0=0;termOrder='FringeOrder'; firstTermIndex=1;kmax=36;
% AreaFactor_pch_name='Weight.pch';
%Read x,y,z of mirror nodes from Hypermesh output files (xxx.rpt)
Face_nodes = surface_byfile(surfacetype, rpt_path_name, Surf_VC, Surf_c, 0, 1);    %Face_nodes(:,1) is nodes ID; Face_nodes(:,2:4) is x,y,z values

%% Define the coordinate And Calculate the transformation matrix %%
%     %% New Method %%
%     %%% Define the global 'coordinate 'g_P & the vertex coordinate 'v_P' %%%
%     g_P=eye(3);    %%%% the global coordinate 'g_P'
%     VC_z = Surf_VC(:,2) - Surf_VC(:,1);
%     VC_z = VC_z/norm(VC_z);
%     VC_y = cross(VC_z, Surf_VC(:,3) - Surf_VC(:,1));
%     VC_y = VC_y/norm(VC_y);
%     VC_x = cross(VC_y,VC_z);
%     VC_x = VC_x/norm(VC_x);
%     v_P=[VC_x,VC_y,VC_z];   %%%% the vertex coordinate 'v_P'
%     %%%%%%%%%% Calculate the coordinate transformation matrix   %%%%%%%%%%%%%%%
%     v_P_g=-Surf_VC(:,1);       %%%%%%%%%%%  平动矩阵-> 'v_P_g'
%     c_P=v_P-repmat(v_P_g,1,3);   %%%% 中转坐标系-> 'c_P'
%     v_R_g=c_P/g_P;     %%%%% 转动矩阵->‘v_R_g’
%     Tranform_gc2vc=[v_R_g,v_P_g;zeros(1,3),eye(1)];   %%%%% 全局转局部矩阵 Tranform_gc2vc
    %% Old Method %%
    VC_z = Surf_VC(:,2) - Surf_VC(:,1);
    VC_z = VC_z/norm(VC_z);
    VC_y = cross(VC_z, Surf_VC(:,3) - Surf_VC(:,1));
    VC_y = VC_y/norm(VC_y);
    VC_x = cross(VC_y,VC_z);
    VC_x = VC_x/norm(VC_x);
    %global coordinate system To vertex coordinate system
    Rotate_gc2vc = eye(4);
    Rotate_gc2vc(1:3,1:3) = eye(3)/[VC_x, VC_y, VC_z];
    Translate_gc2vc = eye(4);
    Translate_gc2vc(1:3,4) = -Surf_VC(:,1);
    Tranform_gc2vc = Rotate_gc2vc*Translate_gc2vc;
%%  ************************ END ************************  %%




%%Read displement of mirror nodes from Optistruct output files (xxx.pch)%% 
fid_punch = fopen(pch_path_name);
i_s = 1;                               % number of subcase
subcase_name = cell(1,1);
Dis = cell(1,1);                   % cell of displacements
while ~feof(fid_punch)     % 文件指针达到文件末尾时，表达式值为“假”，否则为“真”
    tline = fgetl(fid_punch);
    if length(tline)>6 && strcmp(tline(1:6) , '$LABEL')
        temp = tline(12:end);
        temp(isspace(temp)) = [];
        subcase_name{i_s} = temp;
    elseif length(tline)>13 && strcmp(tline(1:14) , '$DISPLACEMENTS')
        tline = fgetl(fid_punch);                       % Don`t delete
        tline = fgetl(fid_punch);
        if length(tline)>11 && strcmp(tline(1:11) , '$SUBCASE ID')
            i = 1;
            Dis_i_s = zeros(5,4);
            while i <= length(Face_nodes)
                tline = fgetl(fid_punch);
                gridid = str2double(tline(1:10));
                if ~isnan(gridid)
                    if ismember(gridid, Face_nodes(:,1))                  % 判断gridid是否在Face_nodes(:,1)集合中出现
                        Dis_i_s(i, 1) = gridid;
                        Dis_i_s(i, 2) = str2double(tline(21:40));
                        Dis_i_s(i, 3) = str2double(tline(41:56));
                        Dis_i_s(i, 4) = str2double(tline(57:72));
                        i = i + 1;
                    end
                end
            end
            Dis{i_s} = Dis_i_s;
            i_s = i_s + 1;
        end
    end
end
fclose(fid_punch);
% verify if Dis_i_s compatible with Face_nodes
for i = 1:length(Dis)
    if sum(abs(Dis{i}(:,1) - Face_nodes(:,1))) > 0
        error('error in file reading')
    end
end
%%%%%%%%%% ********* %%%%%%%%%%%



%%%%%%% Original point coordinates and transformed point coordinates  %%%%%%
Undeform_gc_xyz=Face_nodes(:,2:4);
% Transform to Vertices coordinate
Undeform_vc_xyz = Tranform_gc2vc*[transpose(Undeform_gc_xyz); ones(1,length(Undeform_gc_xyz))];     %初始节点转化坐标
Undeform_vc_xyz = transpose(Undeform_vc_xyz); 

Deform_gc_xyz =cell(1,1);
for i_s=1:length(Dis)                                                                                                                                                                                                 
    Deform_gc_xyz{i_s}(:,1)=Dis{i}(:,1);
    Deform_gc_xyz{i_s}(:,2)=Face_nodes(:,2)+Dis{i_s}(:,2);
    Deform_gc_xyz{i_s}(:,3)=Face_nodes(:,3)+Dis{i_s}(:,3);
    Deform_gc_xyz{i_s}(:,4)=Face_nodes(:,4)+Dis{i_s}(:,4);
end

Deform_vc_xyz=cell(1,1);
for i_s=1:length(Dis)
    Deform_vc_xyz{i_s} =transpose(Tranform_gc2vc*[transpose(Deform_gc_xyz{i_s}(:,2:4)); ones(1,length(Deform_gc_xyz{i_s}(:,2:4)))]);     %变形后节点全局坐标系转化顶点坐标系
end
%%%%%% ************ %%%%%%%%



%%%%%%%%%%%%%%% Substraction rigid body motion & Calculate deformation %%%%%%%%%%%%%%%%%
% A_rm=cell(1,1);
% b_rm=cell(1,1);
% d_x=cell(1,1);
% d_y=cell(1,1);
% d_z=cell(1,1);
% 
% % for i_s=1:length(Dis)                               % Dis 表示工况数个数
% %     for j=1:length(Face_nodes)
% %         A_rm{i_s}(6*j-5,1)= Undeform_vc_xyz(j,2)^2+Undeform_vc_xyz(j,3)^2;
% %         A_rm{i_s}(6*j-5,2)=-1*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,2);
% %         A_rm{i_s}(6*j-5,3)=-1*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-5,4)= 0;
% %         A_rm{i_s}(6*j-5,5)=-1*Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-5,6)=Undeform_vc_xyz(j,2);
% % 
% %         A_rm{i_s}(6*j-4,1)=-1*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,2);
% %         A_rm{i_s}(6*j-4,2)=Undeform_vc_xyz(j,1)^2+Undeform_vc_xyz(j,3)^2;
% %         A_rm{i_s}(6*j-4,3)=-1*Undeform_vc_xyz(j,2)*Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-4,4)= Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-4,5)= 0;
% %         A_rm{i_s}(6*j-4,6)=-1*Undeform_vc_xyz(j,1);
% % 
% %         A_rm{i_s}(6*j-3,1)=-1*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-3,2)=-1*Undeform_vc_xyz(j,2)*Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-3,3)=Undeform_vc_xyz(j,1)^2+Undeform_vc_xyz(j,2)^2;
% %         A_rm{i_s}(6*j-3,4)=-1*Undeform_vc_xyz(j,2);
% %         A_rm{i_s}(6*j-3,5)=Undeform_vc_xyz(j,1);
% %         A_rm{i_s}(6*j-3,6)= 0;
% %     
% %         A_rm{i_s}(6*j-2,1)=0;
% %         A_rm{i_s}(6*j-2,2)=Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-2,3)=-1*Undeform_vc_xyz(j,2);
% %         A_rm{i_s}(6*j-2,4)=1;
% %         A_rm{i_s}(6*j-2,5)=0;
% %         A_rm{i_s}(6*j-2,6)=0;
% %     
% %         A_rm{i_s}(6*j-1,1)=-1*Undeform_vc_xyz(j,3);
% %         A_rm{i_s}(6*j-1,2)=0;
% %         A_rm{i_s}(6*j-1,3)=Undeform_vc_xyz(j,1);
% %         A_rm{i_s}(6*j-1,4)=0;
% %         A_rm{i_s}(6*j-1,5)=1;
% %         A_rm{i_s}(6*j-1,6)=0;
% % 
% %         A_rm{i_s}(6*j,1)=Undeform_vc_xyz(j,2);
% %         A_rm{i_s}(6*j,2)=-1*Undeform_vc_xyz(j,1);
% %         A_rm{i_s}(6*j,3)=0;
% %         A_rm{i_s}(6*j,4)=0;
% %         A_rm{i_s}(6*j,5)=0;
% %         A_rm{i_s}(6*j,6)=1;
% % 
% %         b_rm{i_s}(6*j-5,1)=Undeform_vc_xyz(j,2)*Deform_vc_xyz{i_s}(j,3)-Undeform_vc_xyz(j,3)*Deform_vc_xyz{i_s}(j,2);
% %         b_rm{i_s}(6*j-4,1)=Undeform_vc_xyz(j,3)*Deform_vc_xyz{i_s}(j,1)-Undeform_vc_xyz(j,1)* Deform_vc_xyz{i_s}(j,3);
% %         b_rm{i_s}(6*j-3,1)=Undeform_vc_xyz(j,1)*Deform_vc_xyz{i_s}(j,2)-Undeform_vc_xyz(j,2)*Deform_vc_xyz{i_s}(j,1);
% %         b_rm{i_s}(6*j-2,1)=Deform_vc_xyz{i_s}(j,1)-Undeform_vc_xyz(j,1);
% %         b_rm{i_s}(6*j-1,1)=Deform_vc_xyz{i_s}(j,2)-Undeform_vc_xyz(j,2);
% %         b_rm{i_s}(6*j,1)   = Deform_vc_xyz{i_s}(j,3)-Undeform_vc_xyz(j,3);
% %     end
% % end
% 
% for i_s=1:length(Dis)
%     for i=1:6
%         for j=1:6
%             A_rm{i_s}(i,j)=0;
%         end
%          b_rm{i_s}(i,1)=0;
%     end
% end
% 
% % W=Area_factors;
% W = Area_factor_byfiles(AreaFactor_pch_name,Face_nodes);        % 读取镜面面积因子
W=zeros(3,4);           %%% Node weight
for j=1:length(Undeform_vc_xyz)
    W(j,4)=1;
end
% for i_s=1:length(Dis)                               % Dis 表示工况数个数
%     for j=1:length(Face_nodes)
%          
%         A_rm{i_s}(1,1)= W(j,4)*(Undeform_vc_xyz(j,2)^2+Undeform_vc_xyz(j,3)^2)+A_rm{i_s}(1,1);
%         A_rm{i_s}(1,2)= W(j,4)*(-1)*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,2)+A_rm{i_s}(1,2);
%         A_rm{i_s}(1,3)= W(j,4)*(-1)*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,3)+A_rm{i_s}(1,3);
%         A_rm{i_s}(1,4)= 0;
%         A_rm{i_s}(1,5)= W(j,4)*(-1)*Undeform_vc_xyz(j,3)+A_rm{i_s}(1,5);
%         A_rm{i_s}(1,6)= W(j,4)*Undeform_vc_xyz(j,2)+A_rm{i_s}(1,6);
% 
%         A_rm{i_s}(2,1)= W(j,4)*(-1)*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,2)+A_rm{i_s}(2,1);
%         A_rm{i_s}(2,2)= W(j,4)*(Undeform_vc_xyz(j,1)^2+Undeform_vc_xyz(j,3)^2)+A_rm{i_s}(2,2);
%         A_rm{i_s}(2,3)= W(j,4)*(-1)*Undeform_vc_xyz(j,2)*Undeform_vc_xyz(j,3);
%         A_rm{i_s}(2,4)= W(j,4)*Undeform_vc_xyz(j,3)+A_rm{i_s}(2,4);
%         A_rm{i_s}(2,5)= 0;
%         A_rm{i_s}(2,6)= W(j,4)*(-1)*Undeform_vc_xyz(j,1)+A_rm{i_s}(2,6);
% 
%         A_rm{i_s}(3,1)= W(j,4)*(-1)*Undeform_vc_xyz(j,1)*Undeform_vc_xyz(j,3)+A_rm{i_s}(3,1);
%         A_rm{i_s}(3,2)= W(j,4)*(-1)*Undeform_vc_xyz(j,2)*Undeform_vc_xyz(j,3)+A_rm{i_s}(3,2);
%         A_rm{i_s}(3,3)= W(j,4)*(Undeform_vc_xyz(j,1)^2+Undeform_vc_xyz(j,2)^2)+A_rm{i_s}(3,3);
%         A_rm{i_s}(3,4)= W(j,4)*(-1)*Undeform_vc_xyz(j,2)+A_rm{i_s}(3,4);
%         A_rm{i_s}(3,5)= W(j,4)*Undeform_vc_xyz(j,1)+A_rm{i_s}(3,5);
%         A_rm{i_s}(3,6)= 0;
%     
%         A_rm{i_s}(4,1)= 0;
%         A_rm{i_s}(4,2)= W(j,4)*Undeform_vc_xyz(j,3)+A_rm{i_s}(4,2);
%         A_rm{i_s}(4,3)= W(j,4)*(-1)*Undeform_vc_xyz(j,2)+A_rm{i_s}(4,3);
%         A_rm{i_s}(4,4)= W(j,4)*1+A_rm{i_s}(4,4);
%         A_rm{i_s}(4,5)= 0;
%         A_rm{i_s}(4,6)= 0;
%     
%         A_rm{i_s}(5,1)= W(j,4)*(-1)*Undeform_vc_xyz(j,3)+A_rm{i_s}(5,1);
%         A_rm{i_s}(5,2)= 0;
%         A_rm{i_s}(5,3)= W(j,4)*Undeform_vc_xyz(j,1)+A_rm{i_s}(5,3);
%         A_rm{i_s}(5,4)= 0;
%         A_rm{i_s}(5,5)= W(j,4)*1+A_rm{i_s}(5,5);
%         A_rm{i_s}(5,6)= 0;
% 
%         A_rm{i_s}(6,1)= W(j,4)*Undeform_vc_xyz(j,2)+A_rm{i_s}(6,1);
%         A_rm{i_s}(6,2)= W(j,4)*(-1)*Undeform_vc_xyz(j,1)+A_rm{i_s}(6,2);
%         A_rm{i_s}(6,3)= 0;
%         A_rm{i_s}(6,4)= 0;
%         A_rm{i_s}(6,5)= 0;
%         A_rm{i_s}(6,6)= W(j,4)*1+A_rm{i_s}(6,6);
% 
%         b_rm{i_s}(1,1)= W(j,4)*(Undeform_vc_xyz(j,2)*Deform_vc_xyz{i_s}(j,3)-Undeform_vc_xyz(j,3)*Deform_vc_xyz{i_s}(j,2))+b_rm{i_s}(1,1);
%         b_rm{i_s}(2,1)= W(j,4)*(Undeform_vc_xyz(j,3)*Deform_vc_xyz{i_s}(j,1)-Undeform_vc_xyz(j,1)* Deform_vc_xyz{i_s}(j,3))+b_rm{i_s}(2,1);
%         b_rm{i_s}(3,1)= W(j,4)*(Undeform_vc_xyz(j,1)*Deform_vc_xyz{i_s}(j,2)-Undeform_vc_xyz(j,2)*Deform_vc_xyz{i_s}(j,1))+b_rm{i_s}(3,1);
%         b_rm{i_s}(4,1)= W(j,4)*(Deform_vc_xyz{i_s}(j,1)-Undeform_vc_xyz(j,1))+b_rm{i_s}(4,1);
%         b_rm{i_s}(5,1)= W(j,4)*(Deform_vc_xyz{i_s}(j,2)-Undeform_vc_xyz(j,2))+b_rm{i_s}(5,1);
%         b_rm{i_s}(6,1)= W(j,4)*(Deform_vc_xyz{i_s}(j,3)-Undeform_vc_xyz(j,3))+b_rm{i_s}(6,1);
%     end
% end
% 
% for i=1:length(Dis)
%     rigid_motion{i}=Least_square('Pseudoinverse Matrix',A_rm{i},b_rm{i});
% end
% 
% % rigid_motion=cell(1,1);
% % rigid_motion={Sol};


for i_s=1:length(Dis)
     x1=Undeform_vc_xyz (:,1);
     y1=Undeform_vc_xyz(:,2);
     z1=Undeform_vc_xyz(:,3);
     x2=Deform_vc_xyz{i_s}(:,1);
     y2=Deform_vc_xyz{i_s}(:,2);
     z2=Deform_vc_xyz{i_s}(:,3);
      
    [rigid_motion{i_s}(4,1),rigid_motion{i_s}(5,1),rigid_motion{i_s}(6,1),rigid_motion{i_s}(1,1),rigid_motion{i_s}(2,1),rigid_motion{i_s}(3,1)]=Quatstiting_xy( x1,y1,z1,x2,y2,z2); 
end


for i=1:length(Dis)
    d_x{i}=Deform_vc_xyz{i}(:,1)-Undeform_vc_xyz(:,1)+rigid_motion{i}(3,1).*Undeform_vc_xyz(:,2)-rigid_motion{i}(2,1).*Undeform_vc_xyz(:,3)-rigid_motion{i}(4,1);
    d_y{i}=Deform_vc_xyz{i}(:,2)-Undeform_vc_xyz(:,2)+rigid_motion{i}(1,1).*Undeform_vc_xyz(:,3)-rigid_motion{i}(3,1).*Undeform_vc_xyz(:,1)-rigid_motion{i}(5,1);
    d_z{i}=Deform_vc_xyz{i}(:,3)-Undeform_vc_xyz(:,3)+rigid_motion{i}(2,1).*Undeform_vc_xyz(:,1)-rigid_motion{i}(1,1).*Undeform_vc_xyz(:,2)-rigid_motion{i}(6,1);
end
%%%%%%%%%%%%%%%%%%  *****************   %%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%  Saggital height displacement correction  %%%%%%%%%%
cur=Surf_c;                    % Vertex curvature; 
UR=cell(1,1);
DR=cell(1,1);
delta_sagz=cell(1,1);
sagz=cell(1,1);

switch surfacetype        % selcet mirror surface type：1-paraboloid；2-sphere；3-plane
    case 'paraboloid'   
        for i_s=1:length(Dis)
            UR{i_s}=sqrt(Undeform_vc_xyz(:,1).^2 + Undeform_vc_xyz(:,2).^2);
            DR{i_s}=sqrt((Undeform_vc_xyz(:,1)+d_x{i_s}(:,1)).^2 + (Undeform_vc_xyz(:,2)+d_y{i_s}(:,1)).^2);
            delta_sagz{i_s}=0.5.*cur.*DR{i_s}(:,1).^2-0.5.*cur.*UR{i_s}(:,1).^2;
            sagz{i_s}=-delta_sagz{i_s}+d_z{i_s};
        end
    case 'sphere'
        for i_s=1:length(Dis) 
            UR{i_s}=sqrt(Undeform_vc_xyz(:,1).^2 + Undeform_vc_xyz(:,2).^2);
            DR{i_s}=sqrt((Undeform_vc_xyz(:,1)+d_x{i_s}(:,1)).^2 + (Undeform_vc_xyz(:,2)+d_y{i_s}(:,1)).^2);
            delta_sagz{i_s}=cur.*DR{i_s}(:,1).^2./( 1+sqrt(1 - cur.^2*DR{i_s}(:,1).^2))-cur.*UR{i_s}(:,1).^2./( 1+sqrt(1 - cur.^2*UR{i_s}(:,1).^2));
            sagz{i_s}=-delta_sagz{i_s}+d_z{i_s};
        end  
    case 'plane'
        for i_s=1:length(Dis)
            sagz{i_s}=d_z{i_s};
        end
end
%%%%%%%%%%%%%%%%%%%%% ********  %%%%%%%%%%%%%%%%%%%%%%%



%%%%%% output 'Rigid body displacement'& 'PV'&'RMS'   *********
best_para_fit_rms=zeros(1,1);
best_para_fit_pv=zeros(1,1);

res_rpt_1 = [surf_name,'_fitres.rpt'];           % Output rigid body motion into a file 
fid_res_1 = fopen(res_rpt_1,'wt');
fprintf(fid_res_1, '%s\n',surf_name);

% NN=length(Face_nodes); 
% fprintf(fid_res_1, 'NN = %d\n',subcase_name{i}, NN);
    
for i=1:length(Dis)
   
    fprintf(fid_res_1,'SubCase : %s\n',subcase_name{i});
    rigid_motion_NU=rigid_motion{i}(1:6,1);
    rigid_motion{i}(1:3,1)=rigid_motion{i}(1:3,1).*180.*3600./pi;   % 弧度—>角度
    rigid_motion{i}(4:6,1)=rigid_motion{i}(4:6,1).*1000;                 % um
    best_para_fit_rms(i)=std(sagz{i}(:,1),1)*1000000;      % nm
    best_para_fit_pv(i)  =(max(sagz{i}(:,1))-min(sagz{i}(:,1)))*1000000;
    
    fprintf(fid_res_1, 'PV_%s  = %10.3enm (%10.3emm)\n',subcase_name{i}, best_para_fit_pv(i),best_para_fit_pv(i)/1000000);
    fprintf(fid_res_1, 'rms_%s = %10.3enm (%10.3emm)\n',subcase_name{i}, best_para_fit_rms(i),best_para_fit_rms(i)/1000000);
    fprintf(fid_res_1, 'Δθx_%s = %10.3e" (%10.3erad)\n', subcase_name{i},rigid_motion{i}(1), rigid_motion_NU(1,1));
    fprintf(fid_res_1, 'Δθy_%s = %10.3e" (%10.3erad) \n', subcase_name{i},rigid_motion{i}(2),rigid_motion_NU(2,1));
    fprintf(fid_res_1, 'Δθz_%s = %10.3e" (%10.3erad)\n', subcase_name{i},rigid_motion{i}(3),rigid_motion_NU(3,1));
    fprintf(fid_res_1, 'ΔX_%s  = %10.3eum (%10.3emm)\n', subcase_name{i},rigid_motion{i}(4),rigid_motion_NU(4,1));
    fprintf(fid_res_1, 'ΔY_%s  = %10.3eum (%10.3emm)\n', subcase_name{i},rigid_motion{i}(5),rigid_motion_NU(5,1));
    fprintf(fid_res_1, 'ΔZ_%s  = %10.3eum (%10.3emm)\n\n\n', subcase_name{i},rigid_motion{i}(6),rigid_motion_NU(6,1));
end
fclose(fid_res_1); 
%% Rigid body displacement output to Excel %%
RBD=[];
for i=1:length(Dis)
    RBD(i,1)=rigid_motion_NU(4,1);
    RBD(i,2)=rigid_motion_NU(5,1);
    RBD(i,3)=rigid_motion_NU(6,1);
    RBD(i,4)=rigid_motion_NU(1,1);
    RBD(i,5)=rigid_motion_NU(2,1);
    RBD(i,6)=rigid_motion_NU(3,1);
    RBD(i,7)=best_para_fit_rms(i)/1000000;
    RBD(i,8)=best_para_fit_pv(i)/1000000;
end
% Time=0:0.01:Num_Time;
% Time=num2cell(Time');
RBD2=num2cell(RBD);
Title={'Tx(mm)','Ty(mm)','Tz(mm)','Rx(rad)','Ry(rad)','Rz(rad)','RMS(mm)','PV(mm)'};
% Title={'Time(s)','Tx(mm)','Ty(mm)','Tz(mm)','Rx(rad)','Ry(rad)','Rz(rad)'};
Data_RBD_Coffeicients=cell(length(RBD(:,1))+1,8);
Data_RBD_Coffeicients(1,:)=Title;
% Data_RBD(2:end,1)=Time(1:end,1);
Data_RBD_Coffeicients(2:end,1)=RBD2(1:end,1);
Data_RBD_Coffeicients(2:end,2)=RBD2(1:end,2);
Data_RBD_Coffeicients(2:end,3)=RBD2(1:end,3);
Data_RBD_Coffeicients(2:end,4)=RBD2(1:end,4);
Data_RBD_Coffeicients(2:end,5)=RBD2(1:end,5);
Data_RBD_Coffeicients(2:end,6)=RBD2(1:end,6);
Data_RBD_Coffeicients(2:end,7)=RBD2(1:end,7);
Data_RBD_Coffeicients(2:end,8)=RBD2(1:end,8);
filename1 = 'RBD.xlsx';
xlswrite(filename1,Data_RBD_Coffeicients);
res_rpt_4=filename1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ************  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

%% Zernike Fitting and Output cofficients files
% Calculate Zernike cofficents in each subcases.
% X0=0;
% Y0=0;
% radius=400;
WW=W(:,4);         %%%%%%%%改1
% termOrder='FringeOrder';
% firstTermIndex=1;kmax=37;
kArray=firstTermIndex:kmax;
bb_rm=zeros(length(kArray),1);
% W=cell(1,1);
% for i_s=1:length(Dis)
    % Transform to Vertices coordinate
    coefficients=cell(1,1);
    X=Undeform_vc_xyz(:,1);          %%%%%%%%%%%%% 改2
    Y=Undeform_vc_xyz(:,2);          %%%%%%%%%%%%% 改2
    Z=sagz{i_s}(:,1);
    [N,MATRIX,idx]=ZernikeFit(X,Y,Z,X0,Y0,radius,termOrder,firstTermIndex,kArray);   %%%% 改3
    Fai=MATRIX;
    T_Fai=Fai';
    WW=WW(idx);       %%%%%%%%%  改4
    ZZ=Z(idx);        %%%%%%%%%  改4
    for j=1:length(kArray)
        for i=1:N
            W_Fai(i,:)=WW(i,1)* Fai(i,:);            %%%%%%%% 改5
            bb_rm(j,1)=WW(i,1)*ZZ(i,1)*Fai(i,j)+bb_rm(j,1);     %%%%%%%  改5
        end
    end
      AA_rm=T_Fai*W_Fai;
    coefficients = Least_square ('Pseudoinverse Matrix',AA_rm,bb_rm);     % Select different method to slove the least-square problem
    ZZ=Fai*coefficients;
    Zernike_PV =(max(ZZ(:,1))-min(ZZ(:,1)))*1000000;
    Zernike_RMS=std(ZZ(:,1))*1000000; 
    
    coe_name=strcat(surf_name,'_',subcase_name{i_s});                     % save zernike coffients data into xx.dat files
    res_rpt_2 = [coe_name,'_coe.dat'];
    fid_res_2 = fopen(res_rpt_2,'wt');
    fprintf(fid_res_2, '%s\n',coe_name);
    fprintf(fid_res_2, 'Zernike Term Order : %s\n',termOrder);
    
%     fprintf(fid_res_1, 'N = %d\n',subcase_name{i_s}, N);
    fprintf(fid_res_2, 'PV_%s  = %10.3enm\n',subcase_name{i_s},Zernike_PV);
    fprintf(fid_res_2, 'rms_%s = %10.3enm\n',subcase_name{i_s}, Zernike_RMS);
    Zer_Coefficients=[];
    for i=1:length(kArray)
        fprintf(fid_res_2, '第 %d 阶 = %10.3e\n',i,coefficients(i,1));
        Zer_Coefficients(i,1)=i;
        Zer_Coefficients(i,2)=coefficients(i,1);
    end
        Zer_Coefficients(1,3)=Zernike_PV/1000000;
         Zer_Coefficients(1,4)=Zernike_RMS/1000000;
    
%     for i=1:length(Dis)
%         fprintf(fid_res_2, 'PV_%s  = %10.3f \n',kArray,coefficients{i_s});
%     end
    fclose(fid_res_2);
        %% Zer_Coefficient to Excel %%
        % Time=0:0.01:Num_Time;
        % Time=num2cell(Time');
        Zer_Coffeicients2=num2cell( Zer_Coefficients);
        Title={'Order','Coffeicient(mm)','RMS(mm)','PV(mm)'};
        % Title={'Time(s)','Tx(mm)','Ty(mm)','Tz(mm)','Rx(rad)','Ry(rad)','Rz(rad)'};
        Data_Zer_Coffeicients=cell(length(Zer_Coefficients(:,1))+1,4);
        Data_Zer_Coffeicients(1,:)=Title;
        % Data_RBD(2:end,1)=Time(1:end,1);
        Data_Zer_Coffeicients(2:end,1)=Zer_Coffeicients2(1:end,1);
        Data_Zer_Coffeicients(2:end,2)=Zer_Coffeicients2(1:end,2);
        Data_Zer_Coffeicients(2,3)=Zer_Coffeicients2(1,4);
        Data_Zer_Coffeicients(2,4)=Zer_Coffeicients2(1,3);   
        filename2 = 'Zer_Coffeicients.xlsx';
        xlswrite(filename2,Data_Zer_Coffeicients);
        res_rpt_5=filename2;
%%  ******************************  %%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Print the figure of Sag Change in Zernike and FEM  %%%%%%%%%
    F_Zernike = scatteredInterpolant(X(:,1), Y(:,1), Z(:,1));      %%%%%%三维插值得到拟合曲面函数
    ti_x_Zernike = linspace(min(X(:,1)-X0), max(X(:,1)-X0), radius*2);         
    ti_y_Zernike = linspace(min(Y(:,1)-Y0), max(Y(:,1)-Y0), radius*2);
    [qx_Zernike, qy_Zernike] = meshgrid(ti_x_Zernike,ti_y_Zernike);   %%%%%%%%%设置横纵坐标范围
    input=nan(radius*2);
    input(qx_Zernike.^2+qy_Zernike.^2<=radius^2)=1;            %%%%%%%%建立圆投影图
    qz_Zernike = F_Zernike(qx_Zernike,qy_Zernike);      %%%%%%%%%%%%%求在横纵坐标下矢高z坐标
    qz_Zernike_min=min(min(qz_Zernike));
    qz_Zernike_max=max(max(qz_Zernike));
    d_qz_Zernike=qz_Zernike_max-qz_Zernike_min;
    colormap jet;                                       %%%%%%%%%%%%%设置图形显示形式
    Qz_Zernike=input.*qz_Zernike;                       %%%%%%%%%%%%%计算矢高z的坐标
    surf(ti_x_Zernike,ti_y_Zernike,Qz_Zernike);
    NEW_qz_Zernike_min=qz_Zernike_min+d_qz_Zernike*0.1;
    NEW_qz_Zernike_max=qz_Zernike_max-d_qz_Zernike*0.1;
    caxis([NEW_qz_Zernike_min,NEW_qz_Zernike_max]);     %%%%%%%%%%%%%设置颜色标尺显示划分
    view(0,90);                                 %%%%%%%%%%%显示俯视图
    colorbar('vert');                           %%%%%%%%%%%颜色标尺垂直显示
    shading interp;                             %%%%%%%%%%%%去除网格
    saveas(gcf, 'mirror.jpg');
    saveas(gcf, 'mirror.fig');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%  Output opto-mechanical interface file  '.zpl' %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Create and open the '.zpl' file %%%%%%%%%%%%%%%%%%%%
    zer_name=strcat(surf_name,'_',subcase_name{i_s},'_zer');
    res_rpt_3 = [surf_name,'_zmx.zpl'];
    fid_res_3=fopen(res_rpt_3,'wt');
    
%%%%%%%%%%%%%%%%%%%%%%%%   Output the header file    %%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid_res_3,'! %s\n',zer_name);
    fprintf(fid_res_3,'! Optic # = 2    Optic Label=MIRROR\n');  %%%%%%后期可更改为可交互式,例如设定Optic镜面编号，设定Optic Label 名称
    fprintf(fid_res_3,'! RB motion expressed in optics units=mm and degrees\n');
    fprintf(fid_res_3,'! Normalization radii expressed in optics units=mm\n');
    fprintf(fid_res_3,'! Polynomials normalized & numbered by ZEMAX convention\n');
    fprintf(fid_res_3,'! Rigid-body Z motion represented as THIC command\n');
    fprintf(fid_res_3,'! ************************************************************************\n');

%%%%%%%%%%%%%  Output rigid body displacement, specular number and Zernike coefficient  %%%%%%%%%%%%%%  
    fprintf(fid_res_3,'dz = %10.15e\n',rigid_motion_NU(6,1));
    fprintf(fid_res_3,'SURP 002, BOR,            0\n');         %%%%%%顶点曲率变化视为不变，即为0；镜面编号视为'013'，后期可根据镜面编号修正
    fprintf(fid_res_3,'SURP 002, BDX, %10.15e\n',rigid_motion_NU(4,1));
    fprintf(fid_res_3,'SURP 002, BDY, %10.15e\n',rigid_motion_NU(5,1));
    fprintf(fid_res_3,'SURP 002, BTX, %10.15e\n',rigid_motion_NU(1,1)*180/pi);
    fprintf(fid_res_3,'SURP 002, BTY, %10.15e\n',rigid_motion_NU(2,1)*180/pi);
    fprintf(fid_res_3,'SURP 002, BTZ, %10.15e\n',rigid_motion_NU(3,1)*180/pi);
    fprintf(fid_res_3,'SURP 002, APU,            2\n');
    fprintf(fid_res_3,'SURP 001, THIC, THIC(001)+dz\n');
    fprintf(fid_res_3,'SURP 002, THIC, THIC(002)-dz\n');
    switch termOrder
        case 'StandardOrder'            
             fprintf(fid_res_3,'B$="SZERNSAG"\n');
        case 'FringeOrder' 
             fprintf(fid_res_3,'B$="FZERNSAG"\n');
        case 'NStandarOrder'
             fprintf(fid_res_3,'B$="SZERNSAG"\n');
    end
    fprintf(fid_res_3,'SURP 002,TYPE,B$\n');
    fprintf(fid_res_3,'SURP 002,PARM,1,0\n');
     fprintf(fid_res_3,'SURP 002,PARM,%10.15e, 9\n',X0); %离轴镜面中心
       fprintf(fid_res_3,'SURP 002,PARM,%10.15e, 10\n',Y0);
    fprintf(fid_res_3,'SURP 002,EDVA,   %d,  1\n',kmax);
    fprintf(fid_res_3,'SURP 002,EDVA, %10.15e,  2\n',radius);     %%%%%%%%radius代指最大镜面半径
    for i=1:36
        fprintf(fid_res_3, 'SURP 002,EDVA, %10.15e,    %d\n',coefficients(i,1),i+2);
    end
    fclose(fid_res_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
%%%%%%%%%%%%  Save the results files to the results folder   %%%%%%%%%%
mkdir result;
sourcePath=pwd;
targetPath=strcat(sourcePath,'\result');
movefile([sourcePath,'\',res_rpt_1],targetPath);
movefile([sourcePath,'\',res_rpt_2],targetPath);
movefile([sourcePath,'\',res_rpt_3],targetPath);
movefile([sourcePath,'\',res_rpt_4],targetPath);
movefile([sourcePath,'\',res_rpt_5],targetPath);
movefile([sourcePath,'\','mirror.jpg'],targetPath);
movefile([sourcePath,'\','mirror.fig'],targetPath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  end   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Called Function——Zernike Fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,MATRIX,idx]=ZernikeFit(X,Y,Z,X0,Y0,radius,termOrder,firstTermIndex,kArray,normalQ)       %%%%%% 改6
% Usage: coefficients=ZernikeFit(X,Y,Z,radius,termOrder,firstTermIndex,kArray,normalQ(option))
%   coefficients: coefficients of Zernike fit
%   X: X
%   Y: Y
%   Z: Z
%   X0: x value of mirror center
%   Y0: y value of mirror center
%   radius: radius (>=mirror aperture)
%   termOrder: term order ('StandardOrder','ParityOrder','FringeOrder')
%   firstTermIndex: first term index (0,1,...)  zernike起始项，一般为1
%   kArray: array of single index k:一个多项式的序列编号矩阵，选出参与拟合的多项；，通常firstTermIndex:kmax
%   normalQ(option): use normal Zernike (true,false(default))
%   for example:
error(nargchk(5,12,nargin,'struct'));
if nargin==9
    normalQ=false;
end
SX=(X-X0)/radius;                               %直角坐标系下，坐标归一化；
SY=(Y-Y0)/radius; 
[THETA,R]=cart2pol(SX,SY);                %cart2pol(x,y)将笛卡尔坐标转换为对应的极坐标（theta,rho）
idx=(R<=1)&(~isnan(Z));                     %判断是否归一化并且轴向位移不为空，赋予1  
iMax=length(kArray);                          % iMax是zernike的项数
MATRIX=nan(length(R(idx)),iMax);    %生成length(R(idx))*iMax空矩阵
N=length(R(idx));
for i=1:iMax
    [n,l]=SingleIndexToDoubleIndex(termOrder,firstTermIndex,kArray(i));
    MATRIX(:,i)=Zernike(termOrder,n,l,R(idx),THETA(idx),kArray(i),normalQ);
end
   MATRIX=MATRIX;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n,l]=SingleIndexToDoubleIndex(termOrder,firstTermIndex,k)
% Usage: [n,l]=SingleIndexToDoubleIndex(termOrder,firstTermIndex,k)
%   n: degree part of double index
%   l: frequency part of double index
%   termOrder: term order ('StandardOrder','ParityOrder','FringeOrder')
%   firstTermIndex: first term index (0,1,...)
%   k: single index
% Author: Yasushi Iwasaki (Canon Inc.)
error(nargchk(3,3,nargin,'struct'));

% termOrder='FringeOrder';

switch termOrder
    case 'StandardOrder'         %%%%%% 单位幅值
        j=k-firstTermIndex;      %k=kArray(i);
        n=floor((-1+sqrt(8*j+1))/2);
        l=2*j-n*(n+2);
    case 'ParityOrder'
        n=floor((-1+sqrt(8*j+1))/2);
        i=j-n*(n+1)/2+1;
        m=i-mod(n+i,2);
        l=(-1)^k*m;
    case 'FringeOrder'
        j=k-firstTermIndex;      %k=kArray(i);
        t=floor(sqrt(j));
        u=j-t^2;
        v=floor(u/2);
        n=t+v;
        l=(-1)^u*(t-v);
    case 'NStandarOrder'        %%%% 单位RMS
        j=k;
        n=ceil((sqrt(8*j+1)-3)/2);
        Cri=mod(j-n*(n+3)/2,2);
        if Cri==0
            l=n*(n+1)/2-j;
        else
            l=j-n*(n+1)/2-1;
        end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RET=Zernike(termOrder,n,l,R,THETA,k2,normalQ)
% Usage: RET=Zernike(n,l,R,THETA,normalQ(option))
%   RET: Zernike
%   n: degree part of double index
%   l: frequency part of double index
%   R: radial variable
%   THETA: angular variable
%   normalQ(option): use normal Zernike (true,false(default))
% Author: Yasushi Iwasaki (Canon Inc.)
error(nargchk(6,7,nargin,'struct'));      %检查输入参数的个数在【4,5】之间，否则报错
if nargin==6
    normalQ=false;
end

% termOrder='FringeOrder';

switch termOrder
        case 'StandardOrder'
        if l<=0
            RET=RadialPolynomial(n,l,R).*cos(l*THETA);
        else
            RET=RadialPolynomial(n,l,R).*sin(l*THETA);
        end
    case 'FringeOrder'
        if l>=0
            RET=RadialPolynomial(n,l,R).*cos(l*THETA);
        else
            RET=RadialPolynomial(n,l,R).*sin(-l*THETA);
        end
    case 'NStandarOrder'
        m2=abs(l);
        Cri2=mod(k2,2);
        if l==0
            RET=sqrt(n+1)*RadialPolynomial(n,l,R);
        elseif Cri2==0
            RET=sqrt(2*(n+1))*RadialPolynomial(n,l,R).*cos(m2*THETA);
        else
            RET=sqrt(2*(n+1))*RadialPolynomial(n,l,R).*sin(m2*THETA);
        end
end
if normalQ
    RET=RET/ZernikeRMS(n,l);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RET=RadialPolynomial(n,l,R)
%此函数为Zernike多项式的径向关系递推式
% Usage: RET=RadialPolynomial(n,l,R)
%   RET: radial polynomial
%   n: degree part of double index       径向波数
%   m: degree part of double index      环向波数
%   l: frequency part of double index   
%   R: radial variable
error(nargchk(3,3,nargin,'struct'));
m=abs(l);
p=(n-m)/2;
q=(n+m)/2;
RET=0;
for s=0:p
    RET=RET+(-1)^s*factorial(n-s)/(factorial(s)*factorial(p-s)*factorial(q-s))*R.^(n-2*s);    %zernike多项式径向关系式；factorial()表示阶乘
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Called Function——Read Undeformed mirror face-nodes imformation from xx.rpt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Face_nodes = surface_byfile(surface_type, rpt_path_name, Surf_VC, Surf_c, check_flag, check_tolerence)
% surface setting
% Surf_VC = [[1000; 771.273; 1639.748] [1000; 766.273; 1648.408] [1010; 771.273; 1639.748]];
% Surf_c = -1/2400;
% read .rpt
Face_nodes = dlmread(rpt_path_name);
%  check_flag=1;
%  check_tolerence=0;
%  surface_type='sphere';
if check_flag ==1
    switch surface_type
        case 'sphere'
            %%%%%% Define the coordinate And Calculate the transformation matrix %%%%%%%
            %%% Define the global 'coordinate 'g_P & the vertex coordinate 'v_P' %%%
            g_P=eye(3);    %%%% the global coordinate 'g_P'
            VC_z = Surf_VC(:,2) - Surf_VC(:,1);
            VC_z = VC_z/norm(VC_z);
            VC_y = cross(VC_z, Surf_VC(:,3) - Surf_VC(:,1));
            VC_y = VC_y/norm(VC_y);
            VC_x = cross(VC_y,VC_z);
            VC_x = VC_x/norm(VC_x);
            v_P=[VC_x,VC_y,VC_z];   %%%% the vertex coordinate 'v_P'
            %%%%%%%%%% Calculate the coordinate transformation matrix   %%%%%%%%%%%%%%%
            v_P_g=-Surf_VC(:,1);       %%%%%%%%%%%  平动矩阵-> 'v_P_g'
            c_P=v_P-repmat(v_P_g,1,3);   %%%% 中转坐标系-> 'c_P'
            v_R_g=c_P/g_P;     %%%%% 转动矩阵->‘v_R_g’
            Tranform_gc2vc=[v_R_g,v_P_g;zeros(1,3),eye(1)];   %%%%% 全局转局部矩阵 Tranform_gc2vc
            %%%%%%%%%%  ********************************************    %%%%%%%%%%%%%
            Surf_Center = Tranform_gc2vc\[0; 0; 1/Surf_c; 1]; 
            Surf_Center = Surf_Center(1:3);
            err = abs(sqrt(sum((repmat(Surf_Center,1,size(Face_nodes,1)) - transpose(Face_nodes(:, 2:4))).^2))) - abs(1/Surf_c);    %筛选剔除误差比较大的离散点（怀疑这种容差判断的意义）
            Face_nodes = Face_nodes(transpose(err < check_tolerence), :);
         case 'paraboloid'
            %%%%%% Define the coordinate And Calculate the transformation matrix %%%%%%%
            %%% Define the global 'coordinate 'g_P & the vertex coordinate 'v_P' %%%
            g_P=eye(3);    %%%% the global coordinate 'g_P'
            VC_z = Surf_VC(:,2) - Surf_VC(:,1);
            VC_z = VC_z/norm(VC_z);
            VC_y = cross(VC_z, Surf_VC(:,3) - Surf_VC(:,1));
            VC_y = VC_y/norm(VC_y);
            VC_x = cross(VC_y,VC_z);
            VC_x = VC_x/norm(VC_x);
            v_P=[VC_x,VC_y,VC_z];   %%%% the vertex coordinate 'v_P'
            %%%%%%%%%% Calculate the coordinate transformation matrix   %%%%%%%%%%%%%%%
            v_P_g=-Surf_VC(:,1);       %%%%%%%%%%%  平动矩阵-> 'v_P_g'
            c_P=v_P-repmat(v_P_g,1,3);   %%%% 中转坐标系-> 'c_P'
            v_R_g=c_P/g_P;     %%%%% 转动矩阵->‘v_R_g’
            Tranform_gc2vc=[v_R_g,v_P_g;zeros(1,3),eye(1)];   %%%%% 全局转局部矩阵 Tranform_gc2vc
            %%%%%%%%%%  ********************************************    %%%%%%%%%%%%%
            Surf_Center = Tranform_gc2vc\[0; 0; 1/Surf_c; 1];
            Surf_Center = Surf_Center(1:3);
            err = abs(sqrt(sum((repmat(Surf_Center,1,size(Face_nodes,1)) - transpose(Face_nodes(:, 2:4))).^2))) - abs(1/Surf_c);
            Face_nodes = Face_nodes(transpose(err < check_tolerence), :);
        case 'plane'
            Plane_norm = Surf_VC(:,2)-Surf_VC(:,1);
            Plane_norm = Plane_norm/norm(Plane_norm);
            err = abs(sum((repmat(transpose(Surf_VC(:,1)), size(Face_nodes,1), 1) - Face_nodes(:, 2:4)).*repmat(transpose(Plane_norm), size(Face_nodes,1), 1),2));
            Face_nodes = Face_nodes(transpose(err < check_tolerence), :);
        otherwise
            error('something wrong with surface_type')     
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Called Function——Read Area factor from Optistuct output files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Area_factors =Area_factor_byfiles(AreaFactor_pch_name,Face_nodes)
fid_punch = fopen(AreaFactor_pch_name);
i_s = 1;                               % number of subcase
subcase_name = cell(1,1);
Dis = cell(1,1);                   % cell of displacements
while ~feof(fid_punch)     % 文件指针达到文件末尾时，表达式值为“假”，否则为“真”
    tline = fgetl(fid_punch);
    if length(tline)>6 && strcmp(tline(1:6) , '$LABEL')
        temp = tline(12:end);
        temp(isspace(temp)) = [];
        subcase_name{i_s} = temp;
    elseif length(tline)>4 && strcmp(tline(1:5) , '$SPCF')
        tline = fgetl(fid_punch);                       % Don`t delete
        tline = fgetl(fid_punch);
        if length(tline)>11 && strcmp(tline(1:11) , '$SUBCASE ID')
            i = 1;
            Dis_i_s = zeros(5,1);
            while i <= length(Face_nodes)
                tline = fgetl(fid_punch);
                gridid = str2double(tline(1:10));
                if ~isnan(gridid)
                    if ismember(gridid, Face_nodes(:,1))                  % 判断gridid是否在Face_nodes(:,1)集合中出现
                        Dis_i_s(i, 1) = gridid;
%                         Dis_i_s(i, 2) = str2double(tline(21:40));
%                         Dis_i_s(i, 3) = str2double(tline(41:56));
                        Dis_i_s(i, 4) = str2double(tline(57:72));
                        i = i + 1;
                    end
                end
            end
            Dis{i_s} = Dis_i_s;
            i_s = i_s + 1;
        end
    end
end
fclose(fid_punch);

% verify if Dis_i_s compatible with Face_nodes
% for i = 1:length(Dis)
%     if sum(abs(Dis{i}(:,1) - Face_nodes(:,1))) > 0
%         error('error in file reading')
%     end
% end

Area_factors=Dis_i_s;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Called Function——The Over-determined equation solution
function Sol= Least_square (cal_method,A,b)
% A=A_rm;
% b=b_rm;
% Sol=cell(1,1);
% cal_method='Pseudoinverse Matrix';

switch cal_method
  
    case 'Pseudoinverse Matrix'        % 广义逆矩阵法（Pseudoinverse Matrix）
%           for i=1:length(Dis)   
%             Sol=A{i}\b{i};
%           end
            Sol=A\b;
    case 'QR_MGSO'                           % Modified Gram-Schmidt Orthogonalization(CGSO)
            [m,n]=size(A);                           % m是节点个数；n是Zernike的项数    
            Q=zeros(m,n);
            R=zeros(n,n);
            for j=1:n
                V=A(:,j);                               
                for i=1:j-1                                     
                    R(i,j)=Q(:,i)'*A(:,j);          
                    V=V-R(i,j)*Q(:,i);
                end
                R(j,j)=norm(V);
                Q(:,j)=V/R(j,j);
            end
            Sol=R\(Q'*b);     
    case 'QR_Householder'                 % QR decomposition using Householder matrix
            [m,n]=size(A);
            H=eye(m);
            for i=1:n
                s=length(A(i:m,i));
                e=[1;zeros(s-1,1)];
                if A(i:m,i)==norm(A(i:m,i))*e      %x,z same direction
                    HH=eye(s)-2.*(e*e');
                else
                    t=(A(i:m,i)-norm(A(i:m,i))*e)/norm(A(i:m,i)-norm(A(i:m,i))*e);
                    HH=eye(s)-2.*(t*t');
                end
            end
                H_temp=blkdiag(eye(i-1),HH);
                H=H_temp*H;
                A=H_temp*A;
                Q=H';
                R=A;
                Sol=R\(Q'*b);   
    case 'QR_GivenRotations'            % QR decomposition using Givens matrix
            [m,n]=size(A);
            G=eye(m);
            for i=1:1:n
              u=A(i:m,i);
              s=length(u);
              GG=eye(s,s);
              sum=0;
              for j=2:1:s
                  if u(j) ~=0
                      if sum==0
                          sum=u(1)*u(1);
                          sumNew=sum+u(j)*u(j);
                          cos=u(1)/sqrt(sumNew);
                          sin=u(j)/sqrt(sumNew);
                      else
                          sumNew=sum+u(j)*u(j);
                          cos=sqrt(sum)/sqrt(sumNew);
                          sin=u(j)/sqrt(sumNew);
                      end
                      sum=sumNew;
                      G_temp=eye(s,s);
                      G_temp(1,1)=cos;
                      G_temp(1,j)=sin;
                      G_temp(j,1)=-sin;
                      G_temp(j,j)=cos;
                      GG=G_temp*GG;
                  end
              end
              G_temp=blkdiag(eye(i-1),GG);
              G=G_temp*G;
              A=G_temp*A;
           end
           Q=G';
           R=A;
           Sol=R\(Q'*b);   
    case 'SVD'                                   % Singular Value Decomposition
            [U,S,V]=svd(A,0);
            h=U'*b;
            k=S\h;
            Sol=V*k;                
    otherwise
        error('something wrong with surface_type');
end
end



%% Called Function——The Quatstiting_xy equation solution
function [Tx,Ty,Tz,Rx,Ry,Rz] = Quatstiting_xy( x1,y1,z1,x2,y2,z2)         
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明

% %% Initialization of parameters %%
% rpt_path_name1='Flat1.rpt';rpt_path_name2='Flat_RXYZ_Large.rpt';surf_name='RBD-PM';Surf_VC=[0 0 0;0 0 1;1 0 0]';Surf_c=1500;surfacetype='paraboloid';
% %% ****************************** %%
% 
% %% Read the data of the Mirror nodes %%
% Face_nodes1 = surface_byfile(surfacetype, rpt_path_name1, Surf_VC, Surf_c, 0, 1); 
% Face_nodes2 = surface_byfile(surfacetype, rpt_path_name2, Surf_VC, Surf_c, 0, 1); 
% x1=Face_nodes1(:,2);
% y1=Face_nodes1(:,3);
% z1=Face_nodes1(:,4);
% x2=Face_nodes2(:,2);
% y2=Face_nodes2(:,3);
% z2=Face_nodes2(:,4);
% %% ******************************  %%

n=size(x2,1);
% n=length(x2);

x_1=sum(x1(:))/n;
x_2=sum(x2(:))/n;
y_1=sum(y1(:))/n;
y_2=sum(y2(:))/n;
z_1=sum(z1(:))/n;
z_2=sum(z2(:))/n;
x11=x1-x_1;
x21=x2-x_2;
y11=y1-y_1;
y21=y2-y_2; 
z11=z1-z_1;
z21=z2-z_2;

X1=x11.^2;
Y1=y11.^2;
Z1=z11.^2;
X2=x21.^2;
Y2=y21.^2;
Z2=z21.^2;

% z_B=z_1.^2;     %% 为什么不对z_i做类似Xi和Yi的操作？
% z_S=z_2.^2;
% L_B=sum(X1(:))+sum(Y1(:))+sum(z_B(:));
% L_L=sum(X2(:))+sum(Y2(:))+sum(z_S(:));
% Lamd=sqrt(L_B/L_L);  %% 缩放比例
% % Lamd=1;
L_B=sum(X1(:))+sum(Y1(:))+sum(Z1(:));
L_L=sum(X2(:))+sum(Y2(:))+sum(Z2(:));
Lamd=sqrt(L_B/L_L);  %% 缩放比例


Xs=x2;Xb=x1;Ys=y2;Yb=y1;
xxsb=Xs.*Xb;yysb=Ys.*Yb;zz=z1.*z2;
xybs=Xb.*Ys;xysb=Xs.*Yb;
xz_sb=Xs.*z1;xz_bs=Xb.*z2;
yz_sb=Ys.*z1;yz_bs=Yb.*z2;

C1=-2*Lamd*[sum(xxsb)+sum(yysb)+sum(zz)                sum(yz_bs)-sum(yz_sb)                  -sum(xz_bs)+sum(xz_sb)                         sum(xybs)-sum(xysb);
             sum(yz_bs)-sum(yz_sb)                     sum(xxsb)-sum(yysb)-sum(zz)           sum(xybs)+sum(xysb)                             sum(xz_bs)+sum(xz_sb);
             -sum(xz_bs)+sum(xz_sb)                    sum(xybs)+sum(xysb)                   -sum(xxsb)+sum(yysb)-sum(zz)                    sum(yz_bs)+sum(yz_sb);
             sum(xybs)-sum(xysb)                       sum(xz_bs)+sum(xz_sb)                 sum(yz_bs)+sum(yz_sb)                            -sum(xxsb)-sum(yysb)+sum(zz);];
    
    
C2=2*[0                             sum(-Lamd*Xb)+sum(Xs)           sum(-Lamd*Yb)+sum(Ys)       sum(-Lamd*z1)+sum(z2);
      sum(Lamd*Xb)-sum(Xs)            0                               sum(Lamd*z1)+sum(z2)      sum(-Lamd*Yb)-sum(Ys);
      sum(Lamd*Yb)-sum(Ys)          sum(-Lamd*z1)-sum(z2)         0                           sum(Lamd*Xb)+sum(Xs);
      sum(Lamd*z1)-sum(z2)        sum(Lamd*Yb)+sum(Ys)              sum(-Lamd*Xb)-sum(Xs)         0];
       
  
  
  N=((C2'*C2)./(2*n)-C1-C1')./2;
  [x,y]=eig(N);                 %%%特征值和特征向量,x为特征向量，y为特征值
%   y
%     a=(find(max(max(y))));           %%% 此处理会找到两个特征值，采用max(max())的计算及结果有较大差异；
    [a b]=find(y==max(max(y)))    
    p=x(:,b);                     %%% 最大特征值对应的特征向量，两个特征值对应两个特征向量咋处理？
%   p
  q=-C2*p./2./n;                %%% qd  求解依据公式是啥？
  M=Quat2Mat_xy(p,q);           %%% 求解齐次变换矩阵或旋转矩阵
  MR=M;                         %%% 输出旋转矩阵
 
% %%  ‘M=Quat2Mat_xy(p,q)'求解的为旋转矩阵时 %%
%     x1=x1';
%     y1=y1';
%     z1=z1';
%     x2=x2';
%     y2=y2';
%     z2=z2';
%   REAL=[x1; y1; z1];
%   SI=[x2;y2;z2];
%   T= mean(( SI - M*REAL),2);       %%求所有点平移量的平均值
%   AddnanL=zeros(1,3);
% %   T(1:2)=0;                      %% 之前存在的意义是什么？
%   M=[M;AddnanL];
%   T=[T;1];
%   M=[M,T];
% %%  ******************************************  %%

%% 欧拉角Rx、Ry和Rz求解，携带Tx、Ty、Tz  %%
Rx=atan(2*(p(1)*p(2)+p(3)*p(4))/(1-2*(p(2)^2+p(3)^2)));
Ry=asin(2*(p(1)*p(3)-p(4)*p(2)));
Rz=atan(2*(p(1)*p(4)+p(2)*p(3))/(1-2*(p(3)^2+p(4)^2)));
Tx=M(1,4)/2;
Ty=M(2,4)/2;
Tz=M(3,4)/2; 
%% ****************************************** %%

% %% Rigid body displacement output to Excel %%
% Title={'Tx(mm)','Ty(mm)','Tz(mm)','Rx(rad)','Ry(rad)','Rz(rad)'};
% Data_RBD=cell(2,6);
% Data_RBD(1,:)=Title;
% Data_RBD(2,1)=num2cell(Tx);
% Data_RBD(2,2)=num2cell(Ty);
% Data_RBD(2,3)=num2cell(Tz);
% Data_RBD(2,4)=num2cell(Rx);
% Data_RBD(2,5)=num2cell(Ry);
% Data_RBD(2,6)=num2cell(Rz);
% filename = 'RBD.xlsx';
% xlswrite(filename,Data_RBD);
% %%  ******************************  %%
end

function Cnb = Quat2Mat_xy(qnb, p) 
%姿态四元数转化为姿态矩阵  

% Qq=[qnb(1) -qnb(2) -qnb(3) -qnb(4);
%     qnb(2) qnb(1) -qnb(4) qnb(3);
%     qnb(3) qnb(4) qnb(1) -qnb(2);
%     qnb(4) -qnb(3) qnb(2) qnb(1)];
% Q_q=[qnb(1) -qnb(2) -qnb(3) -qnb(4);
%     qnb(2) qnb(1) qnb(4) -qnb(3);
%     qnb(3) -qnb(4) qnb(1) qnb(2);
%     qnb(4) qnb(3) -qnb(2) qnb(1)];
% R=Q_q'*Qq;
% T=2*Q_q'*p;
% R(:,4)=T;
%     Cnb = R;
%     B=qnb(2:4);
%     E=[1 0 0;0 1 0;0 0 1];
% Cnb =(qnb(1)^2-B'*B)*E+2*(B'*B)+qnb(1)*[0 -qnb(4) qnb(3);qnb(4) 0 -qnb(2);-qnb(3) qnb(2) 0];

% 
    q11 = qnb(1)*qnb(1); q12 = qnb(1)*qnb(2); q13 = qnb(1)*qnb(3); q14 = qnb(1)*qnb(4); 
    q22 = qnb(2)*qnb(2); q23 = qnb(2)*qnb(3); q24 = qnb(2)*qnb(4);     
    q33 = qnb(3)*qnb(3); q34 = qnb(3)*qnb(4);  
    q44 = qnb(4)*qnb(4);
% 
        qp01=qnb(1)*p(2); qp10=qnb(2)*p(1); qp23=qnb(3)*p(4); qp32=qnb(4)*p(3);
        qp02=qnb(1)*p(3); qp13=qnb(2)*p(4); qp20=qnb(3)*p(1); qp31=qnb(4)*p(2);
        qp03=qnb(1)*p(4); qp12=qnb(2)*p(3); qp21=qnb(3)*p(2); qp30=qnb(4)*p(1);
%     

%                             齐次变换矩阵形式                                          %%
    Cnb = [ q11+q22-q33-q44,  2*(q23-q14),     2*(q24+q13)        2*(qp01-qp10+qp23-qp32);
            2*(q23+q14),      q11-q22+q33-q44, 2*(q34-q12)        2*(qp02-qp13-qp20+qp31);
            2*(q24-q13),      2*(q34+q12),     q11-q22-q33+q44    2*(qp03+qp12-qp21-qp30);
            0                  0                 0                     1];
% ************************************************************************************  %%

% %%                              仅旋转矩阵                                         %%            
% %     Cnb = [ q11+q22-q33-q44,  2*(q23-q14),     2*(q24+q13)        ;
% %             2*(q23+q14),      q11-q22+q33-q44, 2*(q34-q12)        ;
% %             2*(q24-q13),      2*(q34+q12),     q11-q22-q33+q44    ];
%             Cnb = [ q11+q22-q33-q44,  2*(q23-q14),     2*(q24+q13)        ;
%             2*(q14+q23),      q11-q22+q33-q44, 2*(q34-q12)        ;
%             2*(q24-q13),      2*(q34+q12),     q11-q22-q33+q44    ];
% %% ******************************************************************************  %%
end


