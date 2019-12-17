t_grid=linspace(0,dt*Nmax,Nmax);
figure(1)
plot(t_grid,abs(dx_),t_grid,abs(dy_),t_grid,abs(dz_),t_grid,dr,t_grid,md_r);
legend('dx','dy','dz','dr','3\sigma');
title('relative position error LVLH');
xlabel('time [s]');
ylabel('error [m]');

figure(2)
plot(t_grid,abs(dvx_),t_grid,abs(dvy_),t_grid,abs(dvz_),t_grid,dvr,t_grid,md_vr);
legend('dvx','dvy','dvz','dv','3\sigma');
title('relative velocity error LVLH');
xlabel('time [s]');
ylabel('error [m/s]');

figure(3)
plot(t_grid,abs(dwx_),t_grid,abs(dwy_),t_grid,abs(dwz_),t_grid,dwr,t_grid,md_wr);
legend('dwx','dwy','dwz','dw','3\sigma');
title('relative angular velocity error');
xlabel('time [s]');
ylabel('error [rad/s]');

figure(4)
plot(t_grid,abs(dk_1),t_grid,abs(dk_2),t_grid,dk,t_grid,md_kr);
legend('dk1','dk2','dk','3\sigma');
title('inertial parameters error');
xlabel('time [s]');
ylabel('error []');

figure(5)
plot(t_grid,err);
title('relative orientation error');
xlabel('time [s]');
ylabel('error [rad]');

% figure(6)
% for i=1:length(landmarks_map.Sk)
%     S=landmarks_map.Sk{i};
%     scatter3(S(:,1),S(:,2),S(:,3));
%     hold on;
% end

% %%
% figure(1)
% Nmax=300;
% t_grid=linspace(0,dt*Nmax,Nmax);
% plot(t_grid,abs(dx_),t_grid,abs(dy_),t_grid,abs(dz_),t_grid,dr);
