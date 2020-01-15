function dh_dxin = jacobian_meas_n_points_n(Xin,Yin,Zin,alpha_u,alpha_v,b)
%JACOBIAN_MEAS_N_POINTS_N
%    DH_DXIN = JACOBIAN_MEAS_N_POINTS_N(XIN,YIN,ZIN,ALPHA_U,ALPHA_V,B)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Jan-2020 14:22:12

t2 = sparse(1.0./Zin);
t3 = sparse(t2.^2);
dh_dxin = sparse([1,2,1,2,3],[1,2,3,3,3],[-alpha_u.*t2,-alpha_v.*t2,Xin.*alpha_u.*t3,Yin.*alpha_v.*t3,-alpha_u.*b.*t3],3,3);
