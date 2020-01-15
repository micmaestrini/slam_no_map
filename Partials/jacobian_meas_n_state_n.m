function dh_dxn = jacobian_meas_n_state_n(X0i,Xin,Y0i,Yin,Z0i,Zin,alpha_u,alpha_v,b,fn,scn1,scn2,scn3,sn1,sn2,sn3)
%JACOBIAN_MEAS_N_STATE_N
%    DH_DXN = JACOBIAN_MEAS_N_STATE_N(X0I,XIN,Y0I,YIN,Z0I,ZIN,ALPHA_U,ALPHA_V,B,FN,SCN1,SCN2,SCN3,SN1,SN2,SN3)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Jan-2020 14:22:11

t2 = sparse(cos(fn));
t3 = sparse(sin(fn));
t4 = sparse(X0i.*sn1);
t5 = sparse(Y0i.*sn2);
t6 = sparse(Z0i.*sn3);
t7 = sparse(scn1.^2);
t8 = sparse(scn1.^3);
t9 = sparse(scn2.^2);
t10 = sparse(scn2.^3);
t11 = sparse(scn3.^2);
t12 = sparse(scn3.^3);
t13 = sparse(sn1.^2);
t14 = sparse(sn1.^3);
t15 = sparse(sn2.^2);
t17 = sparse(sn2.^3);
t18 = sparse(sn3.^2);
t20 = sparse(sn3.^3);
t22 = sparse(X0i.*sn3.*2.0);
t23 = sparse(Y0i.*sn3.*2.0);
t24 = sparse(Z0i.*sn1.*4.0);
t25 = sparse(Z0i.*sn2.*4.0);
t26 = sparse(scn1.*scn2.*2.0);
t27 = sparse(scn1.*scn3.*2.0);
t28 = sparse(scn2.*scn3.*2.0);
t29 = sparse(-X0i);
t30 = sparse(-Y0i);
t31 = sparse(-Z0i);
t32 = sparse(1.0./Zin);
t34 = sparse(-scn1);
t35 = sparse(-scn2);
t36 = sparse(-scn3);
t53 = sparse(X0i.*sn2.*sn3.*3.0);
t55 = sparse(Y0i.*sn1.*sn3.*3.0);
t57 = sparse(Z0i.*sn1.*sn2.*3.0);
t16 = sparse(t13.^2);
t19 = sparse(t15.^2);
t21 = sparse(t18.^2);
t33 = sparse(t32.^2);
t37 = sparse(-t22);
t38 = sparse(-t25);
t39 = sparse(-t26);
t40 = sparse(-t27);
t41 = sparse(t4.*t13);
t42 = sparse(t4.*t14);
t46 = sparse(t5.*t15);
t47 = sparse(t5.*t17);
t51 = sparse(t6.*t18);
t52 = sparse(t6.*t20);
t54 = sparse(sn2.*t4.*6.0);
t56 = sparse(sn1.*t5.*6.0);
t58 = sparse(scn1.*t9);
t59 = sparse(scn2.*t7);
t60 = sparse(scn1.*t11);
t61 = sparse(scn3.*t7);
t62 = sparse(scn2.*t11);
t63 = sparse(scn3.*t9);
t64 = sparse(t7.*8.0);
t65 = sparse(t9.*8.0);
t66 = sparse(t11.*8.0);
t67 = sparse(X0i.*t20.*2.0);
t68 = sparse(X0i.*t15.*6.0);
t69 = sparse(Y0i.*t20.*2.0);
t70 = sparse(Y0i.*t13.*6.0);
t71 = sparse(Z0i.*t14.*4.0);
t72 = sparse(Z0i.*t17.*4.0);
t74 = sparse(-t55);
t76 = sparse(t4.*t15);
t77 = sparse(X0i.*sn2.*t20);
t78 = sparse(X0i.*sn3.*t17);
t79 = sparse(t5.*t13);
t80 = sparse(Y0i.*sn1.*t20);
t81 = sparse(Y0i.*sn3.*t14);
t82 = sparse(Z0i.*sn1.*t17);
t83 = sparse(Z0i.*sn2.*t14);
t84 = sparse(sn2.*sn3.*t4.*8.0);
t85 = sparse(sn1.*sn3.*t5.*8.0);
t86 = sparse(sn1.*sn2.*t6.*8.0);
t87 = sparse(sn1.*sn2.*sn3.*t4);
t88 = sparse(sn1.*sn2.*sn3.*t5);
t89 = sparse(sn1.*sn2.*sn3.*t6);
t94 = sparse(t4.*t17.*2.0);
t96 = sparse(t4.*t18.*3.0);
t97 = sparse(t15.*t22);
t98 = sparse(sn1.*sn3.*t4.*6.0);
t100 = sparse(t13.*t23);
t101 = sparse(t5.*t14.*2.0);
t102 = sparse(t5.*t18.*3.0);
t103 = sparse(sn2.*sn3.*t5.*6.0);
t104 = sparse(t15.*t24);
t105 = sparse(t13.*t25);
t106 = sparse(sn1.*sn3.*t6.*4.0);
t107 = sparse(t6.*t13.*4.0);
t108 = sparse(sn2.*sn3.*t6.*4.0);
t109 = sparse(t6.*t15.*4.0);
t111 = sparse(sn2.*t4.*t18.*2.0);
t112 = sparse(sn1.*t5.*t18.*2.0);
t114 = sparse(sn2.*t20.*t29);
t115 = sparse(X0i.*sn3.*t15.*-2.0);
t116 = sparse(sn3.*t17.*t29);
t119 = sparse(Z0i.*sn1.*t15.*-4.0);
t121 = sparse(sn1.*t4.*t18.*2.0);
t122 = sparse(sn2.*t5.*t18.*2.0);
t124 = sparse(t7+t9+t11+1.0);
t125 = sparse(t13+t15+t18+1.0);
t43 = sparse(X0i.*t19);
t44 = sparse(X0i.*t21);
t45 = sparse(Y0i.*t16);
t48 = sparse(Y0i.*t21);
t49 = sparse(Z0i.*t16);
t50 = sparse(Z0i.*t19);
t73 = sparse(-t54);
t75 = sparse(-t56);
t90 = sparse(t19.*t29);
t91 = sparse(-t67);
t92 = sparse(t16.*t30);
t93 = sparse(-t71);
t95 = sparse(sn2.*t41.*2.0);
t99 = sparse(sn1.*t46.*2.0);
t110 = sparse(-t84);
t113 = sparse(-t96);
t117 = sparse(-t102);
t118 = sparse(-t103);
t120 = sparse(-t108);
t123 = sparse(-t87);
t126 = sparse(t64+t65);
t127 = sparse(t64+t66);
t128 = sparse(t65+t66);
t129 = sparse(1.0./t124.^2);
t130 = sparse(1.0./t125.^3);
t131 = sparse(t8+t28+t34+t58+t60);
t132 = sparse(t12+t26+t36+t61+t63);
t133 = sparse(t10+t35+t40+t59+t62);
t134 = sparse(t12+t36+t39+t61+t63);
t135 = sparse(t126.*t129);
t136 = sparse(t127.*t129);
t137 = sparse(t128.*t129);
t141 = sparse(t3.*t129.*t131.*4.0);
t142 = sparse(t2.*t129.*t131.*4.0);
t143 = sparse(t3.*t129.*t133.*4.0);
t144 = sparse(t2.*t129.*t133.*4.0);
t148 = sparse(t4+t5+t41+t46+t53+t74+t76+t79+t80+t81+t88+t107+t109+t113+t114+t116+t117+t123);
t149 = sparse(t23+t29+t38+t42+t44+t68+t69+t72+t75+t90+t99+t100+t101+t105+t110+t112+t118+t120+t121);
t150 = sparse(t24+t30+t37+t47+t48+t70+t73+t85+t91+t92+t93+t94+t95+t98+t106+t111+t115+t119+t122);
t138 = sparse(t135-1.0);
t139 = sparse(t136-1.0);
t140 = sparse(t137-1.0);
t145 = sparse(-t144);
t146 = sparse(t142+t143);
t147 = sparse(t141+t145);
dh_dxn = sparse([1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3],[1,1,1,2,2,2,3,3,3,10,10,10,11,11,11,12,12,12],[alpha_u.*t32.*(t2.*t140+t3.*t129.*t134.*4.0)+Xin.*alpha_u.*t33.*t147,alpha_v.*t32.*(t3.*t139-t2.*t129.*t132.*4.0)+Yin.*alpha_v.*t33.*t147,-alpha_u.*b.*t33.*t147,-alpha_u.*t32.*(t3.*t140-t2.*t129.*t134.*4.0)+Xin.*alpha_u.*t33.*t146,alpha_v.*t32.*(t2.*t139+t3.*t129.*t132.*4.0)+Yin.*alpha_v.*t33.*t146,-alpha_u.*b.*t33.*t146,-Xin.*alpha_u.*t33.*t138-alpha_u.*t32.*t129.*(t10+t27+t35+t59+t62).*4.0,-Yin.*alpha_v.*t33.*t138+alpha_v.*t32.*t129.*(t8-t28+t34+t58+t60).*4.0,alpha_u.*b.*t33.*t138,alpha_u.*t32.*t130.*(t5+t6+t46+t51+t55-t57+t76.*4.0-t79.*3.0+t82+t83-t88+t89-t6.*t13.*3.0+t6.*t15+t4.*t18.*4.0+t5.*t18+sn3.*t14.*t30+sn1.*t20.*t30).*-8.0-Xin.*alpha_u.*t33.*t130.*t150.*4.0,alpha_v.*t32.*t130.*(t31+t50+t52-t86+X0i.*sn2.*2.0+X0i.*t17.*2.0-Y0i.*sn1.*4.0+Y0i.*t14.*4.0+Z0i.*t13.*6.0-sn3.*t4.*6.0+sn3.*t41.*2.0+sn3.*t76.*2.0+t4.*t20.*2.0+t16.*t31+X0i.*sn2.*t18.*2.0+Y0i.*sn1.*t18.*4.0-sn1.*sn2.*t4.*6.0-sn1.*sn2.*t5.*4.0+sn3.*t6.*t15.*2.0).*-4.0-Yin.*alpha_v.*t33.*t130.*t150.*4.0,alpha_u.*b.*t33.*t130.*t150.*4.0,alpha_u.*t32.*t130.*(t31+t49+t52+t86+X0i.*sn2.*4.0-X0i.*t17.*4.0-Y0i.*sn1.*2.0-Y0i.*t14.*2.0+Z0i.*t15.*6.0-sn3.*t5.*6.0+sn3.*t46.*2.0+sn2.*t56+sn3.*t79.*2.0+t5.*t20.*2.0+t19.*t31-X0i.*sn2.*t18.*4.0-Y0i.*sn1.*t18.*2.0+sn1.*sn2.*t4.*4.0+sn3.*t6.*t13.*2.0).*4.0+Xin.*alpha_u.*t33.*t130.*t149.*4.0,alpha_v.*t32.*t130.*(t4+t6+t41+t51-t53+t57-t76.*3.0+t77+t78+t79.*4.0+t87-t89+t6.*t13-t6.*t15.*3.0+t4.*t18+t5.*t18.*4.0+sn2.*t14.*t31+sn1.*t17.*t31).*-8.0+Yin.*alpha_v.*t33.*t130.*t149.*4.0,alpha_u.*b.*t33.*t130.*t149.*-4.0,alpha_u.*t32.*t130.*(t30+t45+t47-t85-X0i.*sn3.*4.0+X0i.*t20.*4.0+Y0i.*t18.*6.0+Z0i.*sn1.*2.0+Z0i.*t14.*2.0-sn2.*t6.*6.0+sn2.*t51.*2.0+sn2.*t79.*2.0+t6.*t17.*2.0+t21.*t30+X0i.*sn3.*t15.*4.0+Z0i.*sn1.*t15.*2.0-sn1.*sn3.*t4.*4.0-sn1.*sn3.*t6.*6.0+sn2.*t6.*t13.*2.0).*-4.0+Xin.*alpha_u.*t33.*t130.*t148.*8.0,alpha_v.*t32.*t130.*(t29+t42+t43+t84+X0i.*t18.*6.0+Y0i.*sn3.*4.0-Y0i.*t20.*4.0-Z0i.*sn2.*2.0-Z0i.*t17.*2.0-sn1.*t6.*6.0+sn1.*t51.*2.0+sn1.*t76.*2.0+t6.*t14.*2.0+t21.*t29-Y0i.*sn3.*t13.*4.0-Z0i.*sn2.*t13.*2.0+sn2.*sn3.*t5.*4.0+sn2.*sn3.*t6.*6.0+sn1.*t6.*t15.*2.0).*4.0+Yin.*alpha_v.*t33.*t130.*t148.*8.0,alpha_u.*b.*t33.*t130.*t148.*-8.0],3,14);
