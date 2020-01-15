function dh_dy0 = jacobian_meas_n_meas_0(Xin,Yin,Zin,alpha_u,alpha_v,b,d0,s01,s02,s03,sn1,sn2,sn3,u,u0,v,v0)
%JACOBIAN_MEAS_N_MEAS_0
%    DH_DY0 = JACOBIAN_MEAS_N_MEAS_0(XIN,YIN,ZIN,ALPHA_U,ALPHA_V,B,D0,S01,S02,S03,SN1,SN2,SN3,U,U0,V,V0)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Jan-2020 14:24:53

t2 = sparse(s01.^2);
t3 = sparse(s01.^3);
t4 = sparse(s02.^2);
t5 = sparse(s02.^3);
t6 = sparse(s03.^2);
t7 = sparse(s03.^3);
t8 = sparse(sn1.^2);
t9 = sparse(sn1.^3);
t10 = sparse(sn2.^2);
t11 = sparse(sn2.^3);
t12 = sparse(sn3.^2);
t13 = sparse(sn3.^3);
t14 = sparse(s01.*s02.*2.0);
t15 = sparse(s01.*s03.*2.0);
t16 = sparse(s02.*s03.*2.0);
t17 = sparse(sn1.*sn2.*2.0);
t18 = sparse(sn1.*sn3.*2.0);
t19 = sparse(sn2.*sn3.*2.0);
t20 = sparse(1.0./Zin);
t22 = sparse(1.0./alpha_v);
t23 = sparse(1.0./d0);
t25 = sparse(-s01);
t26 = sparse(-s02);
t27 = sparse(-s03);
t28 = sparse(-sn1);
t29 = sparse(-sn2);
t30 = sparse(-sn3);
t31 = sparse(-u0);
t32 = sparse(-v0);
t21 = sparse(t20.^2);
t24 = sparse(t23.^2);
t33 = sparse(-t14);
t34 = sparse(-t15);
t35 = sparse(-t16);
t36 = sparse(-t17);
t37 = sparse(-t18);
t38 = sparse(-t19);
t39 = sparse(s01.*t4);
t40 = sparse(s02.*t2);
t41 = sparse(s01.*t6);
t42 = sparse(s03.*t2);
t43 = sparse(s02.*t6);
t44 = sparse(s03.*t4);
t45 = sparse(sn1.*t10);
t46 = sparse(sn2.*t8);
t47 = sparse(sn1.*t12);
t48 = sparse(sn3.*t8);
t49 = sparse(sn2.*t12);
t50 = sparse(sn3.*t10);
t51 = sparse(t2.*8.0);
t52 = sparse(t4.*8.0);
t53 = sparse(t6.*8.0);
t54 = sparse(t8.*8.0);
t55 = sparse(t10.*8.0);
t56 = sparse(t12.*8.0);
t57 = sparse(t31+u);
t58 = sparse(t32+v);
t59 = sparse(t2+t4+t6+1.0);
t60 = sparse(t8+t10+t12+1.0);
t61 = sparse(t51+t52);
t62 = sparse(t51+t53);
t63 = sparse(t52+t53);
t64 = sparse(t54+t55);
t65 = sparse(t54+t56);
t66 = sparse(t55+t56);
t67 = sparse(1.0./t59.^2);
t68 = sparse(1.0./t60.^2);
t69 = sparse(t3+t16+t25+t39+t41);
t70 = sparse(t5+t15+t26+t40+t43);
t71 = sparse(t7+t14+t27+t42+t44);
t72 = sparse(t9+t19+t28+t45+t47);
t73 = sparse(t11+t18+t29+t46+t49);
t74 = sparse(t13+t17+t30+t48+t50);
t75 = sparse(t3+t25+t35+t39+t41);
t76 = sparse(t5+t26+t34+t40+t43);
t77 = sparse(t7+t27+t33+t42+t44);
t78 = sparse(t9+t28+t38+t45+t47);
t79 = sparse(t11+t29+t37+t46+t49);
t80 = sparse(t13+t30+t36+t48+t50);
t81 = sparse(t61.*t67);
t82 = sparse(t62.*t67);
t83 = sparse(t63.*t67);
t84 = sparse(t64.*t68);
t85 = sparse(t65.*t68);
t86 = sparse(t66.*t68);
t98 = sparse(alpha_v.*t20.*t68.*t72.*4.0);
t99 = sparse(alpha_u.*t20.*t68.*t74.*4.0);
t100 = sparse(alpha_u.*t20.*t68.*t79.*4.0);
t101 = sparse(alpha_v.*t20.*t68.*t80.*4.0);
t102 = sparse(Xin.*alpha_u.*t21.*t68.*t73.*4.0);
t103 = sparse(Yin.*alpha_v.*t21.*t68.*t73.*4.0);
t104 = sparse(Xin.*alpha_u.*t21.*t68.*t78.*4.0);
t105 = sparse(Yin.*alpha_v.*t21.*t68.*t78.*4.0);
t115 = sparse(alpha_u.*b.*t21.*t67.*t68.*t71.*t78.*1.6e+1);
t116 = sparse(alpha_u.*b.*t21.*t67.*t68.*t73.*t77.*1.6e+1);
t87 = sparse(t81-1.0);
t88 = sparse(t82-1.0);
t89 = sparse(t83-1.0);
t90 = sparse(t84-1.0);
t91 = sparse(t85-1.0);
t92 = sparse(t86-1.0);
t106 = sparse(-t105);
t119 = sparse(t99+t104);
t120 = sparse(t101+t103);
t93 = sparse(alpha_u.*t20.*t92);
t94 = sparse(alpha_v.*t20.*t91);
t95 = sparse(Xin.*alpha_u.*t21.*t90);
t96 = sparse(Yin.*alpha_v.*t21.*t90);
t107 = sparse(alpha_u.*b.*t21.*t67.*t69.*t90.*4.0);
t108 = sparse(alpha_u.*b.*t21.*t68.*t73.*t89.*4.0);
t109 = sparse(alpha_u.*b.*t21.*t67.*t76.*t90.*4.0);
t110 = sparse(alpha_u.*b.*t21.*t68.*t78.*t88.*4.0);
t126 = sparse(t88.*t119);
t127 = sparse(t89.*t120);
t133 = sparse(t67.*t71.*t119.*4.0);
t134 = sparse(t67.*t77.*t120.*4.0);
t97 = sparse(-t95);
t111 = sparse(-t109);
t112 = sparse(-t110);
t113 = sparse(t93+t102);
t114 = sparse(t96+t98);
t117 = sparse(t94+t106);
t130 = sparse(t67.*t69.*(t95-t100).*-4.0);
t131 = sparse(-t127);
t132 = sparse(t67.*t76.*(t95-t100).*-4.0);
t118 = sparse(t97+t100);
t121 = sparse(t89.*t113);
t122 = sparse(t88.*t117);
t123 = sparse(t67.*t69.*t114.*4.0);
t124 = sparse(t67.*t77.*t113.*4.0);
t125 = sparse(t67.*t76.*t114.*4.0);
t129 = sparse(t67.*t71.*t117.*4.0);
t135 = sparse(t107+t112+t116);
t136 = sparse(t108+t111+t115);
t128 = sparse(-t124);
t137 = sparse(t122+t123+t134);
t138 = sparse(t121+t132+t133);
t140 = sparse(t125+t129+t131);
t139 = sparse(t126+t128+t130);
dh_dy0 = sparse([1,2,3,1,2,3,1,2,3],[1,1,1,2,2,2,3,3,3],[b.*t23.*t138,-b.*t23.*t140,-b.*t23.*t136,alpha_u.*b.*t22.*t23.*(t124-t126+t67.*t69.*(t95-t100).*4.0),alpha_u.*b.*t22.*t23.*t137,-alpha_u.*b.*t22.*t23.*t135,b.*t24.*t57.*t138-alpha_u.*b.*t24.*(t87.*(t95-t100)+t67.*t70.*t113.*4.0+t67.*t75.*t119.*4.0)+alpha_u.*b.*t22.*t24.*t58.*(t124-t126+t67.*t69.*(t95-t100).*4.0),-b.*t24.*t57.*t140-alpha_u.*b.*t24.*(t87.*t114+t67.*t70.*t120.*4.0-t67.*t75.*t117.*4.0)+alpha_u.*b.*t22.*t24.*t58.*t137,-b.*t24.*t57.*t136+alpha_u.*b.*t24.*(alpha_u.*b.*t21.*t87.*t90+alpha_u.*b.*t21.*t67.*t68.*t70.*t73.*1.6e+1+alpha_u.*b.*t21.*t67.*t68.*t75.*t78.*1.6e+1)-alpha_u.*b.*t22.*t24.*t58.*t135],3,3);
