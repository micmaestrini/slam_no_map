function partials_w = STM_omega_w(k1,k2,s1,s2,s3,wcx,wcy,wcz,wx,wy,wz)
%STM_OMEGA_W
%    PARTIALS_W = STM_OMEGA_W(K1,K2,S1,S2,S3,WCX,WCY,WCZ,WX,WY,WZ)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    04-Nov-2019 18:33:09

t2 = sparse(s1.^2);
t3 = sparse(s2.^2);
t4 = sparse(s3.^2);
t5 = sparse(t2+t3+t4+1.0);
t6 = sparse(1.0./t5.^2);
t7 = sparse(t2.*8.0);
t8 = sparse(t3.*8.0);
t9 = sparse(t7+t8);
t10 = sparse(t6.*t9);
t11 = sparse(t10-1.0);
t12 = sparse(t11.*wcz);
t13 = sparse(exp(k2));
t14 = sparse(s2.*t2);
t15 = sparse(s2.*t4);
t37 = sparse(s1.*s3.*2.0);
t60 = sparse(s2.*t3);
t16 = sparse(-s2+t14+t15-t37+t60);
t17 = sparse(t6.*t16.*wcx.*4.0);
t18 = sparse(s2.*s3.*2.0);
t19 = sparse(s1.*t3);
t20 = sparse(s1.*t4);
t21 = sparse(-s1+t18+t19+t20+s1.*t2);
t22 = sparse(exp(-k1));
t23 = sparse(exp(k1));
t24 = sparse(t4.*8.0);
t25 = sparse(t7+t24);
t26 = sparse(t6.*t25);
t27 = sparse(t26-1.0);
t28 = sparse(t27.*wcy);
t29 = sparse(s1.*s2.*2.0);
t30 = sparse(s3.*t2);
t31 = sparse(s3.*t3);
t61 = sparse(s3.*t4);
t32 = sparse(-s3+t29+t30+t31+t61);
t33 = sparse(-s1-t18+t19+t20+s1.*t2);
t34 = sparse(t6.*t33.*wcz.*4.0);
t36 = sparse(t6.*t32.*wcx.*4.0);
t35 = sparse(t28+t34-t36);
t59 = sparse(t6.*t21.*wcy.*4.0);
t38 = sparse(t12+t17-t59);
t39 = sparse(t23.*t38);
t40 = sparse(t13.*wz);
t41 = sparse(s2.*wcx.*4.0);
t42 = sparse(s1.*t2.*wcy.*4.0);
t43 = sparse(t2.^2);
t44 = sparse(t43.*wcz);
t45 = sparse(t4.*wcz.*2.0);
t46 = sparse(t3.^2);
t47 = sparse(t46.*wcz);
t48 = sparse(t4.^2);
t49 = sparse(t48.*wcz);
t50 = sparse(s1.*t3.*wcy.*4.0);
t51 = sparse(s1.*t4.*wcy.*4.0);
t52 = sparse(t2.*t3.*wcz.*2.0);
t53 = sparse(t2.*t4.*wcz.*2.0);
t54 = sparse(t3.*t4.*wcz.*2.0);
t55 = sparse(s1.*s3.*wcx.*8.0);
t56 = sparse(s2.*s3.*wcy.*8.0);
t57 = sparse(t41+t42+t44+t45+t47+t49+t50+t51+t52+t53+t54+t55+t56+wcz-s1.*wcy.*4.0-t2.*wcz.*6.0-t3.*wcz.*6.0-s2.*t2.*wcx.*4.0-s2.*t3.*wcx.*4.0-s2.*t4.*wcx.*4.0);
t58 = sparse(t6.*t13.*t57);
t62 = sparse(t8+t24);
t63 = sparse(t6.*t62);
t64 = sparse(t63-1.0);
t65 = sparse(t64.*wcx);
t66 = sparse(-s3-t29+t30+t31+t61);
t67 = sparse(t6.*t66.*wcy.*4.0);
t68 = sparse(t23.*t35);
t69 = sparse(t13.*t35);
t70 = sparse(exp(-k2));
t71 = sparse(-s2+t14+t15+t37+s2.*t3);
t72 = sparse(t6.*t71.*wcz.*4.0);
t73 = sparse(s3.*wcy.*4.0);
t74 = sparse(t2.*wcx.*2.0);
t75 = sparse(t43.*wcx);
t76 = sparse(t46.*wcx);
t77 = sparse(t48.*wcx);
t78 = sparse(s2.*t3.*wcz.*4.0);
t79 = sparse(s2.*t2.*wcz.*4.0);
t80 = sparse(s2.*t4.*wcz.*4.0);
t81 = sparse(t2.*t3.*wcx.*2.0);
t82 = sparse(t2.*t4.*wcx.*2.0);
t83 = sparse(t3.*t4.*wcx.*2.0);
t84 = sparse(s1.*s2.*wcy.*8.0);
t85 = sparse(s1.*s3.*wcz.*8.0);
t86 = sparse(t73+t74+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+wcx-s2.*wcz.*4.0-t3.*wcx.*6.0-t4.*wcx.*6.0-s3.*t2.*wcy.*4.0-s3.*t3.*wcy.*4.0-s3.*t4.*wcy.*4.0);
t87 = sparse(-s3-t29+t30+t31+s3.*t4);
t88 = sparse(t6.*t87.*wcy.*4.0);
partials_w = sparse([2,3,1,3,1,2],[1,1,2,2,3,3],[t12+t17+t39+t40+t58-t59-t23.*wz,-t70.*(-t28-t34+t36+t68+t69+wy-t23.*wy),-t22.*(t12+t17+t39+t40+t58-wz-t6.*t21.*wcy.*4.0),t70.*(t65-t72+t88-wx+t23.*wx+t13.*(t65-t72+t88)+t6.*t23.*t86),t22.*(-t28-t34+t36+t68+t69+wy-t13.*wy),-t65-t67+t72+t13.*wx-t23.*wx-t13.*(t65+t67-t6.*wcz.*(-s2+t14+t15+t37+t60).*4.0)-t6.*t23.*t86],3,3);
