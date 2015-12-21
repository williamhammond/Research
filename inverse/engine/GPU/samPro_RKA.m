function [u,v,tag]=samPro_RKA(u,v,par,sti,dt_mat,C)
MAX_RK = 10;
v0 = v;
v_init = v;

u0 = u;
u_init = u;

[u0,v0]=fp(u0,v0,sti,C,par);

a=dt_mat/6.0;
u= u+ u0*a;
v= v+ v0*a;

if max(max(abs(u)))>MAX_RK ||max(max(abs(v)))>MAX_RK
    fprintf('Step size too large\n');
    u = u_init;
    v = v_init;
    tag = 1;
    return
end

a = dt_mat/2.0;
u0 = u0*a;
u0 = u0 + u_init;
v0 = v0*a;
v0 = v0 + v_init;

[u0,v0]=fp(u0,v0,sti,C,par);
a=dt_mat/3.0;
u= u+ u0*a;
v= v+ v0*a;

if max(max(abs(u)))>MAX_RK ||max(max(abs(v)))>MAX_RK
    fprintf('Step size too large\n');
    u = u_init;
    v = v_init;
    tag = 1;
    return
end

a = dt_mat/2.0;
u0 = u0*a;
u0 = u0 + u_init;
v0 = v0*a;
v0 = v0 + v_init;

[u0,v0]=fp(u0,v0,sti,C,par);
a=dt_mat/3.0;
u= u+ u0*a;
v= v+ v0*a;

if max(max(abs(u)))>MAX_RK ||max(max(abs(v)))>MAX_RK
    fprintf('Step size too large\n');
    u = u_init;
    v = v_init;
    tag = 1;
    return
end

u0 = u0*dt_mat;
u0 = u0 + u_init;
v0 = v0*dt_mat;
v0 = v0 + v_init;

[u0,v0]=fp(u0,v0,sti,C,par);

a=dt_mat/6.0;
u= u+ u0*a;
v= v+ v0*a;

if max(max(abs(u)))>MAX_RK ||max(max(abs(v)))>MAX_RK
    fprintf('Step size too large\n');
    u = u_init;
    v = v_init;
    tag = 1;
    return
end
tag = 0;
end