from numba import njit

@njit
def para(i):
    atab=[0.319,0.3191,0.3326,0.3325,0.3557,0.356] #nanomet
    # cac gia tri hoping
    e1tab=[0.683,0.717,0.684,0.728,0.588,0.697]
    e2tab=[1.707,1.916,1.546,1.655,1.303,1.38]
    t0tab=[-0.146,-0.152,-0.146,-0.146,-0.226,-0.109]
    t1tab=[-0.114,-0.097,-0.13,-0.124,-0.234,-0.164]
    t2tab=[0.506,0.59,0.432,0.507,0.036,0.368]
    t11tab=[0.085,0.047,0.144,0.117,0.4,0.204]
    t12tab=[0.162,0.178,0.117,0.127,0.098,0.093]
    t22tab=[0.073,0.016,0.075,0.015,0.017,0.038]
    r0tab=[0.06,0.069,0.039,0.036,0.003,-0.015]
    r1tab=[-0.236,-0.261,-0.209,-0.234,-0.025,-0.209]
    r2tab=[0.067,0.107,0.069,0.107,-0.169,0.107]
    r11tab=[0.016,-0.003,0.052,0.044,0.082,0.115]
    r12tab=[0.087,0.109,0.06,0.075,0.051,0.009]
    u0tab=[-0.038,-0.054,-0.042,-0.061,0.057,-0.066]
    u1tab=[0.046,0.045,0.036,0.032,0.103,0.011]
    u2tab=[0.001,0.002,0.008,0.007,0.187,-0.013]
    u11tab=[0.266,0.325,0.272,0.329,-0.045,0.312]
    u12tab=[-0.176,-0.206,-0.172,-0.202,-0.141,-0.177]
    u22tab=[-0.15,-0.163,-0.15,-0.164,0.087,-0.132]
    lamtab=[0.073,0.211,0.091,0.228,0.107,0.237]
    a=atab[i]
    e1=e1tab[i]
    e2=e2tab[i]
    t0=t0tab[i]
    t1=t1tab[i]
    t2=t2tab[i]
    t11=t11tab[i]
    t12=t12tab[i]
    t22=t22tab[i]
    r0=r0tab[i]
    r1=r1tab[i]
    r2=r2tab[i]
    r11=r11tab[i]
    r12=r12tab[i]
    u0=u0tab[i]
    u1=u1tab[i]
    u2=u2tab[i]
    u11=u11tab[i]
    u12=u12tab[i]
    u22=u22tab[i]
    lam=lamtab[i]
    return a,e1,e2,t0,t1,t2,t11,t12,t22,r0,r1,r2,r11,r12,u0,u1,u2,u11,u12,u22,lam