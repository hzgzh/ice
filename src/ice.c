#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*鍐呯噧鏈哄伐浣滆繃绋嬫暟鍊艰绠楀強浼樺寲*/
/*  S 娲诲琛岀▼ 
    D 娲诲鐩村緞
    r 杩炴潌鏇插緞姣?
    */
const double Pi = 3.1416, e = 8.2, D = 0.093, S = 0.102, L = 0.152, r = 0.335, n = 3600, R = 287.08, Hu = 42500000, gb = 4 * 0.00001, L0 = 14.4, mp = 2, o = 1.0;
const double ca_begin = 208;
const double ca_lead_iv = 12;
const double ca_lag_iv = 28;
const double ca_lead_ev = 28;
const double ca_lag_ev = 12;
const double ca_lead_burn = 40;
const double ca_during_burn = 60;
const double tr = 320;

const double V(double ca) /*浣撶Н鍑芥暟 v*/

{
    /* e 鍘嬬缉姣? */
    double v, x, y;
    double lamda=L/r;
    x = S / (e - 1.0);
    y = S / 2 * ((1 + 1 / lamda) - cos(Pi / 180 * ca) - 1 / lamda * sqrt(1 - lamda * lamda * pow(sin(Pi / 180 * ca), 2)));

    v = (Pi * D * D / 4.0) * (x + y);

    return (v);
}

double AG(double v, double m, double T, double de) /*鎹㈢儹绯绘暟*/
{
    double ag;

    ag = 0.205 * (1 + 0.09) * pow((m * R * T * 1e-6 / v), 0.7) * pow((S * n / 30.0), 0.7) / pow(T, 0.2) / pow(de, 0.3);

    return (ag);
}

double AG1(double ca, double T, double m, double pa, int stage, int burnshape) /*鎹㈢儹绯绘暟*/
{
    double ag;
    double k = 1.458 - 1.628e-4 * T + 4.139e-8 * T * T;
    double pm = pow(V(ca_begin) / V(ca), k) * pa;
    double vh = Pi * D * D / 4 * S;
    double cm = n / 30 * S;
    double c1 = 0, c2 = 0;
    double x = 101350 * V(ca_begin) / 300;
    double p = m * R * T / V(ca);
    if (stage==1)
        c1 = 6.18; //杩涙帓姘旈樁娈?
    else
        c1 = 2.28; //鍘嬬缉鑶ㄨ儉闃舵
    if (burnshape==1)
        c2 = 3.24e-3; //鐩村柗寮忕噧鐑у
    else
        c2 = 6.22e-3; //鍒嗗壊寮忕噧鐑у 棰勭噧寮忕噧鐑у
    double c = c1 * cm + c2 * vh / x * (p - pm);
    ag = 820 * pow(p * 1e-6, 0.8) * pow(T, -0.53) * pow(D, -0.2) * pow(c, 0.8);

    return (ag);
}

double Qw_ca(double ag, double v, double T) /*姘旂几瀵瑰鐣岀殑鐬椂浼犵儹dQw/dca*/
{
    double qw_ca, Tw1, Tw2, Tw3, Pe = 18.4 / Pi;

    Tw1 = 120.0 + 30.0 * Pe + 273.0;
    Tw2 = 100.0 + 7.0 * Pe + 273.0;
    Tw3 = 100.0 + 4.0 * Pe + 273.0;
    qw_ca = (ag / (6* n)) * ((Pi * D * D / 4.0) * (T - Tw1) + (Pi * D * D / 4.0) * (T - Tw2) + (4.0 / D) * v * (T - Tw3));

    return (qw_ca);
}

double V_ca(double ca, double v) /******姘旂几瀹圭Н鐨勭灛鏃跺彉鍖栫巼 dv/d?******/
{
    double v_ca, x, y;
    double lamda=L/r;
    /*r 鏇茶酱鍗婂緞涓庤繛鏉嗛暱搴︿箣姣?*/
    x = sin(ca * Pi / 180.0);
    y = r * sin(ca * Pi / 180.0) * cos(ca * Pi / 180.0) / sqrt(1.0 - lamda * lamda * sin(ca * Pi / 180.0) * sin(ca * Pi / 180.0));
    v_ca = (Pi * Pi * D * D * S / 8.0 / 180) * (x + y);

    return (v_ca);
}
/*灏ゆ柉钂傚畾瀹规瘮鐑叕寮?*/
double U_T(double T, double a) /*du/dT*/
{
    double u_T, x, y;

    x = -3.0 * (0.0957 + 0.0485 / pow(a, 0.75)) * (T - 273.0) * (T - 273) * 0.000001;
    y = 2.0 * (7.768 + 3.36 / pow(a, 0.8)) * (T - 273.0) * 0.0001;
    u_T = 144.55 * (x + y + (489.6 + 46.4 / pow(a, 0.93)) * 0.01);

    return (u_T);
}
/*灏ゆ柉钂傚唴鑳藉叕寮忕悊鎯虫皵浣?*/
double U(double a, double T) /*宸ヨ川鐨勫唴鑳? u*/
{
    double x, y, z, u;

    x = -(0.0975 + 0.0485 / pow(a, 0.75)) * (T - 273.0) * (T - 273.0) * (T - 273.0) * 0.000001;
    y = (7.768 + 3.36 / pow(a, 0.8)) * (T - 273.0) * (T - 273.0) * 0.0001;
    z = (489.6 + 46.4 / pow(a, 0.93)) * (T - 273.0) * 0.01;
    u = 144.55 * (x + y + z + 1356.8);

    return (u);
}

double QB_ca(double ca, double ca_start_burn, double ca_during_burn) /*鐕冩枡鐕冪儳鏃剁殑鐬椂鏀剧儹鐜? dQB/dca*/
{
    /*mp 褰㈢姸鎸囨暟*/
    double qb_ca, y;

    y = (ca - ca_start_burn) / ca_during_burn; //?
    qb_ca = Hu * gb * 6.908 * pow(y, mp) * exp(-6.905 * pow(y, (1.0 + mp))) * ((1.0 + mp) / ca_during_burn);
    if (y < 0 && y > 1)
        qb_ca = 0;
    return (qb_ca);
}

double U_a(double a, double T) /*du/da a杩囬噺绌烘皵绯绘暟*/
{
    double x, y, z, u_a;

    x = 0.75 * (0.0485 / pow(a, 1.75)) * (T - 273.0) * (T - 273.0) * (T - 273.0) * 0.000001;
    y = -0.8 * (3.36 / pow(a, 1.8)) * (T - 273.0) * (T - 273.0) * 0.0001;
    z = -0.93 * (46.4 / pow(a, 1.93)) * (T - 273.0) * 0.01;
    u_a = 144.55 * (x + y + z);

    return (u_a);
}
/*
杩囬噺绌烘皵绯绘暟瀵硅浆瑙掑井鍒?
*/
double A_ca(double ca, double mL, double mB) /*da/dca */
{
    double a_ca;

    a_ca = -mL * QB_ca(ca, tr, ca_during_burn) / (L0 * mB * mB * Hu);

    return (a_ca);
}

/*鎷夋牸鏈楁棩鎻掑�?*/
double lgr(x, y, n, t) int n;
double t, x[], y[];
{
    int i, j, k, m;
    double z, s;
    z = 0.0;
    if (n < 1)
        return (z);
    if (n == 1)
    {
        z = y[0];
        return (z);
    }
    if (n == 2)
    {
        z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
        return (z);
    }
    i = 0;
    while ((x[i] < t) && (i < n))
        i = i + 1;
    k = i - 4;
    if (k < 0)
        k = 0;
    m = i + 3;
    if (m > n - 1)
        m = n - 1;
    for (i = k; i <= m; i++)
    {
        s = 1.0;
        for (j = k; j <= m; j++)
            if (j != i)
                s = s * (t - x[j]) / (x[i] - x[j]);
        z = z + s * y[i];
    }
    return (z);
}
double cd(double ca, double Lv, double Dv)
{
    /*瑙掑害涓庤繘鎺掓皵闃�琛岀▼%鍏崇郴*/
    double xca[15] = {0, 16.5, 32.2, 47.9, 63.6, 78.57, 94.29, 110, 125.7, 142.2, 157.1, 172.86, 188.6, 204.3, 220};
    double vlf[15] = {0, 5.94, 24.4, 48.8, 70.6, 87.46, 97.03, 100, 97.4, 87.1, 70.96, 48.85, 24.1, 6.27, 0};
    /*闈炴爣骞哄寲鐨勯榾闂ㄨ绋嬩笌闃�闂ㄩ�氬緞鍏崇郴Lv/Dv*/
    double vld[15] = {0, 0.0213, 0.0387, 0.0609, 0.0812, 0.0995, 0.12, 0.14, 0.159, 0.181, 0.202, 0.221, 0.242, 0.262, 0.280};
    double cdf[15] = {0, 0.0437, 0.09, 0.141, 0.190, 0.244, 0.282, 0.351, 0.413, 0.445, 0.470, 0.483, 0.492, 0.493, 0.494};
    double v, c;
    v = lgr(xca, vlf, 2, ca) / 100;
    c = lgr(vld, cdf, 2, v * Lv / Dv);
    if (ca < 0 || ca > 220)
        c = 0;
    return (c);
}
double Ke(double T)
{
    return 1.458 - 1.628e-4 * T + 4.139e-8 * T * T;
}
double mv_ca(double ca, double p, double pe, double T, double Lv, double Dv)
{
    double k = Ke(T);
    double Av = Pi * Dv * Dv / 4;
    double c = cd(ca, Lv, Dv);
    double f = 0;

    if (pe / p < (pow(2 / (k + 1), k / (k - 1))))
        f = c * Av * p / (6* n) / sqrt(R * T) / sqrt(k) / (1 - pow(2 / (k + 1), (k + 1) / 2 / (k - 1)));
    else
        f = c * Av * p / (6*n) / sqrt(R * T) * pow(pe / p, 1 / k) * sqrt((2 * k) / (k - 1) * (1 - pow(pe / p, (k - 1) / k)));
    if (ca < 0 || ca > 220)
        f = 0;
    return f;
}
void rkt1(double t, double y[], int n, double h, int k, double z[], void (*f)(double, double[], int, double *))
{
    int i, j, l;
    double a[4], tt, *b, *d;
    b = malloc(n * sizeof(double));
    d = malloc(n * sizeof(double));
    a[0] = h / 2.0;
    a[1] = a[0];
    a[2] = h;
    a[3] = h;
    for (i = 0; i <= n - 1; i++)
        z[i * k] = y[i];
    for (l = 1; l <= k - 1; l++)
    {
        (*f)(t, y, n, d);
        for (i = 0; i <= n - 1; i++)
            b[i] = y[i];
        for (j = 0; j <= 2; j++)
        {
            for (i = 0; i <= n - 1; i++)
            {
                y[i] = z[i * k + l - 1] + a[j] * d[i];
                b[i] = b[i] + a[j + 1] * d[i] / 3.0;
            }
            tt = t + a[j];
            (*f)(tt, y, n, d);
        }
        for (i = 0; i <= n - 1; i++)
            y[i] = b[i] + h * d[i] / 6.0;
        for (i = 0; i <= n - 1; i++)
            z[i * k + l] = y[i];
        t = t + h;
    }
    free(b);
    free(d);
    return;
}

void fun(double ca, double y[], int n, double *d)
{

    double mL, mB, W, Qh, Qb, Qex, v, ag, u_T, u_a, a, u, h, de;
    double Div, Dev, Liv, Lev;
    double mev_ca, miv_ca, T_ca, qw_ca, qb_ca, v_ca, a_ca, u_ca, mL_ca, mB_ca;
    double T, P;

    Div = 0.43 * D;
    Dev = 0.35 * D;
    Liv = 0.25 * Div;
    Lev = 0.25 * Dev;
    mL = y[0];
    mB = y[1];
    T = y[2];
    W = y[3];
    Qh = y[4];
    Qb = y[5];
    Qex = y[6];
    v = V(ca);
    h = 4 * v / (Pi * D * D);
    de = 2 * D * h / (D + 2 * h);
    v_ca = V_ca(ca, v);
    P = (mL + mB) * R * T / v;

    if ((ca >= ca_begin) && (ca < tr)) /*鍘嬬缉鏈?*/
    {
        u = U(1000, T);
        u_T = U_T(T, 1000);
        u_ca = 0;
        mB_ca = 0;
        mL_ca = 0;
        miv_ca = 0;
        mev_ca = 0;
        a_ca = 0;
        qb_ca = 0;
        ag=AG1(ca,T,mB+mL,101350,0,1);
    }

    if ((ca >= tr) && (ca < tr + ca_during_burn)) /*鐕冪儳鏈?*/
    {
        if (ca == tr)
        {
            a = 1000.0, mB = 1e-8;
            mB_ca = 0;
        }
        else
        {
            mB_ca = qb_ca / Hu;
            a = (mL) / (L0 * (mB + 1e-8));
        }
        u_T = U_T(T, a);
        u = U(a, T);
        qb_ca = QB_ca(ca, tr, ca_during_burn);
        u_a = U_a(a, T);
        a_ca = A_ca(ca, mL, (mB + 1e-8));
        u_ca = u_a * a_ca;
        mL_ca = 0;
        miv_ca = 0;
        mev_ca = 0;
        ag=AG1(ca,T,mB+mL,101350,0,1);
    }

    if ((ca >= tr + ca_during_burn) && (ca < 540 - ca_lead_ev)) /*鑶ㄨ儉鏈?*/
    {
        u_T = U_T(T, a);
        u = U(a, T);
        qb_ca = 0;
        u_ca = 0;
        mB_ca = 0;
        mL_ca = 0;
        a_ca = 0;
        miv_ca = 0;
        mev_ca = 0;
        ag=AG1(ca,T,mB+mL,101350,0,1);
    }

    /*杩涙皵鍜屾帓姘旀祦閲?*/
    double Pe = 101350;
    double an = 720 - ca_lead_iv;
    double Te = 300;
    double hi, ho, ei, eo;
    if ((ca >= 540 - ca_lead_ev) && (ca < 900 + ca_lag_iv)) /*鎹㈡皵鏈?*/
    {
        ho = U(a, T) + R * T;
        hi = U(1000, Te) + R * Te;
        if (P > Pe)
        {
            mev_ca = -mv_ca(ca - 540 + ca_lead_ev, P, Pe, T, Lev, Dev);
            miv_ca = -mv_ca(ca - an, P, Pe, T, Liv, Div);
            eo = ho * mev_ca;
            ei = ho * miv_ca;
            mB_ca = mB / (mB + mL) * (mev_ca + miv_ca);
            a_ca = 1 / (L0 * mB) * (mev_ca + miv_ca);
            mL_ca = mL / (mB + mL) * (mev_ca + miv_ca);
        }
        else
        {
            mev_ca = mv_ca(ca - 540 + ca_lead_ev, Pe, P, Te, Lev, Dev);
            miv_ca = mv_ca(ca - an, Pe, P, Te, Liv, Div);
            eo = hi * mev_ca;
            ei = hi * miv_ca;
            mB_ca = 0;
            a_ca = 1 / (L0 * mB) * (mev_ca + miv_ca);
            mL_ca = mev_ca + miv_ca;
        }

        a = mL / L0 / mB;
        u_ca = U_a(a, T) * a_ca;
        u_T = U_T(T, a);
        u = U(a, T);
        qb_ca = 0;
        ag=AG1(ca,T,mB+mL,101350,1,1);
    }
    // ag = AG(v, (mL + mB), T, de);
    qw_ca = Qw_ca(ag, v, T);
    T_ca = (-u * (mL_ca + mB_ca) + qb_ca - qw_ca - (mL + mB) * R * T / V(ca) * V_ca(ca, V(ca)) + ei + eo - (mL + mB) * u_ca) / (mL + mB) / u_T;

    d[0] = mL_ca;
    d[1] = mB_ca;
    d[2] = T_ca;
    d[3] = -(mL + mB) * R * T / V(ca) * V_ca(ca, V(ca));
    d[4] = qw_ca;
    d[5] = qb_ca;
    d[6] = eo;
}
int main()
{
    FILE *fp;

    int i, j;
    double t,max=0;
    double y[7], n, ca, h, z[7 * 7200], P = 101350;
    fp = fopen("file.txt", "w+");
    n = 4;
    ca = 208;
    h = 1;
    y[2] = 300;
    y[0] = P * V(ca) / R / y[2];
    y[1] = 0;
    y[2] = 300;
    y[3] = 0;
    y[4] = 0;
    y[5] = 0;
    y[6] = 0;
    rkt1(ca, y, 7, h, 720 / h, z, fun);
    for (i = 0; i < (int)(720 / h); i = i + 1)
    {
        t = 208.0 + i * h;
        fprintf(fp, "ca=(%f) \n", t);

        double ml = z[i];
        double mb = z[(int)(720 / h * 1) + i];
        double T = z[(int)(720 / h) * 2 + i];
        double P= (ml + mb) * R * T / V(t) / 1e6;
        if(P>max) max=P;
        fprintf(fp, "mL=%e ", z[i]);
        fprintf(fp, "mB=%e ", z[(int)(720 / h) + i]);
        fprintf(fp, "T=%e ", z[(int)(720 / h) * 2 + i]);
        fprintf(fp, "P=%e \n", P);
        fprintf(fp, "W=%e ", z[(int)(720 / h) * 3 + i]);
        fprintf(fp, "Qh=%e ", z[(int)(720 / h) * 4 + i]);
        fprintf(fp, "Qb=%e ", z[(int)(720 / h) * 5 + i]);
        fprintf(fp, "Qex=%e \n", z[(int)(720 / h) * 6 + i]);
    }
    double w=z[(int)((720/h)*3+719/h)];
    printf("单杠每循环的指示功 W %5.2f J\n",w);
    double pmi=-w*4/Pi/D/D/S/1e6;
    printf("平均指示压力Pmi %f MPa\n",pmi);
    double pm=0.0062+0.0016*max+0.00003*n;
    double eta=(pmi-pm)/pmi;
    printf("机械效率 eta %5.2f %\n",eta*100);
    printf("平均有效压力 Pme %5.2f MPa\n",pmi*eta);
    double bi=3.6e5*gb/Hu/w;
    printf("指示油耗率 bi %5.2f MPa\n",bi);
    printf("有效油耗率 be %5.2f %\n",bi/eta);
    double pe=Pi/4*D*D*S*1000*n*4*pmi*1e6*eta/30/4;
    printf("整机有效功率 Pe %5.2f kW\n",pe);
    double etai=3.6e6/Hu/bi;
    printf("指示效率 etai %5.2f %\n",etai*100);
    double etae=eta*etai;
    printf("有效效率 etae %5.2f %\n",etae*100);

    fclose(fp);
}


