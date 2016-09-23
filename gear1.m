% This is a try to calculate the gears on octave 
# Setting the variables of on stage of the gear train,
# every stage must be calculate separately:
# VARIAVEIS ::
filename="engrenagem.txt"; 
id = fopen(filename, "w");
#t = strftime ("%r (%Z) %A %e %B %Y", localtime (time ()))
val = save_header_format_string ();
tm = strftime (val, localtime (time ())); 
fprintf(id, "%s \n", tm)
#=============================================================
#PASSO NORMAL 
pn = 0.5; % in 
# ANGULO DE HELICE
PSI = pi/6; % rad => 30 graus
#ANGULO DE PRESSAO 
phi = pi/9; %rad => 20 graus
# NUMERO DE DENTES DO PINHAO
Np = 28 %dentes
# NUMERO DE DENTES DA ENGRENAGEM
Ng = 56	%dentes
# TORQUE DO PINHAO[lbf*in]
Tp = 14338.215
# VELOCIDADE ANGULAR DO PINHAO [rpm]
omega_p = 3.3333;
#RAZAO DE TRASMISSAO
mv = 3; % razao de transmissao
# Maxima temperatura de operacao em celsius
T = 200;
# LARGURA DA FACE
F = 2.8;  %in
# NUMERO DE CICLOS [VIDA INFINITA]
N = 1E7
#***************************************************************** 
# DEFININDO FATORES DE PARA TENSAO DE FLEXAO 
# FATOR DE APLICACAO (Ka)
Ka = 1; % suave 
# FATOR DE TAMANHO (Ks)
Ks = 1; % AGMA
# FATOR DE ESPESSURA DA BORDA (Kb)
Kb = 1; % Engrenagem solida
# FATOR DE CICLO DE CARGA (Ki)
Ki = 1; % Engrenagem n solta
# FATOR DE FABRICACAO
Qv = 8;
#******************************************************************
#CORECAO PARA A RESISTENCIA A FADIGA DE FLEXAO
# FATOR DE CONFIABILIDADE(Kr)
Kr = 1; % Confiabilidade de 99% 
# AGMA GRAU DO MATERIAL UTILIZADO
Gr = 3;
# DUREZA BRINNEL DO MATERIAL
# SE utilizado ligas de aco deve se definir Alloy =1 e a HB sera a 
#dureza da parte nao tratada superficialmente 
# Veriicar a norma AGMA 2001-D04 para info adicional  
HB = 325
Alloy = 1
#******************************************************************
# FATORES PARA A TENSAO DE SUPERFICIE
#FATOR DE ACABAMENTO DE SUPERFICIE (Cf)
Cf = 1; % metodos de fabricacao convencionais
#****************************************************************** 
#CORRECAO PARA A RESISTENCIA A FADIGA DE SUPERFICIE
# FATOR DE RAZAO DE DUREZA(Ch)
Ch = 1; % Mesmo material do pinhao e engrenagem
#*******************************************************************
#PROPRIEDADES DOS MATERIAIS 
# Poison do pinhao e da engrenagem
poison_p = 0.28;
poison_g = 0.28;
# Modulo de Young
Ep = 30E6;
Eg = 30E6;
#===================================================================
# DO NOT MODIFY BELOW HERE!!!
#===================================================================
# FATOR GEOMETRICO DE FLEXAO(J)
# FATOR GEOMETRICO DE FLEXAO PARA phi = 20 e psi = 30
# profundidade completa
Jf = [14, 14, 0.39, 0.39; ...
    14, 17, 0.39, 0.41; ...
    14, 21, 0.40, 0.43; ...
    14, 26, 0.41, 0.44; ...
    14, 35, 0.41, 0.46; ...
    14, 55, 0.42, 0.49; ...
    14, 135, 0.43, 0.51; ...
    17, 17, 0.41, 0.41; ...
    17, 21, 0.42, 0.43; ...
    17, 26, 0.43, 0.45; ...
    17, 35, 0.43, 0.47; ...
    17, 55, 0.44, 0.49; ...
    17, 135, 0.45, 0.52; ...
    21, 21, 0.44, 0.44; ...
    21, 26, 0.45, 0.46; ...
    21, 35, 0.45, 0.48; ...
    21, 55, 0.46, 0.50; ...
    21, 135, 0.47, 0.53; ...
    26, 26, 0.46, 0.46; ...
    26, 35, 0.47, 0.48; ...
    26, 55, 0.48, 0.50; ...
    26, 135, 0.49, 0.53; ...
    35, 35, 0.49, 0.49; ...
    35, 55, 0.50, 0.51; ...
    35, 135, 0.51, 0.54; ...
    55, 55, 0.52, 0.52; ...
    55, 135, 0.53, 0.55; ...
    135, 135, 0.56, 0.56];
#
x1 = Jf(:,1);
y1 = Jf(:,2);
jp = Jf(:,3);
jg = Jf(:,4);
J_p = griddata(x1,y1,jp,Np,Ng, "linear")
J_g = griddata(x1, y1, jg, Np, Ng, "linear")
#
# FATOR DE CARGA (Km)
Lm = [2, 6, 9, 20];
Kn = [1.6, 1.7, 1.8, 2.0]; 
Km = interp1(Lm, Kn, F, "linear") 
#==================================================================
# Calculating.
# PASSO TRANSVERSAL;
pt = pn/cos(PSI)  % in
#PASSO AXIAL
px = pn/sin(PSI) % in
# PASSO DIAMETRAL;
pd = pi/pt % dentes/ in 
# PASSO DIAMETRAL NO PLANO NORMAL;
pnd = pd /cos(PSI)  % dentes/in
#DIAMETRO DO PINHAO
dp = Np/pd  %in
dg = Ng/pd  %in
rp = dp/2;  %in
rg = dg/2;  %in
# ANGULO DE PRESSAO NORMAL E TANGENCIAL
phi_n = atan(tan(phi)*cos(PSI));
phi_t = phi;
# ANGULO DE HELICE DA BASE (PSI_b)
PSI_b = acos(cos(PSI)*(cos(phi_n)/cos(phi)));
# 
Wt = Tp/rp  %lbf
# 
Wr = Wt*tan(phi); %lbf
#
W = Wt/(cos(PSI)*cos(phi_n)); %lbf
#
# RAIO DA ELIPSE
re_p = rp/(cos(PSI))^2; %in
re_g = rg/(cos(PSI))^2; %in
#
#NUMERO VIRTUAL DE DENTES
Ne_p = 2*pi*re_p/pn;    
Ne_g = 2*pi*re_g/pn;
#
#RAZAO DE CONTATO AXIAL (mf)
mf = F/px;  
# Checando a razao de contato da engrenagem
    if (mf< 1)
    printf("Engrenagem LACR \n")
    elseif(mf>(1.15))
    printf("Engrenagem Convencional \n")
    endif
#
# CALCULO DAS TENSOES; 
#
# Vt deve estar em ft/min entao a formula esta corrigida
Vt = (rp/12)*(omega_p)*(2*pi); % ft/min
#FATOR DINAMICO(Kv)
if (Qv>5)
    B = ((12 - Qv)^(2/3))/4;
    A = 50 + 56*(1-B);    
    Kv = (A/(A+sqrt(Vt)))^B;
    Vtmax = (A + (Qv - 3))^2;
    elseif(Qv_p<=5)
    Kv = (50/(50+sqrt(Vt)));
endif 
#
#
sigma_bp = (Wt*pd*Ka*Km*Ks*Kb*Ki)/F*J_p*Kv  %psi
sigma_bg = (Wt*pd*Ka*Km*Ks*Kb*Ki)/F*J_g*Kv  %psi
#
# Calculo do adendo, dedendo, folga e Distancia entre centros
ap = 1/pd; %in 
ag = ap; %in
bp = 1.25/pd; %in 
bg = bp; %in
ht = ap + bp; %in 
c = bp - ap; %in
C = rp + rg %in 
# COMPRIMENTO DE ACAO
Z = sqrt((rp + ap)^2 - (rp*cos(phi))^2) + ...
    sqrt((rg + ag)^2 - (rg*cos(phi))^2) - C*sin(phi); %in  
#
mp = (pd*Z)/pi*cos(phi);
mF = (F*pd*tan(PSI))/pi;
#
# COMPRIMENTO MINIMO DAS LINHAS DE CONTATO
nr = abs(mp - floor (mp));
na = abs(mF - floor (mF));
#
x = 1 - nr;
#
    if (na > x)
    Lmin =  (mp*F - (1 - na)*(1 - nr)* px)/cos(PSI_b);
    elseif(na <= x)
    Lmin = (mp*F - na*nr*px)/cos(PSI_b);
    endif 
#
#RAIOS DE CURVATURA
pho_p= sqrt((0.5*((rp + ap) + (C - rg - ag)))^2 - (rp*cos(phi))^2);
pho_g = C*sin(phi) - pho_p;
#
#RAZAO DE DIVISAO DE CARGA
mN = F/Lmin 
# 
#FATOR DE GEOMETRIA DA SUPERFICIE 
# engrenagens externas 
Ic = cos(phi)/(((1/pho_p) + (1/pho_g))*dp*mN); 
#
#TENSAO DE SUPERFICIE 
#FATORES DE TENSAO DE SUPERFICIE 
Ca = Ka;
Cm = Km;
Cv = Kv;
Cs = Ks;
#
#COEFICIENTE ELASTICO
Cp = sqrt(1/(pi*(((1-poison_p^2)/Ep) + ((1-poison_g^2)/Eg))))
#
sigma_cp = Cp * sqrt((Wt*Ca*Cm*Cs*Cf)/(F*Ic*dp*Cv))
sigma_cg = Cp * sqrt((Wt*Ca*Cm*Cs*Cf)/(F*Ic*dg*Cv)) 
#
#FATOR DE VIDA PARA RESISTENCIA A FLEXAO
Kl = 1.6831*N^(-0.0323);
# FATOR DE TEMPERATURA (KT)
T_F = (9/5)*(T + 273.15) - 459.67; 
KT = (460 + T_F)/620;
#RESISTENCIA A FADIGA DE FLEXAO NAO CORRIGIDA
    if (Gr == (1) & Alloy == 0)
        sBF = -274 + 167*HB - 0.152*HB^2;  %psi
    elseif(Gr == (1) &  Alloy==1)
        sBF = 105.2*HB + 9280 
    elseif(Gr == 2 & Alloy ==1)
        sBF = 105.2*HB + 22280
    elseif((Gr == 3) & (Alloy ==1))
        sBF = 105.2*HB + 29280
    endif 
# RESISTENCIA A FADIGA DE FLEXAO CORRIGIDA 
Sbf = (Kl*sBF)/KT*Kr  %psi
#
# CORRECAO PARA A RESISTENCIA A FADIGA DE SUPERFICIE
#
CT = KT;, Cr = Kr;
#FATOR DE RESISTENCIA A FADIGA DE SUPERFICIE 
Cl = 2.466*N^(-0.056)
# 
    if(Gr==(1) & Alloy == 0)
        sCF = 26000 + 327*HB %psi
    elseif(Gr==(2) & Alloy==0)
        sCF = 27000 + 364*HB %psi
    elseif(Gr==1 & Alloy==1)
        sCF = 176000
    elseif(Gr==2 & Alloy==1)
        sCF = 196000
    elseif(Gr==3 & Alloy==1)
        sCF = 216000
    endif 
Scf = Cl*Ch*sCF/CT*Cr   %psi
#
#COEFICIENTE DE SEGURANCA DE FLEXAO (Nb)
# pinhao 
Nbp = Sbf/sigma_bp
# engrenagem 
Nbg = Sbf/sigma_bg
# 
#COEFICIENTE DE SEGURANCA DE SUPERFICIE(Nc)
# pinhao
Ncp = (Scf/sigma_cp)^2
# engrenagem 
Ncg = (Scf/sigma_cg)^2
#=============================================================
fprintf(id, "\n")
fprintf(id, "   Programming by:                 Paulo Conci\n")
fprintf(id, "\n\n")
fprintf(id, "   DADOS INICIAIS: \n")
fprintf(id, "   Passo normal = %f [in]\n", pn)
fprintf(id, "   Angulo de helice = %g [rads] \n", PSI)
fprintf(id, "   Angulo de pressao = %g [rads] \n", phi)
fprintf(id, "   Numero de dentes do pinhao = %g \n", Np)
fprintf(id, "   Numero de dentes da engrenagem = %g \n", Ng)
fprintf(id, "   Torque aplicado ao pinhao = %g [lbf*in]\n", Tp)
fprintf(id, "   Velocidade angular do pinhao = %g [rpm]\n", omega_p)
fprintf(id, "   Maxima temperatura de operacao = %g [Celsius]\n", T)
fprintf(id, "   Largura da face = %g [in] \n", F)
fprintf(id, "   Numero de ciclos estimados de vida = %g \n", N)
fprintf(id, "   Material utilizado = Nitroalloy 135M Gr 2\n")
fprintf(id, "   Dureza Core do material utilizado = %g [HB]\n", HB)
fprintf(id, "\n\n")
fprintf(id, "   DADOS CALCULADOS :\n")
fprintf(id, "\n")
fprintf(id, "   Passo transversal = %g [in]\n", pt)
fprintf(id, "   Passo axial = %g [in] \n", px)
fprintf(id, "   Passo diametral = %g [dentes/in] \n", pd)
fprintf(id, "   Diametro de referencia do pinhao = %g [in]\n", dp)
fprintf(id, "   Diametro de referencia da coroa = %g [in]\n", dg)
fprintf(id, "   Angulo de pressao normal = %g [rads]\n", phi_n)
fprintf(id, "   Angulo de helice da base = %g [rads]\n", PSI_b)
fprintf(id, "   Adendo = %g [in]\n", ap)
fprintf(id, "   Dedendo = %g [in]\n", bp)
fprintf(id, "   Altura total = %g [in]\n", ht)
fprintf(id, "   Folga = %g [in]\n", c)
fprintf(id, "   Distancia entre centros = %g [in]\n", C)
fprintf(id, "   Comprimento de acao = %g [in]\n", Z)
fprintf(id, "   Comprimento minimo das linhas de contato = %g [in]\n", Lmin)
fprintf(id, "   Raio de curvatura do pinhao = %g [in]\n", pho_p)
fprintf(id, "   Raio de curvatura da coroa = %g [in]\n", pho_g)
fprintf(id, "   Razao de contato axial = %g [dentes]\n", mf)
fprintf(id, "   Razao da divisao de carga = %g \n", mN)
fprintf(id, "\n")
fprintf(id, "   Carga tangencial aplicada no pinhao = %g [lbf]\n", Wt)
fprintf(id, "   Carga radial aplicada no pinhao = %g [lbf]\n", Wr)
fprintf(id, "   Carga total aplicada = %g [lbf]\n", W)
fprintf(id, "   Velocidade tangencial do pinhao = %g [ft/min]\n", Vt)
fprintf(id, "\n")
fprintf(id, "   Tensao de flexao aplicada ao pinhao = %g [psi]\n", sigma_bp)
fprintf(id, "   Tensao de flexao aplicada a coroa = %g [psi]\n", sigma_bg)
fprintf(id, "   Tensao de superficie no pinhao = %g [psi]\n", sigma_cp)
fprintf(id, "   Tensao de superficie na coroa = %g [psi]\n", sigma_cg)
fprintf(id, "\n")
fprintf(id, "   Resistencia a fadiga de flexao nao corrigida = %g [psi]\n", sBF)
fprintf(id, "   Resistencia a fadiga de flexao corrigida = %g [psi]\n", Sbf)
fprintf(id, "\n")
fprintf(id, "   Resistencia a fadiga de superficie nao corrigida = %g [psi]\n", sCF)
fprintf(id, "   Resistencia a fadiga de superficie corrigida = %g [psi]\n", Scf)
fprintf(id, "\n")
fprintf(id, "   Coeficientes de seguranca\n")
fprintf(id, "   Flexao do pinhao = %g \n", Nbp)
fprintf(id, "   Flexao da coroa = %g \n", Nbg)
fprintf(id, "   Superficie do pinhao = %g \n", Ncp)
fprintf(id, "   Superficie da coroa = %g \n",Ncg)
fclose(id);
