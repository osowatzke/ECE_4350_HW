%% problem 1

% waveguide dimensions (in)
a = 0.9;
b = 0.4;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% number of modes desired
num_modes = 10;

% create row vector of all m and n values
m = zeros(1,(num_modes+1)^2);
n = zeros(1,(num_modes+1)^2);

% loop initializes row vectors
count = 1;
for i = 0:num_modes
    for j = 0:num_modes
        n(count) = i;
        m(count) = j;
        count = count + 1;
    end
end

% permeability
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% determine vector of cutoff frequencies
fc = 1./(2.*pi.*sqrt(u.*e)).*((m.*pi./a).^2+(n.*pi./b).^2).^(1/2);

% sort frequencies and indices
[fc,I] = sort(fc);
m = m(I);
n = n(I);

% loop through indices to find first 10 modes
modes_found = 0;
index = 2;
modes = strings(1,num_modes);
modes_fc = zeros(1,num_modes);
while (modes_found < 10)
    modes_found = modes_found + 1;
    modes(modes_found) = strcat("TE",int2str(m(index)),...
            int2str(n(index))); 
    modes_fc(modes_found) = fc(index)*1e-9;
    if m(index) ~= 0 && n(index) ~= 0
        modes_found = modes_found + 1;
        modes(modes_found) = strcat("TM",int2str(m(index)),...
            int2str(n(index)));
        modes_fc(modes_found) = fc(index)*1e-9;
    end
    index = index+1;
end

% output table of modes
fprintf("(1)\n\n");
T = table(modes',modes_fc','VariableNames',{'mode','fc (GHz)'});
disp(T);

%% problem 2

% waveguide dimensions in inches
a = 5.1;
b = 2.55;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% determine cutoff wavelength
Lco = 2*a;

% output cutoff wavelength
fprintf("2a) Cutoff Wavelength: %.4f m\n", Lco);

% permeability 
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% determine cutoff frequency
fco = 1/(2*a*sqrt(u*e));

% output cutoff frequency
fprintf("2a) Cutoff Frequency: %.4f GHz\n",fco*1e-9);

% frequency
f = 10e9;

% determine wavenumber
k = 2*pi*f*sqrt(u*e);

% determine propagation constant
B = k*sqrt(1-(fco/f)^2);

% output propation constant
fprintf("2b) Progation Constant: %.4f rad/m\n",B);

% determine wavelength
L = (1/sqrt(u*e))/f;

% determine characteristic impedance of dielectric
n = sqrt(u/e);

% Transverse-wave impedance
ZTE = n*(1-(L/(2*a))^2)^(-1/2);

% output transverse-wave impedance
fprintf("2b) Transverse-Wave Impedance: %.4f Ohms\n", ZTE);

% determine TEM propagation constant
B = k;

% determine TEM transverse wave impedance
ZTEM = n;

% output propogation constant
fprintf("2c) TEM Progation Constant: %.4f rad/m\n",B);

% output transverse wave impedance
fprintf("2c) TEM Transverse-Wave Impedance: %.4f Ohms\n",ZTEM);

%% problem 3

% waveguide dimensions in inches
a = 5.1;
b = 2.55;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% permeability 
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% determine cutoff frequency
fco = 1/(2*pi*sqrt(u*e))*sqrt((pi/a)^2+(pi/b)^2);

% determine cutoff wavelength
Lco = (1/sqrt(u*e))/fco;

% output cutoff wavelength
fprintf("3a) Cutoff Wavelength: %.4f m\n", Lco);

% output cutoff frequency
fprintf("3a) Cutoff Frequency: %.4f GHz\n",fco*1e-9);

% frequency
f = 10e9;

% determine wavenumber
k = 2*pi*f*sqrt(u*e);

% determine propagation constant
B = k*sqrt(1-(fco/f)^2);

% determine characteristic impedance of dielectric
n = sqrt(u/e);

% Transverse-wave impedance
ZTM = n*sqrt(1-(fco/f)^2);

% output propation constant
fprintf("3b) Progation Constant: %.4f rad/m\n",B);

% output transverse-wave impedance
fprintf("3b) Transverse-Wave Impedance: %.4f Ohms\n", ZTM);

% determine TEM propagation constant
B = k;

% determine TEM transverse wave impedance
ZTEM = n;

% output propogation constant
fprintf("3c) TEM Progation Constant: %.4f rad/m\n",B);

% output transverse wave impedance
fprintf("3c) TEM Transverse-Wave Impedance: %.4f Ohms\n",ZTEM);

%% problem 4

% waveguide dimensions in inches
a = 5.1;
b = 2.55;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% frequency
f = 10e9;
w = 2*pi*f;

% permeability
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% determine cutoff frequency
fco = 1/(2*a*sqrt(u*e));

% determine wavenumber
k = 2*pi*f*sqrt(u*e);

% determine propagation constant
B = k*sqrt(1-(fco/f)^2);

% determine amplitude of Hz wave
kx = pi/a;
A = kx*10/(1i*B);

% output results
fprintf("(4)\tEx=Ez=Hy=0\n");
fprintf("\tHz = -j%.4fcos(%.4fx) A/m\n",abs(imag(A)),kx);
fprintf("\tEy = %.4fsin(%.4fx) V/m\n",real(-1i*w*u*A/kx),kx);
fprintf("\tHx = 10sin(%.4fx) A/m\n",kx);

%% problem 5

% waveguide dimensions in inches
a = 0.75;
b = 0.375;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% permeability 
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% characteristic impedance of the dielectric
n = sqrt(u/e);

% cutoff frequency
fc = 1/(2*a*sqrt(u*e));

% frequency
f = 12e9;

% determine transverse wave impedance
ZTE = n/sqrt(1-(fc/f)^2);

% maximum electric field
E0 = 2e6;

% maximum power transfer
PT = E0^2*b*a/(4*ZTE);

% output maximum power transfer
fprintf("(5) Maximum Power Transfer: %.4f kW\n", PT*1e-3);

%% problem 6

% waveguide dimensions in inches
a = 0.75;
b = 0.375;

% convert waveguide dimensions to meters
a = a*0.0254;
b = b*0.0254;

% note both materials are nonmagnetic and share same permeability
u = 4*pi*1e-7;

% dielectric permittivity
e = 2.1*8.854e-12;

% conductor conducitity
s = 5.8e7;

% frequency 
f = 12e9;

% conductivity of conductor 
Rs = sqrt(pi*f*u/s);

% characteristic impedance of dielectric
n = sqrt(u/e);

% wavelength
L = 1/(f*sqrt(u*e));

% conductor loss
ac = Rs/(b*n*sqrt(1-(L/(2*a))^2))*(1+2*b/a*(L/(2*a))^2);

% output conductor loss
fprintf("(6)\tConductor Loss: %.4f Np/m = %.4f dB/m\n", ac,...
    20*log10(exp(ac)));

% loss tangent
loss_tan = 5e-4;

% wave number
k = 2*pi*f*sqrt(u*e);

% dielectric loss
ad = k*loss_tan/(2*sqrt(1-(L/(2*a))^2));

% output dielectric loss
fprintf("\tDielectric Loss: %.4f Np/m = %.4f dB/m\n", ad,...
    20*log10(exp(ad)));

%% problem 7

% plate separation in mm
a = 15;

% plate seperation in meters
a = a*1e-3;

% symbolic variable to represent mode and frequecy
syms n;
syms f;

% permeability 
u = 4*pi*1e-7;

% permittivity
e = 8.854e-12;

% display cutoff frequency expression
fc = n/(2*a*sqrt(u*e));
fprintf("7a)\tfc = ");
disp(vpa(fc,4));

% determine modal characteristic impedance expression
ZTE = sqrt(u/e)/sqrt(1-(fc/f)^2);
fprintf("\tZTE = ");
disp(vpa(ZTE,4));

% for TM modes frequency expression is the same
fprintf("7b)\tfc = ");
disp(vpa(fc,4));

% determine modal characteristic impedance expression
ZTM = sqrt(u/e)*sqrt(1-(fc/f)^2);
fprintf("\tZTM = ");
disp(vpa(ZTE,4));

% determine number of modes that would propagate if f = 33Ghz
n = floor(double(solve(fc == 33e9)));

% output results
fprintf("7c) %d TE modes would propagate.\n",n);

%% problem 8

% plate spacing in cm
a = 1.5;

% plate spacing in m
a = a*1e-2;

% loss tangent is zero for air
% => air has no losses
fprintf("(8)\tWith Air Dielectric\n");
fprintf("\tDielectric Loss: 0 dB/m\n");

%frequency
f = 12e9;

% permeability (same for all cases and materials)
u = 4*pi*1e-7;

% permittivity of air
e = 8.854e-12;

% conducitivity
s = 5.8e7;

% surface resistance
Rs = sqrt(pi*f*u/s);

% characteristic impedance of air dielectric
n = sqrt(u/e);

% cutoff frequency with air dielectric
fc = 1/(2*a*sqrt(u*e));

% conductor loss;
ac = 2*Rs/(n*a*sqrt(1-(fc/f)^2));

% output conductor loss
fprintf("\tConductor Loss: %.4f dB/m\n\n", 20*log10(exp(ac)));

% permittivity of glass dielectric
e = 4*8.854e-12;

% characteristic impedance of glass dielectric
n = sqrt(u/e);

% cutoff frequency with glass dielectric
fc = 1/(2*a*sqrt(u*e));

% conductor loss
ac = 2*Rs/(n*a*sqrt(1-(fc/f)^2));

% loss tangenent
loss_tan = 2e-3;

% dielectric loss
ad = 2*pi*f*sqrt(u*e)*loss_tan/(2*sqrt(1-(fc/f)^2));

% output losses
fprintf("\tWith Glass Dielectric\n");
fprintf("\tDielectric Loss: %.4f dB/m\n", 20*log10(exp(ad)));
fprintf("\tConductor Loss: %.4f dB/m\n", 20*log10(exp(ac)));













