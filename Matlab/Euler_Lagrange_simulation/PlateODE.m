function dy = PlateODE(t,y,v0,dv0,M,K,Ma,Ia,Q,cycles,frot)
 w =([[(5*pi^2*cos(20*pi*t))/3], [0], [(192000*pi^4*t^3)/(8000*pi^3*t^3 + 10) - (1152000000*pi^7*t^6)/(8000*pi^3*t^3 + 10)^2]]);
 dw=([[-(100*pi^3*sin(20*pi*t))/3], [0], [(576000*pi^4*t^2)/(8000*pi^3*t^3 + 10) - (11520000000*pi^7*t^5)/(8000*pi^3*t^3 + 10)^2 + (55296000000000*pi^10*t^8)/(8000*pi^3*t^3 + 10)^3]]);
 Q0 = w(2); 
 P0 = w(1); 
 Ca = funcCa(v0,w); 
 dy(1:6,1) = y(7:end); 
 dy(7:12,1) = mldivide(M,(-Q*dy(1:6) - K*y(1:6)+Ia*Ca*transpose(w)+(Q0^2+P0^2)*M*y(1:6)-Ma*[transpose(dv0);transpose(dw)]   )); 
 t;