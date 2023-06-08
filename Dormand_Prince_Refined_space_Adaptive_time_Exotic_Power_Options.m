% 1. DPC-Loc2.........Gamma-3.......

clc
clear all
close all

format long
m = 0.5;    % power term of the options. 
            % Could be any arbitrary rational number

rd = 0.05;   % risk free interest for each regime
sigm = 0.2;  % constant volatility for each regime
sigma = m*sigm; % new volatility depending on m

K = 100;      % asset price at expiration time
T = 3;        % expiration time
Tol = 10^-2;  % Tolerance

h = 0.05; dt = 0.0025; x_max = 3; 
nx = x_max/h; nt = T/dt;   % rectangular domain of x and t (mesh)
                     
sfm(1) = K; % allocating first column to be the initial
                        % optimal exercise boundary  

% begin.........adding few grid points locally to the left boundary                        
haa = h/4;
hx(1) = h/4; hx(2) = h/4; hx(3) = h/4; hx(4) = h/4; 
hx(5) = h/2; hx(6) = h/2;  
for i = 7:nx+4
    hx(i)  = h;
end

xx(1) = 0;
for i = 2:nx+5
    xx(i) = xx(i-1)+hx(i-1);
end
% end.........adding few grid points locally to the left boundary 

% begin.........intializing variables 
for i = 1:nx+5
    U_old(i) = 0;  
    W_old(i) = 0;
    Y_old(i) = 0;
end
% end.........intializing variables

% begin.........formulating the non-uniform compact differencing
for i = 1:nx+3
    ab1(i) = (hx(i+1)/(hx(i)+hx(i+1)))*((((hx(i))^2)+hx(i)*hx(i+1)-...
        ((hx(i+1))^2))/(((hx(i))^2)+3*hx(i)*hx(i+1)+((hx(i+1))^2)));
    ab3(i) = (hx(i)/(hx(i)+hx(i+1)))*((((hx(i+1))^2)+hx(i)*hx(i+1)-...
        ((hx(i))^2))/(((hx(i))^2)+3*hx(i)*hx(i+1)+((hx(i+1))^2)));
    
    bc1(i) = (hx(i+1)/(hx(i)+hx(i+1)))*(12/(((hx(i))^2)+3*hx(i)*hx(i+1)...
        +((hx(i+1))^2)));
    bc2(i) = -(12/(((hx(i))^2)+3*hx(i)*hx(i+1)+((hx(i+1))^2)));
    bc3(i) = (hx(i)/(hx(i)+hx(i+1)))*(12/(((hx(i))^2)+3*hx(i)*hx(i+1)...
        +((hx(i+1))^2)));
end
% end.........formulating the non-uniform compact differencing

% begin.........Discrete matrix system
for i=1:nx+3
    if i == 1
        B(i,i) = -24/haa^2; B(i,i+1) = 12/haa^2;
        N(i,i) = 6; 
        A(i,i) = -120/haa^2; A(i,i+1) = -30/(2*haa^2);
        M(i,i) = 10;
    elseif i == nx+3
        B(i,i-1) = bc1(i); B(i,i) = bc2(i); 
        N(i,i-1) = ab1(i); N(i,i) = 1; 

        A(i,i-1) = bc1(i); A(i,i) = bc2(i); 
        M(i,i-1) = ab1(i); M(i,i) = 1; 
    else
        B(i,i-1) = bc1(i); B(i,i) = bc2(i); B(i,i+1) = bc3(i);
        N(i,i-1) = ab1(i); N(i,i) = 1; N(i,i+1) = ab3(i);

        A(i,i-1) = bc1(i); A(i,i) = bc2(i); A(i,i+1) = bc3(i);
        M(i,i-1) = ab1(i); M(i,i) = 1; M(i,i+1) = ab3(i);
    end
end
% end.........Discrete matrix system

% begin.........further initialization of parameters
kk = 1; t = 0; rr = 0; dta = dt; rf = 0;
mb1 = 0.5*sigma^2; mc1 = rd; nu = m*(rd-rf-0.5*sigm^2);

a20 = 2/sigma^2; 
a21 = rd-m*rd+(m*sigm^2)/2; 
a22 = a20*a21;

a1 = (85/18)+(2*rd*haa^2)/sigma^2; 
a2 = (85/18)+(11/3)*haa+a22*haa^2;  
     
b0 = 85/18; b1 = 6; 
b2 = 3/2; b3 = 2/9; 
% end.........further initialization of parameters

tic
while t<T
    
    % begin.................................................
    nn = 0; 
    if t+dt>T
      dt = T-t;  % ensure that the t+dt<=T
    end
    % end.................................................
    
    while(1)
        
        % begin.................................................
        U1 = (2/haa^2)*(U_old(1)-2*U_old(2)+U_old(3))-(1/(2*haa))*...
            (W_old(3)-W_old(1));
        U2 = (2/haa^2)*(U_old(2)-2*U_old(3)+U_old(4))-(1/(2*haa))*...
            (W_old(4)-W_old(2));
        U3 = (2/haa^2)*(U_old(3)-2*U_old(4)+U_old(5))-(1/(2*haa))*...
            (W_old(5)-W_old(3));

        Uh1 = (0.5*sigma^2)*(b1*U1-b2*U2+b3*U3);
        Uh2 = mc1*(b1*U_old(2)-b2*U_old(3)+b3*U_old(4));
        Wh1 = nu*(b1*W_old(2)-b2*W_old(3)+b3*W_old(4));
    
        h0 = a2+(1/sfm(kk))*(b1*W_old(2)-b2*W_old(3)+b3*W_old(4));
    
        Ra1 =  -(Uh1-Uh2+Wh1)/h0; 
        s1 = Ra1/sfm(kk)+nu;
     
        Hu1(1) = (12/haa^2)*(K-sfm(kk))-(3/haa)*(W_old(3)+sfm(kk)); 
        Hu1(nx+3) = 0; Hu1 = Hu1';

        Hw1(1) = ((15*sfm(kk))/haa^2)+(150/(2*haa^3))*(U_old(3)-...
            (K-sfm(kk))); 
        Hw1(nx+3) = 0; Hw1 = Hw1'; 


        y0 = (2*rd*K)/sigma^2-a22*sfm(kk);
        Hy1(1) = (150/(2*haa^3))*(W_old(3)+(sfm(kk)))-((15*y0)...
            /haa^2);
        Hy1(nx+3) = 0; Hy1 = Hy1';

        Ru1 = (mb1)*(N\(B*U_old(2:nx+4)'+Hu1))+s1*W_old(2:nx+4)'-mc1*...
            U_old(2:nx+4)'; 
        Ua1 = U_old(2:nx+4)'+(dt/5)*Ru1; 
    
        sfa = -(6*Ua1(1)-1.5*Ua1(2)+(2/9)*Ua1(3)-a1*K)/a2;    
        Rb1 = -(6*Ru1(1)-1.5*Ru1(2)+(2/9)*Ru1(3))/a2; 
            
        Rw1 = (mb1)*(M\(A*W_old(2:nx+4)'+Hw1))+s1*(N\(B*U_old...
            (2:nx+4)'+Hu1))-mc1*W_old(2:nx+4)';
        Wa1 = W_old(2:nx+4)'+(dt/5)*Rw1;  

        Ry1 = (mb1)*(M\(A*Y_old(2:nx+4)'+Hy1))+s1*(M\(A*W_old(2:nx+4)'...
            +Hw1))-mc1*Y_old(2:nx+4)';
        Ya1 = Y_old(2:nx+4)'+(dt/5)*Ry1; 
        % end.................................................
    
        % begin.................................................
        s1 = Rb1/sfa+nu;
  
        Iu1(1) = (12/haa^2)*(K-sfa)-(3/haa)*(Wa1(2)+sfa); 
        Iu1(nx+3) = 0; Iu1 = Iu1';

        Iw1(1) = ((15*sfa)/haa^2)+(150/(2*haa^3))*(Ua1(2)-(K-sfa));  
        Iw1(nx+3) = 0; Iw1 = Iw1'; 
 
        y0 = (2*rd*K)/sigma^2-a22*sfa;
        Iy1(1) = (150/(2*haa^3))*(Wa1(2)+(sfa))-((15*y0)/haa^2);
        Iy1(nx+3) = 0; Iy1 = Iy1';
    
        Su1 = (mb1)*(N\(B*Ua1+Iu1))+s1*Wa1-mc1*Ua1;
        Ub1 = U_old(2:nx+4)'+dt*((3/40)*Ru1+(9/40)*Su1); 
      
        sfb = -(6*Ub1(1)-1.5*Ub1(2)+(2/9)*Ub1(3)-a1*K)/a2; 
        Rc1 = -(6*Su1(1)-1.5*Su1(2)+(2/9)*Su1(3))/a2; 
    
        Sw1 = (mb1)*(M\(A*Wa1+Iw1))+s1*(N\(B*Ua1+Iu1))-mc1*Wa1;
        Wb1 = W_old(2:nx+4)'+dt*((3/40)*Rw1+(9/40)*Sw1);

        Sy1 = (mb1)*(M\(A*Ya1+Iy1))+s1*(M\(A*Wa1+Iw1))-mc1*Ya1;
        Yb1 = Y_old(2:nx+4)'+dt*((3/40)*Ry1+(9/40)*Sy1);
        % end.................................................
    
        % begin.................................................
        s1 = Rc1/sfb+nu;
    
        Ju1(1) = (12/haa^2)*(K-sfb)-(3/haa)*(Wb1(2)+sfb); 
        Ju1(nx+3) = 0; Ju1 = Ju1';

        Jw1(1) = ((15*sfb)/haa^2)+(150/(2*haa^3))*(Ub1(2)-(K-sfb));  
        Jw1(nx+3) = 0; Jw1 = Jw1';

        y0 = (2*rd*K)/sigma^2-a22*sfb;
        Jy1(1) = (150/(2*haa^3))*(Wb1(2)+(sfb))-((15*y0)/haa^2);
        Jy1(nx+3) = 0; Jy1 = Jy1';

        Tu1 = (mb1)*(N\(B*Ub1+Ju1))+s1*Wb1-mc1*Ub1;    
        Uc1 = U_old(2:nx+4)'+dt*((44/45)*Ru1-(56/15)*Su1+(32/9)*Tu1); 
    
        sfc = -(6*Uc1(1)-1.5*Uc1(2)+(2/9)*Uc1(3)-a1*K)/a2; 
        Rd1 = -(6*Tu1(1)-1.5*Tu1(2)+(2/9)*Tu1(3))/a2; 
    
        Tw1 = (mb1)*(M\(A*Wb1+Jw1))+s1*(N\(B*Ub1+Ju1))-mc1*Wb1;    
        Wc1 = W_old(2:nx+4)'+dt*((44/45)*Rw1-(56/15)*Sw1+(32/9)*Tw1); 

        Ty1 = (mb1)*(M\(A*Yb1+Jy1))+s1*(M\(A*Wb1+Jw1))-mc1*Yb1;
        Yc1 = Y_old(2:nx+4)'+dt*((44/45)*Ry1-(56/15)*Sy1+(32/9)*Ty1);
        % end.................................................
    
        % begin.................................................
        s1 = Rd1/sfc+nu;
    
        Ku1(1) = (12/haa^2)*(K-sfc)-(3/haa)*(Wc1(2)+sfc); 
        Ku1(nx+3) = 0; Ku1 = Ku1';

        Kw1(1) = ((15*sfc)/haa^2)+(150/(2*haa^3))*(Uc1(2)-(K-sfc));  
        Kw1(nx+3) = 0; Kw1 = Kw1';

        y0 = (2*rd*K)/sigma^2-a22*sfc;
        Ky1(1) = (150/(2*haa^3))*(Wc1(2)+(sfc))-((15*y0)/haa^2);
        Ky1(nx+3) = 0; Ky1 = Ky1';

    
        Uu1 = (mb1)*(N\(B*Uc1+Ku1))+s1*Wc1-mc1*Uc1;    
        Ud1 = U_old(2:nx+4)'+dt*((19372/6561)*Ru1-(25360/2187)*Su1+...
            (64448/6561)*Tu1-(212/729)*Uu1); 
    
        sfd = -(6*Ud1(1)-1.5*Ud1(2)+(2/9)*Ud1(3)-a1*K)/a2; 
        Re1 = -(6*Uu1(1)-1.5*Uu1(2)+(2/9)*Uu1(3))/a2; 
    
        Uw1 = (mb1)*(M\(A*Wc1+Kw1))+s1*(N\(B*Uc1+Ku1))-mc1*Wc1;    
        Wd1 = W_old(2:nx+4)'+dt*((19372/6561)*Rw1-(25360/2187)*Sw1+...
            (64448/6561)*Tw1-(212/729)*Uw1);

        Uy1 = (mb1)*(M\(A*Yc1+Ky1))+s1*(M\(A*Wc1+Kw1))-mc1*Yc1;
        Yd1 = Y_old(2:nx+4)'+dt*((19372/6561)*Ry1-(25360/2187)*Sy1+...
            (64448/6561)*Ty1-(212/729)*Uy1);
        % end.................................................
       
        % begin.................................................
        s1 = Re1/sfd+nu;
    
        Lu1(1) = (12/haa^2)*(K-sfd)-(3/haa)*(Wd1(2)+sfd);  
        Lu1(nx+3) = 0; Lu1 = Lu1';

        Lw1(1) = ((15*sfd)/haa^2)+(150/(2*haa^3))*(Ud1(2)-(K-sfd));  
        Lw1(nx+3) = 0; Lw1 = Lw1';

        ya = (11/(3*haa))*sfd-(85/(18*haa^2))*(K-sfd);
        yb = (1/haa^2)*(6*Ud1(1)-1.5*Ud1(2)+(2/9)*Ud1(3)); 
        y0 = (2*rd*K)/sigma^2-a22*sfd;

        Ly1(1) = (150/(2*haa^3))*(Wd1(2)+(sfd))-((15*y0)/haa^2);
        Ly1(nx+3) = 0; Ly1 = Ly1';
    
        Vu1 = (mb1)*(N\(B*Ud1+Lu1))+s1*Wd1-mc1*Ud1;    
        Ue1 = U_old(2:nx+4)'+dt*((9017/3168)*Ru1-(355/33)*Su1+(46732/...
            5247)*Tu1+(49/176)*Uu1-(5103/18656)*Vu1); 
    
        sfe = -(6*Ue1(1)-1.5*Ue1(2)+(2/9)*Ue1(3)-a1*K)/a2; 
        Rf1 = -(6*Vu1(1)-1.5*Vu1(2)+(2/9)*Vu1(3))/a2; 
    
        Vw1 = (mb1)*(M\(A*Wd1+Lw1))+s1*(N\(B*Ud1+Lu1))-mc1*Wd1;    
        We1 = W_old(2:nx+4)'+dt*((9017/3168)*Rw1-(355/33)*Sw1+(46732/...
            5247)*Tw1+(49/176)*Uw1-(5103/18656)*Vw1);

        Vy1 = (mb1)*(M\(A*Yd1+Ly1))+s1*(M\(A*Wd1+Lw1))-mc1*Yd1;
        Ye1 = Y_old(2:nx+4)'+dt*((9017/3168)*Ry1-(355/33)*Sy1+(46732/...
            5247)*Ty1+(49/176)*Uy1-(5103/18656)*Vy1);
        % end.................................................
        
        % begin.................................................
        s1 = Rf1/sfe+nu;
    
        Mu1(1) = (12/haa^2)*(K-sfe)-(3/haa)*(We1(2)+sfe);  
        Mu1(nx+3) = 0; Mu1 = Mu1';

        Mw1(1) = ((15*sfe)/haa^2)+(150/(2*haa^3))*(Ue1(2)-(K-sfe)); 
        Mw1(nx+3) = 0; Mw1 = Mw1';

        ya = (11/(3*haa))*sfe-(85/(18*haa^2))*(K-sfe);
        yb = (1/haa^2)*(6*Ue1(1)-1.5*Ue1(2)+(2/9)*Ue1(3)); 
        y0 = (2*rd*K)/sigma^2-a22*sfe;
      
        My1(1) = (150/(2*haa^3))*(We1(2)+(sfe))-((15*y0)/haa^2);
        My1(nx+3) = 0; My1 = My1';
    
        Wu1 = (mb1)*(N\(B*Ue1+Mu1))+s1*We1-mc1*Ue1;    
        Uf1 = U_old(2:nx+4)'+dt*((35/384)*Ru1+(500/1113)*Tu1+(125/192)...
            *Uu1-(2187/6784)*Vu1+(11/84)*Wu1); 
    
        sff = -(6*Uf1(1)-1.5*Uf1(2)+(2/9)*Uf1(3)-a1*K)/a2; 
        Rg1 = -(6*Wu1(1)-1.5*Wu1(2)+(2/9)*Wu1(3))/a2; 
    
        Ww1 = (mb1)*(M\(A*We1+Mw1))+s1*(N\(B*Ue1+Mu1))-mc1*We1;    
        Wf1 = W_old(2:nx+4)'+dt*((35/384)*Rw1+(500/1113)*Tw1+(125/192)...
            *Uw1-(2187/6784)*Vw1+(11/84)*Ww1);

        Wy1 = (mb1)*(M\(A*Ye1+My1))+s1*(M\(A*We1+Mw1))-mc1*Ye1;
        Yf1 = Y_old(2:nx+4)'+dt*((35/384)*Ry1+(500/1113)*Ty1+(125/192)...
            *Uy1-(2187/6784)*Vy1+(11/84)*Wy1);
        
        ab = isreal(sum(Uf1)); ac = isreal(sum(Wf1)); ad = isreal(sff);
        if ab == 1 && ac == 1 && ad == 1
            break;
        else
            dt = 0.5*dt;    % establising maximum time step
            Hu1 = Hu1';  Hw1 = Hw1'; Hy1 = Hy1';
            Iu1 = Iu1';  Iw1 = Iw1'; Iy1 = Iy1';
            Ju1 = Ju1';  Jw1 = Jw1'; Jy1 = Jy1';
            Ku1 = Ku1';  Kw1 = Kw1'; Ky1 = Ky1';
            Lu1 = Lu1';  Lw1 = Lw1'; Ly1 = Ly1'; 
            Mu1 = Mu1';  Mw1 = Mw1'; My1 = My1';
        end       
    % end.................................................
    end    
    
    % begin.................................................    
    s1 = Rg1/sff+nu;
    
    Nu1(1) = (12/haa^2)*(K-sff)-(3/haa)*(Wf1(2)+sff);  
    Nu1(nx+3) = 0; Nu1 = Nu1';

    Nw1(1) = ((15*sff)/haa^2)+(150/(2*haa^3))*(Uf1(2)-(K-sff)); 
    Nw1(nx+3) = 0; Nw1 = Nw1';
 
    y0 = (2*rd*K)/sigma^2-a22*sff;
      
    Ny1(1) = (150/(2*haa^3))*(Wf1(2)+(sff))-((15*y0)/haa^2);
    Ny1(nx+3) = 0; Ny1 = Ny1';
    
    Xu1 = (mb1)*(N\(B*Uf1+Nu1))+s1*Wf1-mc1*Uf1;    
    Un1 = Uf1;
    sfnew = -(6*Un1(1)-1.5*Un1(2)+(2/9)*Un1(3)-a1*K)/a2;
    
    Un = U_old(2:nx+4)'+dt*((5179/57600)*Ru1+(7571/16695)*Tu1+...
      (393/640)*Uu1-(92097/339200)*Vu1+(187/2100)*Wu1+(1/40)*Xu1);
    sfne = -(6*Un(1)-1.5*Un(2)+(2/9)*Un(3)-a1*K)/a2;
    
    Xw1 = (mb1)*(M\(A*Wf1+Nw1))+s1*(N\(B*Uf1+Nu1))-mc1*Wf1;    
    Wn1 = Wf1;
    
    Wn = W_old(2:nx+4)'+dt*((5179/57600)*Rw1+(7571/16695)*Tw1+...
      (393/640)*Uw1-(92097/339200)*Vw1+(187/2100)*Ww1+(1/40)*Xw1);

    Xy1 = (mb1)*(M\(A*Yf1+Ny1))+s1*(M\(A*Wf1+Nw1))-mc1*Yf1;
    Yn1 = Yf1;
    
    Yn = Y_old(2:nx+4)'+dt*((5179/57600)*Ry1+(7571/16695)*Ty1+...
      (393/640)*Uy1-(92097/339200)*Vy1+(187/2100)*Wy1+(1/40)*Xy1);
    % end.......................................................
    
    % begin: defining error threshold
    Reu(1) = max(abs(Un1-Un)); Reu(2) = max(abs(Wn1-Wn)); 
    Reu(3) = abs(sfnew-sfne); Reu(4) = max(abs(Yn1-Yn)); 
    aw = abs(max(Reu));
    % end: defining error threshold
   
   % begin: transfer of variables.
   if aw<=Tol
       U_iter1(1) = K-sfnew; U_iter1(nx+5) = 0; 
       for i = 2:nx+4
           U_iter1(i) = Un1(i-1);
       end
       
       W_iter1(1) = -sfnew; W_iter1(nx+5) = 0;  
       for i = 2:nx+4
          W_iter1(i) = Wn1(i-1);
       end

       Y_iter1(1) = (2*rd*K)/sigma^2-a22*sfnew; Y_iter1(nx+5) = 0;  
       for i = 2:nx+4
          Y_iter1(i) = Yn1(i-1);
       end
       % end: transfer of variables.
       sfm(kk+1) = sfnew; sf=sfnew  % Storing the early exercise boundary
       fff(kk) = dt; kk = kk+1;
        
        for i = 1:nx+5
            U_old(i) = U_iter1(i); 
            W_old(i) = W_iter1(i); % assigning values
            Y_old(i) = Y_iter1(i);
        end       
        Hu1 = Hu1';  Hw1 = Hw1'; Hy1 = Hy1';
        Iu1 = Iu1';  Iw1 = Iw1'; Iy1 = Iy1';
        Ju1 = Ju1';  Jw1 = Jw1'; Jy1 = Jy1';
        Ku1 = Ku1';  Kw1 = Kw1'; Ky1 = Ky1';
        Lu1 = Lu1';  Lw1 = Lw1'; Ly1 = Ly1'; 
        Mu1 = Mu1';  Mw1 = Mw1'; My1 = My1';
        Nu1 = Nu1';  Nw1 = Nw1'; Ny1 = Ny1';
        t = t+dt;  
        dt = 0.9*dt*((Tol/abs(max(Reu)))^(0.25));
        % end: transfer values
    else
        % begin: select the optimal step size
        dt = 0.9*dt*((Tol/abs(max(Reu)))^(0.2));
        Hu1 = Hu1';  Hw1 = Hw1'; Hy1 = Hy1';
        Iu1 = Iu1';  Iw1 = Iw1'; Iy1 = Iy1';
        Ju1 = Ju1';  Jw1 = Jw1'; Jy1 = Jy1';
        Ku1 = Ku1';  Kw1 = Kw1'; Ky1 = Ky1';
        Lu1 = Lu1';  Lw1 = Lw1'; Ly1 = Ly1'; 
        Mu1 = Mu1';  Mw1 = Mw1'; My1 = My1';
        Nu1 = Nu1';  Nw1 = Nw1'; Ny1 = Ny1';
        rr = rr+1;
    end    % change the column
end

% begin.......visualization procedure

ay(1) = 0;
for i = 2:length(fff)+1
ay(i) = ay(i-1)+fff(i-1);
end
                                % plot parameter
ayy(1) = 0;
for i = 2:length(fff)
ayy(i) = ay(i-1)+fff(i-1);
end

for i = 1:nx+5
    S1(i) = exp(xx(i))*sfnew; % computing asset value
end

for i = 1:nx+6
    if i==1
        Smain1(i) = 0;  
        Smain3(i) = 0;
        Vmain1(i) = K - Smain1(i); 
        Vmain3(i) = K - Smain3(i);  % option value
    else
        Smain1(i) = S1(i-1); 
        Smain3(i) = K+(i-2);       
        Vmain1(i) = U_iter1(i-1);
        Vmain3(i) = 0;
    end
end

Wmain1(1) = -1;
for i = 2:nx+6
    Wmain1(i) = W_iter1(i-1)/S1(i-1);  % delta sensitivity
end

Ymain1(1) = 0;
for i = 2:nx+6
    if i==2
        Ymain1(i) = 0;
    else
        Ymain1(i) = (Y_iter1(i-1)-W_iter1(i-1))/... % Gamma
        (S1(i-1)*S1(i-1)); 
    end
end

% start: interpolated data
Tableu(1) = interp1(Smain1,Vmain1,80); 
Tableu(2) = interp1(Smain1,Vmain1,90,'spline');  
Tableu(3) = interp1(Smain1,Vmain1,100,'spline');  
Tableu(4) = interp1(Smain1,Vmain1,110,'spline'); 
Tableu(5) = interp1(Smain1,Vmain1,120,'spline');     

Tablew(1) = interp1(Smain1,Wmain1,80); 
Tablew(2) = interp1(Smain1,Wmain1,90,'spline'); 
Tablew(3) = interp1(Smain1,Wmain1,100,'spline'); 
Tablew(4) = interp1(Smain1,Wmain1,110,'spline'); 
Tablew(5) = interp1(Smain1,Wmain1,120,'spline');

Tabley(1) = interp1(Smain1,Ymain1,80); 
Tabley(2) = interp1(Smain1,Ymain1,90,'pchip'); 
Tabley(3) = interp1(Smain1,Ymain1,100,'pchip'); 
Tabley(4) = interp1(Smain1,Ymain1,110,'pchip'); 
Tabley(5) = interp1(Smain1,Ymain1,120,'pchip');
% end: interpolated data

% start: visualization
figure(1)
plot(ay,sfm,'LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Optimal Exercise Boundary','FontSize',14);
get(o)
p=xlabel('Time','FontSize',14);
get(p)
grid on
axis([0,0.5,75,100])

figure(2)
subplot(1,2,1)
plot(Smain1,Vmain1,'b','LineWidth',2)
hold on
plot(Smain3,Vmain3,'k','LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Asset Option','FontSize',14);
get(o)
p=xlabel('Asset Price','FontSize',14);
get(p)
grid on
axis([0,200,0,100])

subplot(1,2,2)
plot(Smain1,Wmain1,'LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Delta Option','FontSize',14);
get(o)
p=xlabel('Asset Price','FontSize',14);
get(p)
grid on
axis([0 200 -1 0])

figure(3)
plot(Smain1,Ymain1,'LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Gamma Option','FontSize',14);
get(o)
p=xlabel('Asset Price','FontSize',14);
get(p)
grid on
axis([0,200,0,0.1])
% end: visualization
% end.......visualization procedure
toc