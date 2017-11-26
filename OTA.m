% Optimises tracking error to get minimum value using Lagrange's Multiplier Method

W = [8.24; 14.75; 2.14; 20.21; 19.35; 4.85; 15.43; 2.64; 3.94; 8.46];
Rp = [36.74; 20.58; 25.89; 36.24; 34.36; 29.67; 27.28; 18.45; 8.49; 10.02];
Ri = [40.56; 22.68; 22.27; 33.21; 38.74; 37.63; 26.23; 22.73; 6.49; 8.75];


I = ones(10, 1);


r = [17.02, 33.11, 30.04, 32.20, 35.37, 26.19, 69.22, 17.93, 52.47, 43.85;
  8.20, 7.57, 32.43, 15.78, 17.93, 16.00, 50.22, 7.86, 32.33, 14.03;
  31.23, 32.22, 11.11, 13.75, 29.11, 7.42, 54.23, 19.23, 27.75, 32.85;
  10.11, 25.75, 41.12, 24.56, 57.75, 23.47, 27.85, 29.04, 44.79, 52.11;
  34.56, 19.25, 56.12, 12.11, 19.97, 24.56, 32.31, 40.45, 45.25, 14.97;
  23.45, 28.93, 45.65, 12.75, 20.56, 9.95, 45.23, 15.25, 50.35, 20.75;
  47.89, 17.96, 48.22, 22.37, 35.89, 44.67, 60.34, 30.03, 40.44, 10.23;
  3.72, 5.67, 12.33, 23.45, 14.57, 11.12, 1.76, 2.34, 9.87, 2.04;
  40.02, 25.75, 27.31, 11.65, 21.56, 17.35, 2.33, 7.89, 13.76, 21.34;
  7.20, 14.53, 2.35, 18.97, 11.45, 7.89, 20.21, 4.45, 7.02, 6.22];


V = zeros(10, 10);


for i = 1:10
  for j = 1:10
    if i == j
      xbar = sum(r(i,:)) / 10;
      V(i, j) = sum((r(i,:).-xbar).^2) / 10 ;

    else
      xbar = sum(r(i,:)) / 10;
      ybar = sum(r(j,:)) / 10;
      V(i, j) = ((r(i,:).-xbar) * (r(j,:).-ybar)') / 10;
    end
  end
end;


% Covariance matrix
V(:,:)


A = Rp'*inv(V)*I;
B = Rp'*inv(V)*Rp;
C = I'*inv(V)*I;


Rp1 = Rp' * W;
M = [B, A; A, C];
lambda = zeros(2, 1);
lambda = inv(M) * [Rp1; 1];

W1 = inv(V)*(Rp*lambda(1, 1) + I*lambda(2, 1))

meanRp = sum(Rp)/10;
meanRi = sum(Ri)/10;

stdRp = sqrt(sum((Rp.-meanRp).^2) / 9);
stdRi = sqrt(sum((Ri.-meanRi).^2) / 9);

% Objective function
volatility = (W1)' * V * (W1);   % sigma square

Zxp = (Rp.-meanRp)/stdRp;
Zxi = (Ri.-meanRi)/stdRi;

Zf = Zxp .* Zxi;
roh = sum(Zf)/9;

% Tracking error in %
Tracking_err = sqrt(volatility) * sqrt(1-roh^2) / 10  % Sigma
%printf(Tracking_err);
