clear 
small = 1e-6
[A1,B1,C1,D1] = tf2ss(50*[1 202 401 200],[1 25 250 1500 5000]);

tau = 0.0001;
Ti = 0.11;
Tt=0.11;
K = 0.07;

A = [A1 B1 -(K/tau)*B1; zeros(1,5) -K/(Ti*tau) ; C1 0 -1/(tau)];

Bp = [-B1;-1/Tt;0];
Bw = [K*B1;K/Ti;0];

Cq = [zeros(1,4) 1 -K/tau]
Cz = [C1  zeros(1,2)]

Dqw = K;
Dqp = 0;
Dzw = -1;
Dzp = 0;

a = 2
d = 1 

P = sdpvar(6,6);
a = sdpvar(1)
b = sdpvar(1)
gamma = sdpvar(1);
IQC = [a 0; 0 -b];
IQCp = [gamma^2 0; 0 -1];

Bigmat = [eye(6) zeros(6,2); A Bp Bw; Cq Dqp Dqw; zeros(1,6) 1 0;  Cz Dzp Dzw; zeros(1,7) 1];

BigP = [zeros(6,6) P zeros(6,4); P zeros(6,10);
    zeros(2,12) IQC zeros(2,2); zeros(2,14) IQCp]

Constraints = [P>=small*eye(6), Bigmat'*BigP*Bigmat<=-small*eye(8), a-b>=small]

Cost = gamma
options = sdpsettings('solver','sdpt3','verbose',1);
sol = optimize(Constraints,Cost,options);
