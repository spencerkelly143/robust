clear 
small = 1e-6
[A1,B1,C1,D1] = tf2ss(50*[1 202 401 200],[1 25 250 1500 5000]);

tau = sdpvar(1);
Ti = 0.11;
Tt=0.11;
K = 0.07;

A = [A1 B1 -(K/tau)*B1; zeros(1,5) -K/(Ti*tau) ; C1 0 -1/(tau)];

B = [-B1;-1/Tt;0];

C = [zeros(1,4) 1 -K/tau]

D = 0;

P = sdpvar(6,6);
gamma =sdpvar(1)
IQC = [1 0; 0 1];

BigLmi = [eye(6) zeros(6,1); A B]'*[zeros(6,6) P; P zeros(6,6)]*[eye(6) zeros(6,1); A B] %+ [C D; zeros(1,6) 1]'*IQC*[C D; zeros(1,6) 1];
Constraints = [P>=small*eye(6), tau>=0.0001]
%, BigLmi <= -small
Cost = 0
options = sdpsettings('solver','sdpt3','verbose',1);
sol = optimize(Constraints,Cost,options);

