%% Define Problem Data
[A1,B1,C1,D1] = tf2ss(50*[1 202 401 200],[1 25 250 1500 5000]);

Ti = 0.11;
Tt = 0.11;
K = 0.07;

tol = 1e-6;

%% Circle Criterion, MQP block

tau = 0.001;
del = 0;
while del < 1
    tau = tau + tol;

    Aqp = [A1, B1, -B1*K/tau; zeros(1,4), 0, -K/(tau * Ti); C1, 0, -1/tau];
    Bqp = [-B1; -1/Tt; 0];
    Cqp = [zeros(1,4), 1, -K/tau];
    Dqp = 0;
    
    [num, dem] = ss2tf(Aqp, Bqp, Cqp, Dqp);
    [re,im,wout] = nyquist(tf(num, dem));
    del = max(re);    
end

nyquist(tf(num, dem))
disp(tau)
disp(del)