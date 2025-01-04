var c l k z;
varexo e;

model;
    z = .95*z(-1)+.007*e;
    # lam = 1/c;
    # lamplus = 1/c(+1);
    lam = .97*lamplus*(.33*exp(z(+1))*(k/l(+1))^(-.67) + .9);
    l = lam*.67*exp(z)*(k/l)^.33;
    k = exp(z)*k(-1)^.33*l^.67 + .9*k(-1) - c;
end;

initval;
    z = 0;
    k = 3.7612;
    l = 0.9465;
    c = 1.1161;
end;

shocks;
    var e = 1;
end;

steady;
stoch_simul(order=3, pruning);