
fig1 = figure(1);

tvec = [];
runtimes = linspace(0.01,50,100);
for runtime = runtimes
    tic
    [s,x,o,t] = HRAP_Function('TestConfigg',150,runtime,0.002,false,  0.098,0.014,0.762,0.8,101325,  0.013,0.98,293.15,  0.03,0.056,0.95,0.95,  1,0.01,0.36);
    % clear('fig1')
    % plot(o.t,o.F_thr)
    % xlim([0 15])
    % ylim([0,3000])
    % drawnow
    tvec = [tvec,toc];
end

plot(runtimes,tvec);