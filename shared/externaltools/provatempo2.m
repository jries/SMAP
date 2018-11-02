clc;clear;


Lmax    = 2048;
Lmin    = 3;

FFTrv  = zeros(Lmax,1);
FFTiv  = zeros(Lmax,1);
IFFTiv = zeros(Lmax,1);

% If you want to modify/add new data simply comment/decomment the
% corresponding lines.
% load fftexecutiontimes

tempomin   = 3;
incremento = 10;
tsoglia    = 1;
volte0     = 1000000;
aggiunta   = 5000;
for L=Lmin:Lmax
    disp('---');
    disp('Volte0');
    disp(volte0);
    disp('---');
    a = rand(L,1);
    b = rand(L,1)+i*rand(L,1);
    disp('**************************************************************');
    %----------------------------------------------------------------------
    volte = volte0;
    tempo = 0;
    while tempo<tempomin
        disp('C');
        tic;
        for ii=1:volte
            fft(a);
        end
        tempo = toc;
        volteeffettive = volte;
        % volte updating
        if tempo<tsoglia
            volte = volte*incremento;
        else
            volte = round(volte*tempomin/tempo)+aggiunta;
        end
    end
    tmedio               = tempo/volteeffettive;
    FFTrv(L) = tmedio;
    disp('Length FFT real');
    disp(L);
    disp('Avg time');
    disp(tmedio);
    disp('Total time');
    disp(tempo);
    %--------------
    voltesalva = round(volteeffettive*tempomin/tempo);
    %----------------------------------------------------------------------
    volte = volte0;
    tempo = 0;
    while tempo<tempomin
        disp('C');
        tic;
        for ii=1:volte
            fft(b);
        end
        tempo = toc;
        volteeffettive = volte;
        % volte updating
        if tempo<tsoglia
            volte = volte*incremento;
        else
            volte = round(volte*tempomin/tempo)+aggiunta;
        end
    end
    tmedio                   = tempo/volteeffettive;
    FFTiv(L) = tmedio;
    disp('Length FFT complex');
    disp(L);
    disp('Avg time');
    disp(tmedio);
    disp('Total time');
    disp(tempo);
    %----------------------------------------------------------------------
    volte = volte0;
    tempo = 0;
    while tempo<tempomin
        disp('C');
        tic;
        for ii=1:volte
            ifft(b);
        end
        tempo = toc;
        volteeffettive = volte;
        % volte updating
        if tempo<tsoglia
            volte = volte*incremento;
        else
            volte = round(volte*tempomin/tempo)+aggiunta;
        end
    end
    tmedio                    = tempo/volteeffettive;
    IFFTiv(L) = tmedio;
    disp('Length IFFT complex');
    disp(L);
    disp('Avg time');
    disp(tmedio);
    disp('Total time');
    disp(tempo);
    %----------------------------------------------------------------------
    if (mod(L,100)==0) ||(L==Lmax)
        save fftexecutiontimes FFTrv FFTiv IFFTiv
    end
    volte0 = round(voltesalva*1.25);
    if volte0<2000
        volte0 = 2000;
    end
end