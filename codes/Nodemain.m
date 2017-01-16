%Delay Node Anaylsis Main Application code
%We have AX=Y, find the best solution of X based on poission's distribtuion
%COnsidering a toy example here
clc;
display('==============================================================');
display('Program to find MLE Node Delays Numerically using Poision distribution');
display('==============================================================');
A = [1 1 0 1
    1 0 1 0];
Y = [10;
    7
    ];
display(A);display(Y);
display('Processing...');
EpVal = power(10, -4);
X = NodeDelay(A, Y, EpVal);
if any(X)
    display(X);
else
    display('No Valid Solutions Found');
end    