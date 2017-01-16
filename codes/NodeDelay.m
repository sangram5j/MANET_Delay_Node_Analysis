

function [X] = NodeDelay(A, Y, EpVal)
%Find all the X[] Matix solutions for AX=Y
%-->Get all the 
exitflag = 1;
ROW_IND = 1; COL_IND = 2;
No_unk = size(A, COL_IND);
Kn_EqCnt = size(A, ROW_IND);
RemEqCnt = abs(No_unk - Kn_EqCnt);

if RemEqCnt > 0
    MMax = zeros(1, Kn_EqCnt);
    %Sort the Matix A row wise for Maximal Solutions into A_Tmp
    for iLoop = 1:Kn_EqCnt
        ARowSum = sum(A(iLoop, :));
        if ARowSum > 0
            MMax(iLoop) = No_unk * Y(iLoop) / ARowSum;
        else
            error('A Matrix is singular');
            return;
        end    
    end

    MMaxVal = max(max(Y), ceil(max(MMax)));
    MMinVal = 1; %Currently setting it to worst case scenario


    MRemArr = ones(RemEqCnt, 1);MRemArr(1) = 0;
    %MRemArr = 0:RemEqCnt - 1;
    IncRemArr = zeros(RemEqCnt, 1);IncRemArr(1) = MMinVal;
    BaseVal = MMaxVal;%power(2, No_unk) - 1;
    %XiCnt = funcnPr(power(2, No_unk) - 1, RemEqCnt);
    TotSols = 0;
    LastOne = true;
    Fin_SolArr = zeros(1, No_unk);
    
    %Generate (Xi, Yi) Solutions for the Node topology
    while(LastOne)
        %Populate the remaining Yi values to make the A,Y matrices as square matrix 
        [Y_Tmp, MRemArr, LastOne] = GenYiValues(Y, MRemArr, IncRemArr, MMaxVal + 1);

        %Populate the remaining Xi values to make the A,Y matrices as square matrix
        XiIndexArr = zeros(1, RemEqCnt);EndSolv = false;
        while(LastOne && EndSolv == false) 
            [A_Tmp, XiIndexArr, EndSolv] = GenXiValues(A, XiIndexArr, BaseVal, No_unk);
            %Using Cramer's Rule to find the values of Xi
            [X_Tmp, ValidFlg] = CramsSols(A_Tmp, Y_Tmp, EpVal);
            if ValidFlg && isempty(find(X_Tmp < 0, 1)) && any(X_Tmp)
                if DistinctSol(Fin_SolArr, X_Tmp', EpVal)
                    display(X_Tmp');
                    TotSols = TotSols + 1;
                    Fin_SolArr(TotSols, :) = X_Tmp; 
                end    
            end
        end
    end

else
    [X, DenVal] = CramsSols(A, Y);
    if (DenVal > 0) && size(find(X < 0), COL_IND) > 0
        TotSols = 1;
        Fin_SolArr = X;
    end
end

%Calculate the expectance value
if any(Fin_SolArr) > 0
    Fin_SolArr = FilterDeciSols(Fin_SolArr, EpVal);
    
    display(sprintf('\nPossible solution set for {X1, X2...Xn} are as follows :\nWaiting...'));
    display(Fin_SolArr);
    %Fin_SolArr = [100 400 600; 300 500 700];
    TotSols = size(Fin_SolArr, ROW_IND);
    thetaArr = ones(No_unk, 1);
    EpArr = ones(No_unk, 1) * EpVal;
    while (exitflag)
        ProbArr = ones(TotSols, 1);
        ExpVal =  zeros(No_unk, 1);
        ProbSum = 0;
        for iLoop = 1 : TotSols
            for jLoop = 1 : No_unk
                ProbArr(iLoop) = ProbArr(iLoop) * ProbPart1(thetaArr(jLoop), Fin_SolArr(iLoop, jLoop));
            end
            ProbSum = ProbSum + ProbArr(iLoop);
        end

        for iLoop = 1 : TotSols
            ProbArr(iLoop) = ProbArr(iLoop) / ProbSum;
        end    

        for iLoop = 1 : TotSols
            ExpVal = ExpVal + Fin_SolArr(iLoop, :)' * ProbArr(iLoop);
        end
        %M step: ?istr+1 ? argmaxQ(?, ?i)
        %for iCurr = 1:TotSols
            %ExpVal(iCurr) = max(ExpVal(iCurr), thetaArr(iCurr));
        %end

        %Taking the consecutive iteration diffrences and comparin to Ep Arr for
        %exit condition.
        TmpArr = find(abs(thetaArr - ExpVal) <= EpArr);
        %display(ExpVal);
        if(size(TmpArr, ROW_IND) ==  size(EpArr, ROW_IND))
            exitflag = 0;
        end
        thetaArr = ExpVal;
    end
    X = thetaArr;
else
    X = zeros(No_unk, 1);
end


end

function [SArr, NextCarr] =   BaseAdd(A, B, CarrVal, BaseVal)
ROW_IND = 1; COL_IND = 2;
SArr = zeros(1, size(A, COL_IND));
NextCarr = CarrVal;
A = floor(A); B = floor(B); BaseVal = floor(BaseVal);
if BaseVal >= 1
    if size(A) == size(B)
        for iLoop = 1:size(A, COL_IND)
            SArr(iLoop) = mod(A(iLoop) +  B(iLoop) + NextCarr, BaseVal);
            NextCarr = floor((A(iLoop) +  B(iLoop) + NextCarr) /BaseVal);
        end    
    else        
        error('Numeric Addition : Dimension Error');
    end    
else
    error('Numeric Addition : Base is invalid');
end
end

function [PartVal] = ProbPart1(theta, xVal)
    PartVal = power(theta, xVal) / factorial(xVal);
end

function [XArr] = CopyXSols(SolArr, LocInd_Arr, LocCount, No_unk)
    XArr = zeros(No_unk, 1);
    for iLoop = 1:LocCount
        XArr(iLoop) = SolArr(LocInd_Arr(iLoop));
    end    
end

function [XArr, ValidFlg] = CramsSols(AArr, YArr, EpVal)
    ROW_IND = 1; COL_IND = 2; ValidFlg = true; 
    XArr = zeros(size(AArr, ROW_IND), 1);
    
    DenVal = det(AArr);
    if(not(DenVal <= EpVal))
        if size(AArr, COL_IND) == size(YArr, ROW_IND)
            for iLoop = 1:size(AArr, ROW_IND)
                Ytmp = AArr;
                Ytmp(:, iLoop) = YArr;
                XArr(iLoop) = det(Ytmp) / DenVal;
                if XArr(iLoop) < 0
                    ValidFlg = false;
                    return;
                end    
            end    
        else
            ValidFlg = false;
            error('Cramers Solution function: Dimension Error');
        end
    end    
end

function [A_Tmp, XiIndexArr, EndSolv]= GenXiValues(A, XiIndexArr, BaseVal, No_unk)
    ROW_IND = 1; COL_IND = 2;    
    IncIndexArr = zeros(1, size(XiIndexArr, COL_IND)); IncIndexArr(1) = 1;
    NextCarr = 0;EndSolv = true;
    while(NextCarr == 0)
        [XiIndexArr, NextCarr] = BaseAdd(XiIndexArr, IncIndexArr, 0, BaseVal);
        if DistinctCheck(XiIndexArr) && NextCarr == 0
            SolArr = GenBinArr(XiIndexArr, No_unk, 2);
            A_Tmp = cat(1, A, SolArr);
            EndSolv = false;
            return;
        elseif NextCarr == 1
            A_Tmp = zeros(No_unk, size(A, COL_IND));
            return;
        end   
    end    
end
 
function [Y_Tmp, MRemArr, CarryRst]= GenYiValues(Y, MRemArr, IncRemArr, MMaxVal)
    NextCarr = 0;
    while(NextCarr == 0)
        [MRemArr, NextCarr] = BaseAdd(MRemArr', IncRemArr', 0, MMaxVal);
        if find(MRemArr > 0)
            break;
        end    
    end
    MRemArr = MRemArr';
    Y_Tmp = cat(1, Y, MRemArr);
    if NextCarr > 0
        CarryRst = false;
    else
        CarryRst = true;
    end    
end
    


function [Val] = funcnPr(nVal, rVal)
    [nr] = floor([nVal rVal]);
    Val = factorial(nr(1))/factorial(nr(1) - nr(2));
end

function [Val] = DistinctCheck(Arr)
    ROW_IND = 1; COL_IND = 2;
    noElem = size(Arr, COL_IND); Val = 1;
    if find(Arr == 0)
       Val = 0; return;
    end    
    for iLoop = 1:noElem - 1
        for jLoop = iLoop + 1:noElem
            if find(Arr(iLoop) == Arr(jLoop:noElem)) > 0
                Val = 0;
                return;
            end    
        end    
    end
end

function [ChkFlg] = DistinctSol(Ain, Bchk, EpVal)
    ROW_IND = 1; COL_IND = 2;
    ChkFlg = true;
    for iLoop = 1:size(Ain, ROW_IND)
        if size(find(abs(Ain(iLoop, :) - Bchk) < EpVal), COL_IND) == size(Bchk, COL_IND)
            ChkFlg = false;
            return;
        end    
    end
end
function [ValArr] = GenBinArr(Arr, veclen, baseval)
    ROW_IND = 1; COL_IND = 2;
    noElem = size(Arr, COL_IND); Val = 1;
    ValArr = zeros(noElem, veclen);
    for iLoop = 1:noElem
        Tmp = Arr(iLoop);
        for jLoop = 1:veclen
            ValArr(iLoop, jLoop) = mod(Tmp, baseval);Tmp = floor(Tmp / baseval);                
        end    
    end
end

%Filter the decimal solutions
function [NewArr] = FilterDeciSols(Arr, EpVal)
    ROW_IND = 1; COL_IND = 2; jLoop = 1;
    TmpArr = [];
    for iLoop = 1:size(Arr, ROW_IND)
        Tmp = (ceil(Arr(iLoop, :)) - Arr(iLoop, :)) > EpVal; 
        if (any(Tmp) == 0)
            TmpArr(jLoop,:) = Arr(iLoop, :);
            jLoop = jLoop + 1;
        end
    end
    NewArr = double(uint64(TmpArr));
end