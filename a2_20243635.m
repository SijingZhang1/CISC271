function [rmsvars lowndx rmstrain rmstest] = a2_20243635
% [RMSVARS LOWNDX RMSTRAIN RMSTEST]=A3 finds the RMS errors of
% linear regression of the data in the file "GOODS.CSV" by treating
% each column as a vector of dependent observations, using the other
% columns of the data as observations of independent varaibles. The
% individual RMS errors are returned in RMSVARS and the index of the
% smallest RMS error is returned in LOWNDX. For the variable that is
% best explained by the other variables, a 5-fold cross validation is
% computed. The RMS errors for the training of each fold are returned
% in RMSTEST and the RMS errors for the testing of each fold are
% returned in RMSTEST.
%
% INPUTS:
%         none
% OUTPUTS:
%         RMSVARS  - 1xN array of RMS errors of linear regression
%         LOWNDX   - integer scalar, index into RMSVALS
%         RMSTRAIN - 1x5 array of RMS errors for 5-fold training
%         RMSTEST  - 1x5 array of RMS errors for 5-fold testing
% smallest RMS error is returned in LOWNDX. For the variable that is
% best explained by the other variables, a 5-fold cross validation is
% computed. The RMS errors for the training of each fold are returned
% in RMSTEST and the RMS errors for the testing of each fold are
% returned in RMSTEST.
%
% INPUTS:
%         none
% OUTPUTS:
%         RMSVARS  - 1xN array of RMS errors of linear regression
%         LOWNDX   - integer scalar, index into RMSVALS
%         RMSTRAIN - 1x5 array of RMS errors for 5-fold training
%         RMSTEST  - 1x5 array of RMS errors for 5-fold testing

    filename = 'A3.csv';
    [rmsvars lowndx] = a2q1(filename);
    [rmstrain rmstest] = a2q2(filename, lowndx)

end
function [rmsvars, lowndx] = a2q1(filename)
% [RMSVARS, LOWNDX]=A2Q1(FILENAME) finds the RMS errors of
% linear regression of the data in the file FILENAME by treating
% each column as a vector of dependent observations, using the other
% columns of the data as observations of independent varaibles. The
% individual RMS errors are returned in RMSVARS and the index of the
% smallest RMS error is returned in LOWNDX. 
%
% INPUTS
%         FILENAME - character string, name of file to be processed;
%                    assume that the first row describes the data variables
% OUTPUTS:
%         RMSVARS  - 1xN array of RMS errors of linear regression
%         LOWNDX   - integer scalar, index into RMSVALS

    % Read the test data from CSV file; find the size of the data
    % Start reading at row 2 and column 2
    Amat = csvread(filename, 1,1);
    % Finding data size
    [rows, cols] = size(Amat);

    % Standardize the data
    standardizedData = zscore(Amat);

    % Pre-allocate array for RMS errors
    rmsvars = zeros(1, cols);

    for i = 1:cols
        % Extract the current column as the dependent variable
        dependentVariable = standardizedData(:,i);
        % Extract all other columns as independent variables
        independentVariables = standardizedData(:,1:cols ~= i);

        % Calculate weight vector for the standardized data
        weights = independentVariables\dependentVariable;
        % Calculate predicted dependent variable
        predictedDependentVariable = independentVariables*weights;
        % Calculate RMS error
        rmsvars(i) = rms(dependentVariable - predictedDependentVariable);
    end
    
    % Find the index of the smallest RMS error
    [val, lowndx] = min(rmsvars);
    rmsvars
    % Find the regression on your choice of standardized
    % or unstandardized variables
    astand = standardizedData(:, [1:lowndx - 1, lowndx + 1:end]);
    cstand = standardizedData(:, lowndx);
    wstand = astand \ cstand;
    cpred = astand * wstand;
    rmsstand = rms(cpred - cstand);
    disp("rmsstand")
    disp(rmsstand);
end





function [rmstrain, rmstest] = a2q2(filename,lowndx)
% [RMSTRAIN RMSTEST]=A3Q2(LOWNDX) finds the RMS errors of 5-fold
% cross-validation for the variable LOWNDX of the data in the file
% FILENAME. The RMS errors for the training of each fold are returned
% in RMSTEST and the RMS errors for the testing of each fold are
% returned in RMSTEST.
%
% INPUTS:
%         FILENAME - character string, name of file to be processed;
%                    assume that the first row describes the data variables
%         LOWNDX   - integer scalar, index into the data
% OUTPUTS:
%         RMSTRAIN - 1x5 array of RMS errors for 5-fold training
%         RMSTEST  - 1x5 array of RMS errors for 5-fold testing

    % Read the test data from a CSV file; find the size of the data
    % Create Xmat and yvec from the data and the input parameter,
    % accounting for no standardization of data
    % %
    % % STUDENT CODE GOES HERE: REMOVE THIS COMMENT
    % % THEN ASSIGN THE VARIABLES FROM THE DATASET
    % %

    % Compute the RMS errors of 5-fold cross-validation
    % %
    % % STUDENT CODE GOES HERE: REMOVE THE NEXT 2 LINES AND THIS COMMENT
    % % THEN PERFORM THE COMPUTATIONS
    % %
% Read the test data from a CSV file; find the size of the data
dataraw = csvread(filename,1,1);
[rows, cols] = size(dataraw);

% Standardize the data
Amat = zscore(dataraw);

% Create Xmat and yvec from the data and the input parameter,
% accounting for no standardization of data
Xmat = dataraw(:, setdiff(1:cols, lowndx));
yvec = dataraw(:, lowndx);

% [RMSTRAIN,RMSTEST]=MYKFOLD(XMAT,yvec,K) performs a k-fold validation
% of the least-squares linear fit of yvec to XMAT. If K is omitted,
% the default is 5.
%
% INPUTS:
%         XMAT     - MxN data vector
%         yvec     - Mx1 data vector
%         K        - positive integer, number of folds to use
% OUTPUTS:
%         RMSTRAIN - 1xK vector of RMS error of the training fits
%         RMSTEST  - 1xK vector of RMS error of the testing fits
    Xmat = Amat(:, [1:lowndx-1, lowndx+1:cols]);
    % Problem size
    yvec = Amat(:, lowndx);

    function [rmstrain,rmstest]=mykfold(Xmat, yvec, k_in)
% [RMSTRAIN,RMSTEST]=MYKFOLD(XMAT,yvec,K) performs a k-fold validation
% of the least-squares linear fit of yvec to XMAT. If K is omitted,
% the default is 5.
%
% INPUTS:
%         XMAT     - MxN data vector
%         yvec     - Mx1 data vector
%         K        - positive integer, number of folds to use
% OUTPUTS:
%         RMSTRAIN - 1xK vector of RMS error of the training fits
%         RMSTEST  - 1xK vector of RMS error of the testing fits

    % Problem size
    M = size(Xmat, 1);

    % Set the number of folds
    if nargin >= 3 & ~isempty(k_in)
        k = max(min(round(k_in), M-1), 2)
    else
        k = 3;
    end

    % Initialize the return variables
    rmstrain = zeros(1, k);
    rmstest  = zeros(1, k);

    % Determine the number of rows per experiment
    fold_size = floor(M/k)
    test_indices = zeros(k, fold_size);
    train_indices = zeros(k, M-fold_size);

    rng('default') %makes sure get same resuts each time (consistent seed)
    randidx = randperm(linspace(1,size(Xmat,1),1)); %random permutation of rows
    
    Xmatperm = Xmat(randidx, :); %rearrange order of rows in A
    yvecperm = yvec(randidx, :); %rearrange order of rows in C

    % Process each fold
    for ix = 1:k
        % Define start and end indices for the test set
        test_start = (ix-1)*fold_size + 1;
        test_end = test_start + fold_size - 1;
        % Create an array of indices for the test set and training set
        test_indices(ix,:) = test_start:test_end;
        train_indices(ix,:) = [1:test_start-1, test_end+1:M];

        % Split the data into training and test sets
        xmat_train = Xmatperm(train_indices(ix,:),:);
        yvec_train = yvecperm(train_indices(ix,:));

        % Compute "wvec" for the training data
        wvec = (xmat_train'*xmat_train)\(xmat_train'*yvec_train);

        % Define testing set data
        xmat_test = Xmatperm(test_indices(ix,:),:);
        yvec_test = yvecperm(test_indices(ix,:));
        
        % Compute the RMS error for each fold
        rmstrain(ix) = rms(xmat_train*wvec - yvec_train);
        rmstest(ix) = rms(xmat_test*wvec - yvec_test);
    end
end

    [rmstrain,rmstest] = mykfold(Xmat, yvec, 3);

end
