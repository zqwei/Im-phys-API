function explVar = dpcaEV(Xfull, W, varargin)

% explVar = dpca_explainedVariance(X, W, V) computes various measures of
% explained variance and returns them in a structure explVar. X is the data
% matrix, W is the decoder matrix, V is the encoder matrix. Returned values:
%
%  * explVar.totalVar             - total variance
%  * explVar.totalMarginalizedVar - total variance in each marginalization
%  * explVar.dPCAcomponentVar     - variance of each component (%)
%  * explVar.dPCAmargVar          - variance of each component in each marginalization (%)
%  * explVar.PCAcomponentVar      - variance of each component (%)
%  * explVar.PCAmargVar           - variance of each component in each marginalization (%)


% default input parameters
options = struct('combinedParams', []);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end

% centering
X = Xfull(:,:);
Xfull = bsxfun(@minus, Xfull, mean(X,2));
X = bsxfun(@minus, X, mean(X,2));

% marginalizing
Xmargs = dpca_marginalize(Xfull, 'combinedParams', options.combinedParams, 'ifFlat', 'yes');

% total variance
explVar.totalVar = sum(sum(X.^2));

% total marginalized variance
for i=1:length(Xmargs)
    explVar.totalMarginalizedVar(i) = sum(Xmargs{i}(:).^2);
end

% variance of each component
for i=1:length(Xmargs)
    explVar.dPCAmargVar(i,:) = sum((W' * Xmargs{i}).^2, 2)' / explVar.totalVar * 100;
end
explVar.dPCAcomponentVar = sum(explVar.dPCAmargVar);

% PCA explained variance
[~, ~, Wpca] = svd(X');

for i=1:length(Xmargs)
    explVar.PCAmargVar(i,:) = sum((Wpca' * Xmargs{i}).^2, 2)' / explVar.totalVar * 100;
end
explVar.PCAcomponentVar = sum(explVar.PCAmargVar);


% % PCA explained variance
% [~,S,Wpca] = svd(X');
% S = diag(S);
% S = S(1:size(W,2));
% explVar.cumulativePCA = cumsum(S.^2'/ explVar.totalVar * 100);

% % dPCA explained variance
% Z = W'*X;
% for i=1:size(W,2)
%     explVar.cumulativeDPCA(i) = 100 - sum(sum((X - V(:,1:i)*Z(1:i,:)).^2)) / explVar.totalVar * 100;    
% end
